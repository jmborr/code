'''
This script will try to find the elastic line and the adjacent peak
for a certain number of pixels.
Edit the below variables, then run the script
'''

import sys
sys.path.append('/home/jmborr/devel/mantidproject/mantid/Code/debug/bin')
#from mantid.simpleapi import *
from mantid.simpleapi import FindPeaks, CloneWorkspace, Fit, ExtractSingleSpectrum, Rebin, DeleteWorkspace, Load, Max
from pdb import set_trace as trace
from copy import deepcopy

#DT=100
dt=4.0  # bin spacing (in microseconds)
minHeight=200 # minimum peak height
minIntensity=100 # minimum intensity for the second pass

def FindHighestPeak(sp,PeakFunction='Gaussian', tofRange=(2000,5000), dx=1.0):
  """Use FindPeaks algorithm with a variable FWHM to find all possible peaks
  (2000,5000) is the expected TOF range where the elastic line can be found
  dx is the bin width (in microseconds)
  This thing is an overkill and takes a lot of time, but it should find the peak
  """

  def compareLines(line,peaks,targetWidth,tofRange):
    """ Parse the results from FindPeaks in order to retrieve the best peak candidate,
    then compare to current best peak """
    nRows = peaks.rowCount()
    for iRow in range(0,nRows):
      if not line or line['height'] < peaks.row( iRow )['height']:
        width = peaks.row( iRow )['width']  #is the peak within width range ?
        if 0.5*targetWidth<width and width<2*targetWidth:
          centre = peaks.row( iRow )['centre']  #is the peak within range ?
          if tofRange[0]<centre and centre<tofRange[1]:
            line = deepcopy( peaks.row( iRow ) )
    return line
  line = None

  #print 'scan assumming position of the peak'
  tofList=range(int(tofRange[0]),int(tofRange[1]),int(dt))
  for fwhm in range (40,200,20):
    print 'fwhm=',fwhm,'HighBackground=True'
    peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=True,BackgroundType='Linear',PeakPositions=tofList)
    line = compareLines(line,peaks,dx*fwhm,tofRange)
    print 'fwhm=',fwhm,'HighBackground=False'
    peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=False,BackgroundType='Linear',PeakPositions=tofList)
    line = compareLines(line,peaks,dx*fwhm,tofRange)
  #print '\nFindPeaks returns=',line
  return line

def findElasticLinePosition(spectrum,tofRange=(0,5000)):
  """Fit the elastic line with a Gaussian, and return the center position  """
  peak = FindHighestPeak(spectrum,tofRange=tofRange,dx=dt) #find the strongest peak in the TOF range
  #Fit the line to a Gaussian. Initial guess from output from FindHighestPeak
  try:
    centre = peak['centre']; sigma = peak['width']/2.35 #From FWHM to standard deviation
  except TypeError:
    print 'could not find peak for spectrum index ',spectrum.getDetector(0).getID()
    raise
  startX = centre-sigma;  endX = centre+sigma #set the fitting boundaries
  funcStr = 'name=Gaussian, Height=%f, PeakCentre=%f, Sigma=%f'%(peak['height'], centre, sigma)
  gfit = Fit(funcStr, spectrum.getName(), StartX=startX, EndX=endX, CreateOutput='1')
  #Retrieve the optimized fitting parameters
  fittedParams = { 'Height':0, 'PeakCentre':0, 'Sigma':0}
  for iRow in range( gfit[3].rowCount() - 1):
    row = gfit[3].row( iRow )
    fittedParams[ row['Name'] ] = row['Value']
  return fittedParams['PeakCentre'], fittedParams['Height'],fittedParams['Sigma']

def readIn(datafile):
  """ Read in the output data file from the first pass of guessElasticLine.
  The output structure is a dictionary whose keys are the tube ID's and
  whose values are the list of detectors in each tube """
  fps={}
  for line in open(datafile,'r').readlines():
    if line[0]=='#': continue
    items=line.split()
    entry={'w_index':int(items[0]),
           'det_id':int(items[1]),
           'peak_center':float(items[2]),
           'peak_height':float(items[3]),
           'peak_width':float(items[4]),
           'nevents':int(items[5]),
           'bank':int(items[6]),
           'tube_id':int(items[7]),
           'line':line}
    if entry['tube_id']  not in fps.keys():
      fps[entry['tube_id']] = [entry,]
    else:
      fps[entry['tube_id']].append(entry)
  return fps

def linearExtrapolate(fp, iX, direction):
  """ Find extrapolated fit """
  a= (fp[iX+1]['peak_center']-fp[iX-1]['peak_center']) / (fp[iX+1]['w_index']-fp[iX-1]['w_index']) #slope
  center = a*(fp[iX+2*direction]['w_index']-fp[iX]['w_index']) + fp[iX]['peak_center'] # estimated peak_center at fp[iX+2*direction]['w_index']
  #a= (fp[iX+1]['peak_height']-fp[iX-1]['peak_height']) / (fp[iX+1]['w_index']-fp[iX-1]['w_index']) #slope
  #height = a*(fp[iX+2*direction]['w_index']-fp[iX]['w_index']) + fp[iX]['peak_height'] # estimated peak_height at fp[iX+2*direction]['w_index']
  #a= (fp[iX+1]['peak_width']-fp[iX-1]['peak_width']) / (fp[iX+1]['w_index']-fp[iX-1]['w_index']) #slope
  #width = a*(fp[iX+2*direction]['w_index']-fp[iX]['w_index']) + fp[iX]['peak_width'] # estimated peak_width at fp[iX+2*direction]['w_index']
  #return [center,height,width]
  return [ center,fp[iX+2*direction]['peak_height'],fp[iX+2*direction]['peak_width'] ]

def isIntense(spectrum, center, neighborhood):
  """ Confirm if there is strong peak in the neighborhood of center """
  results=Max(InputWorkspace=spectrum.getName(),OutputWorkspace='isIntense.results',RangeLower=center-neighborhood,RangeUpper=center+neighborhood)
  maxY = results.readY(0)[0]
  DeleteWorkspace('isIntense.results')
  if maxY > minIntensity:
    return True
  return False

def replacementLine( fp, ix, center, width, neighborhood):
  """ Find the peak within the constraints """
  print 'ix=',ix, 'w_index',fp[ix]['w_index'], 'det_id=',fp[ix]['det_id'], 'center=',center, 'width=',width, 'neighborhood=',neighborhood
  sp = ExtractSingleSpectrum(InputWorkspace='events',WorkspaceIndex=fp[ix]['w_index'],OutputWorkspace='events_single')
  spRebinned = Rebin( sp, OutputWorkspace='events_histo', params=(max(0,center-10*width),dt,center+10*width), PreserveEvents=0)
  #check if we have enough intensity in the region of TOF to investigate
  if not isIntense(spRebinned, center, neighborhood):
    return None
  center,height,width=findElasticLinePosition(spRebinned,tofRange=(center-neighborhood,center+neighborhood))
  fp[ix]['peak_center']=center # Update/correct the peak info, which will serve for the next extrapolation
  fp[ix]['peak_height']=height
  fp[ix]['peak_width']=width
  for suffix in ('_single','_histo', '_histo_peaks','_histo_NormalisedCovarianceMatrix','_histo_Parameters','_histo_Workspace'):
    DeleteWorkspace('events'+suffix)
  line='%04d %05d %4.0f %6.0f %6.0f %7d %2d %2d\n'%(fp[ix]['w_index'],fp[ix]['det_id'],center,height,width,fp[ix]['nevents'],fp[ix]['bank'],fp[ix]['tube_id'])
  return line

if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='guess elastic line for each pixel')
  parser.add_argument('--nexusFile',help='event nexus file')
  parser.add_argument('--bankNumber',help='number of bank to load')
  parser.add_argument('--firstPass',help='data file, the output of guessElasticLinePosition.py')
  parser.add_argument('--outFile',help='output file')
  args = parser.parse_args()
  tol_center = 40.0 # maximum tolerated change in elastic-line-TOF between consecutive detectors
  pFile = open(args.outFile,'w')
  pFile.write('#1.workspace_index 2.detectorID 3.peakCenter 4.peakHeight 5.peakWidth 6.Nevents 7.bank 8.tube\n')
  ibank=0
  try:
    ibank=int(args.bankNumber)
    ws=Load(Filename=args.nexusFile,OutputWorkspace='events',BankName=r'bank%d'%ibank)
  except:
    ws=Load(Filename=args.nexusFile,OutputWorkspace='events')
  fps=readIn(args.firstPass) # read in the output data file from the first pass of guessElasticLine.
  tube_IDs=fps.keys()
  tube_IDs.sort()
  trace()
  for tube_id in tube_IDs:
    print 'tube_id=',tube_id
    fp=fps[tube_id]
    iL=len(fp)/2 # Central index
    lines=fp[iL-1]['line']+fp[iL]['line']+fp[iL+1]['line']
    # process towards smaller detector ID's from the center of the tube
    while iL>1:
      [center,height,width]=linearExtrapolate(fp,iL,-1)
      if abs( (center-fp[iL-2]['peak_center'])/(fp[iL-2]['w_index']-fp[iL]['w_index']) ) < tol_center:
        # The peak is located where it is supposed to be
        lines=fp[iL-2]['line']+lines # prepend to buffer
      else:
        # The peak is located outside the region where it was expected to be. Find peak in the expected region
        line=replacementLine(fp, iL-2, center, width, tol_center)
        if not line:
          iL=1 # not enough intensity in the neighborhood of center, thus stop
        else:
          lines=line+lines # prepend to buffer
      iL-=1 # "descend" along the tube
    # process towards bigger detector ID's from the center of the tube
    iR=len(fp)/2 # Central index
    while iR<len(fp)-2:
      #print iR,len(fp)
      [center,height,width]=linearExtrapolate(fp,iR,+1)
      if abs( (center-fp[iR+2]['peak_center'])/(fp[iR+2]['w_index']-fp[iR]['w_index']) ) < tol_center:
        lines+=fp[iR+2]['line'] # append to buffer
      else:
        line=replacementLine(fp, iR+2, center, width, tol_center)
        if not line:
          iR=len(fp)-2 # not enough intensity in the neighborhood of center, thus stop
        else:
          lines+=line # append to buffer
      iR+=1 # "ascend" along the tube
    pFile.write(lines) # store the processed tube
  DeleteWorkspace('events')
  pFile.close()