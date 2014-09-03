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
minIntensity=50 # minimum intensity for the second pass
n_tube=8    # eight tubes per pack
n_pixel_per_tube=128 # 128 tubes per pixel
n_pixel = n_tube * n_pixel_per_tube # total number of pixels (1024)

def FindHighestPeak(sp,PeakFunction='Gaussian', tofRange=(2000,5000), dx=1.0,guess_width=100):
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

  #print 'scan assuming position of the peak'
  tofList=range(int(tofRange[0]),int(tofRange[1]),int(dt))
  fwhm=int(guess_width*0.8)
  max_fwhm=int(guess_width*1.2)
  delta_fwhm=max(1,(max_fwhm-fwhm)/5)
  while fwhm <= max_fwhm:
    print 'fwhm=',fwhm,'HighBackground=True'
    peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=True,BackgroundType='Linear',PeakPositions=tofList)
    line = compareLines(line,peaks,dx*fwhm,tofRange)
    print 'fwhm=',fwhm,'HighBackground=False'
    peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=False,BackgroundType='Linear',PeakPositions=tofList)
    line = compareLines(line,peaks,dx*fwhm,tofRange)
    fwhm += delta_fwhm
  #print '\nFindPeaks returns=',line
  return line

def findElasticLinePosition(spectrum,tofRange=(0,5000),guess_width=100):
  """Fit the elastic line with a Gaussian, and return the center position  """
  peak = FindHighestPeak(spectrum,tofRange=tofRange,dx=dt,guess_width=guess_width) #find the strongest peak in the TOF range
  #Fit the line to a Gaussian. Initial guess from output from FindHighestPeak
  try:
    centre = peak['centre']
    sigma = peak['width']/2.35 #From FWHM to standard deviation
  except TypeError:
    print 'could not find peak for spectrum index ',spectrum.getDetector(0).getID()
    return None
  startX = centre-sigma
  endX = centre+sigma #set the fitting boundaries
  funcStr = 'name=Gaussian, Height=%f, PeakCentre=%f, Sigma=%f'%(peak['height'], centre, sigma)
  gfit = Fit(funcStr, spectrum.getName(), StartX=startX, EndX=endX, CreateOutput='1')
  #Retrieve the optimized fitting parameters
  fittedParams = { 'Height':0, 'PeakCentre':0, 'Sigma':0}
  for iRow in range( gfit[3].rowCount() - 1):
    row = gfit[3].row( iRow )
    fittedParams[ row['Name'] ] = row['Value']
  return fittedParams

def readIn(datafile):
  """ Read in the output data file from the first pass of guessElasticLine.
  The output structure is a list on 'detectors' """
  fp=[None,]*n_pixel
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
    i_pixel=entry['w_index']
    fp[i_pixel]=entry
  return fp

def linearExtrapolate(fp, i_pixel, direction):
  """ Find extrapolated linear fit for the peak center
  Return height and width of center pixel as the extrapolated values"""
  a= (fp[i_pixel+1]['peak_center']-fp[i_pixel-1]['peak_center'])/2 # slope for change of center with workspace index
  center = a*2*direction + fp[i_pixel]['peak_center'] # estimated peak_center at i_pixel+2*direction
  return [ center,fp[i_pixel]['peak_height'],fp[i_pixel]['peak_width'] ]

def isIntense(spectrum, center, neighborhood):
  """ Confirm if there is strong peak in the neighborhood of center """
  results=Max(InputWorkspace=spectrum.getName(),OutputWorkspace='isIntense.results',RangeLower=center-neighborhood,RangeUpper=center+neighborhood)
  maxY = results.readY(0)[0]
  DeleteWorkspace('isIntense.results')
  if maxY > minIntensity:
    return True
  return False

def replacementLine( fp, ix, i_pixel, center, width, neighborhood):
  """ Find the peak within the constraints at pixel ix"""
  print 'ix=',ix, 'center=',center, 'width=',width, 'neighborhood=',neighborhood
  sp = ExtractSingleSpectrum(InputWorkspace='events',WorkspaceIndex=ix,OutputWorkspace='events_single')
  spRebinned = Rebin( sp, OutputWorkspace='events_histo', params=(max(0,center-10*width),dt,center+10*width), PreserveEvents=0)
  if not isIntense(spRebinned, center, neighborhood): #check if we have enough intensity in the region of TOF to investigate
    return None
  fittedParams=findElasticLinePosition(spRebinned,tofRange=(center-neighborhood,center+neighborhood),guess_width=width)
  if not fittedParams: # Check if we found a peak
    return None
  center=fittedParams['PeakCentre']
  height=fittedParams['Height']
  width=fittedParams['Sigma']
  if not fp[ix]:
    fp[ix]={}
    fp[ix]['w_index']=ix
    fp[ix]['det_id']=fp[i_pixel]['det_id']+(ix-i_pixel)
    fp[ix]['nevents']=sp.getNumberEvents()
    fp[ix]['bank']=fp[i_pixel]['bank']
    fp[ix]['tube_id']=fp[i_pixel]['tube_id']
  fp[ix]['peak_center']=center # Update/correct the peak info, which will serve for the next extrapolation
  fp[ix]['peak_height']=height
  fp[ix]['peak_width']=width
  for suffix in ('_single','_histo', '_histo_peaks','_histo_NormalisedCovarianceMatrix','_histo_Parameters','_histo_Workspace'):
    DeleteWorkspace('events'+suffix)
  line='%04d %05d %4.0f %6.0f %6.0f %7d %2d %2d\n'%(fp[ix]['w_index'],fp[ix]['det_id'],center,height,width,fp[ix]['nevents'],fp[ix]['bank'],fp[ix]['tube_id'])
  fp[ix]['line']=line
  return line

if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='guess elastic line for each pixel')
  parser.add_argument('--nexusFile',help='event nexus file')
  parser.add_argument('--bankNumber',help='number of bank to load')
  parser.add_argument('--firstPass',help='data file, the output of guessElasticLinePosition.py')
  parser.add_argument('--centerPixels', help='space-sparated list of workpace indexes serving as tube centers. The do not have to coincide with the workspace indexes corresponding to the center of the tubes. If the script is not to do a second pass on a particular tube, then pass zero for the center pixel')
  parser.add_argument('--outFile',help='output file')
  args = parser.parse_args()
  tol_center = 40.0 # maximum tolerated change in elastic-line-TOF between consecutive detectors
  pFile = open(args.outFile,'w')
  pFile.write('#1.workspace_index 2.detectorID 3.peakCenter 4.peakHeight 5.peakWidth 6.Nevents 7.bank 8.tube\n')
  ibank=int(args.bankNumber)
  ws=Load(Filename=args.nexusFile,OutputWorkspace='events',BankName=r'bank%d'%ibank)
  #trace()
  center_pixel=[ int(x) for x in args.centerPixels.split() ]
  fp=readIn(args.firstPass) # read in the output data file from the first pass of guessElasticLine.
  for i_tube in range(n_tube):
    if not center_pixel[i_tube]: continue # do not process this tube
    begin_pixel=i_tube*n_pixel_per_tube
    i_pixel=center_pixel[i_tube] # workspace index from which we start
    lines=fp[i_pixel-1]['line']+fp[i_pixel]['line']+fp[i_pixel+1]['line']
    # process towards smaller workspace indexes from the center of the tube
    while i_pixel>=begin_pixel+2:
      [center,height,width]=linearExtrapolate(fp,i_pixel,-1)
      if fp[i_pixel-2] and abs((center-fp[i_pixel-2]['peak_center'])/2) < tol_center:
        # The peak exists and is located where it is supposed to be
        lines=fp[i_pixel-2]['line']+lines # prepend to buffer
      else:
        # Non-existant peak, or peak located outside the region where it was expected to be. Find peak in the expected region
        line=replacementLine(fp, i_pixel-2, i_pixel, center, width, tol_center)
        if not line: # not enough intensity in the neighborhood of center, thus stop
          i_pixel=begin_pixel
        else:
          lines=line+lines # prepend to buffer
      i_pixel-=1 # "descend" along the tube
    # process towards bigger detector ID's from the center of the tube
    i_pixel=center_pixel[i_tube] # index of center detector for given tube
    end_pixel=begin_pixel+n_pixel_per_tube-1
    while i_pixel<=end_pixel-2:
      #print iR,len(fp)
      [center,height,width]=linearExtrapolate(fp,i_pixel,+1)
      if fp[i_pixel+2] and abs((center-fp[i_pixel+2]['peak_center'])/2) < tol_center:
        lines+=fp[i_pixel+2]['line'] # append to buffer
      else:
        line=replacementLine(fp, i_pixel+2, i_pixel, center, width, tol_center)
        if not line:
          i_pixel=end_pixel # not enough intensity in the neighborhood of center, thus stop
        else:
          lines+=line # append to buffer
      i_pixel+=1 # "ascend" along the tube
    pFile.write(lines) # store the processed tube
  DeleteWorkspace('events')
  pFile.close()