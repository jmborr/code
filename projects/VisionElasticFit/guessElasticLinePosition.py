'''
This script will try to find the elastic line and the adjacent peak
for a certain number of pixels.
Edit the below variables, then run the script
'''

import sys
sys.path.append('/home/jmborr/devel/mantidproject/mantid/Code/debug/bin')
from mantid.simpleapi import *
import argparse
import numpy
from pdb import set_trace as trace

tofmax=36000
DT=100
dt=4.0  # bin spacing (in microseconds)

############################

from math import exp, erfc, sqrt
from copy import deepcopy

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
  #print 'search assuming strong peak signal'
  #for fwhm in range (40,200,5):
  #  peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=True,BackgroundType='Linear')
  #  line = compareLines(line,peaks,dx*fwhm,tofRange)

  #print 'search assuming highbackground'
  #for fwhm in range (40,200,5):
  #  peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=False,BackgroundType='Linear')
  #  line = compareLines(line,peaks,dx*fwhm,tofRange)  

  #print 'scan assumming position of the peak'
  tofList=range(int(tofRange[0]),int(tofRange[1]),int(dt))
  for fwhm in range (40,200,20):
    peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=True,BackgroundType='Linear',PeakPositions=tofList)
    line = compareLines(line,peaks,dx*fwhm,tofRange)
    peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=False,BackgroundType='Linear',PeakPositions=tofList)
    line = compareLines(line,peaks,dx*fwhm,tofRange)
  #print '\nFindPeaks returns=',line
  return line


def findElasticLinePosition(spectrum,tofRange=(0,5000)):
  """Fit the elastic line with a Gaussian, and return the center position  """
  sp = CloneWorkspace(spectrum,OutputWorkspace=spectrum.getName()+'_cloned')
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

def FindRoughGuess(spectrum):
  """Find a rough estimate of the elastic line"""
  spRebinned = Rebin(InputWorkspace=spectrum,OutputWorkspace='FindRoughGuess.rebin', params=(0,dt,tofmax), PreserveEvents=0)
  results=Max(InputWorkspace='FindRoughGuess.rebin',OutputWorkspace='FindRoughGuess.results')
  maxX = results.readX(0)[0]
  maxY = results.readY(0)[0]
  DeleteWorkspace('FindRoughGuess.rebin')
  DeleteWorkspace('FindRoughGuess.results')
  return maxX, maxY
  
if __name__=='__main__':
  import re
  pattern=re.compile('tube(\d+)') # Regular expression to extract the tube ID
  parser = argparse.ArgumentParser(description='guess elastic line for each pixel')
  parser.add_argument('--nexusFile',help='event nexus file')
  parser.add_argument('--bankNumber',help='number of bank to load')
  parser.add_argument('--outFile',help='output file')
  args = parser.parse_args()
  pFile = open(args.outFile,'w')
  pFile.write('#1.workspace_index 2.detectorID 3.peakCenter 4.peakHeight 5.peakWidth 6.Nevents 7.bank 8.tube\n')
  #for ibank in (1,2,4,5,7,8,10,11,13,14,16,17,19,20):
  ibank=0
  try:
    ibank=int(args.bankNumber)
    ws=Load(Filename=args.nexusFile,OutputWorkspace='events',BankName=r'bank%d'%ibank)
  except:
    ws=Load(Filename=args.nexusFile,OutputWorkspace='events')
  min_nevents = 15000 # assessed from script collectNevents.py
  min_height = 100 # required peak height
  lines='' # buffer per bank analyzed
  for ix in range( ws.getNumberHistograms() ):
    sp = ExtractSingleSpectrum(InputWorkspace='events',WorkspaceIndex=ix,OutputWorkspace='events_single')
    n_events=sp.readY(0)[0]
    if n_events < min_nevents:
      DeleteWorkspace('events_single') # require spectrum contains a minimum number of events
      continue
    rough_guess,height=FindRoughGuess('events_single') #rough estimate of the elastic line position
    if height < min_height:
      DeleteWorkspace('events_single') # require elastic line to contain minimum intensity
      continue
    det_id=sp.getDetector(0).getID() # pixel ID
    tube_id=int( pattern.search( sp.getDetector(0).getFullName() ).groups()[0] )
    spRebinned = Rebin( sp, OutputWorkspace='events_histo', params=(rough_guess-5*DT,dt,rough_guess+5*DT), PreserveEvents=0)
    peak_center,height,width=findElasticLinePosition(spRebinned,tofRange=(rough_guess-DT,rough_guess+DT))
    line='%04d %05d %4.0f %6.0f %6.0f %7d %2d %2d\n'%(ix, det_id,peak_center,height,width,n_events,ibank,tube_id)
    sys.stdout.write(line)
    lines+=line
    for suffix in ('_single','_histo','_histo_cloned','_histo_peaks','_histo_NormalisedCovarianceMatrix','_histo_Parameters','_histo_Workspace'):
      DeleteWorkspace('events'+suffix)
  pFile.write(lines)
  DeleteWorkspace('events')
  pFile.close()
