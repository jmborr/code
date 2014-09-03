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

dt=4.0  # bin spacing (in microseconds)
DT=1500  # maximum time of flight around the expected position of the elastic line
minHeight=200 # minimum peak height

############################

from math import exp, erfc, sqrt
from copy import deepcopy

def FindHighestPeak(sp,PeakFunction='Gaussian', tofRange=(2000,5000), dx=1.0):
  """Use FindPeaks algorithm with a variable FWHM to find all possible peaks
  (2000,5000) is the expected TOF range where the elastic line can be found
  dx is the bin width (in microseconds)
  This thing is an overkill and takes a lot of time, but it should find the peak
  """
  from copy import deepcopy
  def compareLines(line,peaks,targetWidth,tofRange):
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
  for fwhm in range (40,200,5):
    peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=True,BackgroundType='Linear')
    line = compareLines(line,peaks,dx*fwhm,tofRange)

  #print 'search assuming highbackground'
  for fwhm in range (40,200,5):
    peaks = FindPeaks( sp, PeaksList=sp.getName()+'_peaks', FWHM=fwhm, PeakFunction=PeakFunction,HighBackground=False,BackgroundType='Linear')
    line = compareLines(line,peaks,dx*fwhm,tofRange)  
    
  #print 'scan assumming position of the peak'
  tofList = [];  delta = ( tofRange[1] - tofRange[0] ) / 50.;  tof = tofRange[0] + delta; tofList.append( tof )
  while tof < tofRange[1]:
    tofList.append( tof );  tof += delta
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
  startX = centre-1.5*sigma;  endX = centre+1.5*sigma #set the fitting boundaries
  funcStr = 'name=Gaussian, Height=%f, PeakCentre=%f, Sigma=%f'%(peak['height'], centre, sigma)
  gfit = Fit(funcStr, spectrum.getName(), StartX=startX, EndX=endX, CreateOutput='1')
  #Retrieve the optimized fitting parameters
  fittedParams = { 'Height':0, 'PeakCentre':0, 'Sigma':0}
  for iRow in range( gfit[3].rowCount() - 1):
    row = gfit[3].row( iRow )
    fittedParams[ row['Name'] ] = row['Value']
  return fittedParams['PeakCentre'], fittedParams['Height']

def FindAverageTOF(spectrum,dt=100.0,Max=40000.0):
  """Find a rough estimate of the elastic line"""
  spRebinned = Rebin(InputWorkspace=spectrum,OutputWorkspace='FindAverageTOF.rebin', params=(0,dt,Max), PreserveEvents=0)
  L=len(spRebinned.readY(0))
  avTOF=numpy.dot(spRebinned.readX(0)[:L],spRebinned.readY(0))/numpy.sum(spRebinned.readY(0))
  DeleteWorkspace('FindAverageTOF.rebin')
  return avTOF


if __name__=='__main__':
  parser = argparse.ArgumentParser(description='guess elastic line for each pixel')
  parser.add_argument('--nexusFile',help='event nexus file')
  parser.add_argument('--bankNumber',help='number of bank to load')
  parser.add_argument('--outFile',help='output file')
  args = parser.parse_args()
  pFile = open(args.outFile,'w')
  pFile.write('#1.detectorID 2.peakCenter 3.peakHeight 4.Nevents 5.bank\n')
  #for ibank in (1,2,4,5,7,8,10,11,13,14,16,17,19,20):
  try:
    ibank=int(args.bankNumber)
    ws=Load(Filename=args.nexusFile,OutputWorkspace='events',BankName=r'bank%d'%ibank)
  except:
    ws=Load(Filename=args.nexusFile,OutputWorkspace='events')
  min_nevents = 20000 # assessed from script collectNevents.py
  min_heigth = 125 # required peak height
  lines='' # buffer per bank analyzed
  for ix in range( ws.getNumberHistograms() ):
    sp = ExtractSingleSpectrum(InputWorkspace='events',WorkspaceIndex=ix,OutputWorkspace='events_single')
    n_events=sp.readY(0)[0]
    if n_events < min_nevents:
      DeleteWorkspace('events_single') # require spectrum contains a minimum number of events
      continue
    det_id=sp.getDetector(0).getID() # pixel ID
    avTOF=FindAverageTOF('events_single') #rough estimate of the elastic line position
    spRebinned = Rebin( sp, OutputWorkspace='events_histo', params=(max(0,avTOF-2*DT),dt,avTOF+2*DT), PreserveEvents=0)
    peak_center,height=findElasticLinePosition(spRebinned,tofRange=(max(0,avTOF-DT),avTOF+DT))
    if height<minHeight:
      DeleteWorkspace('events_single') # require the elastic line to have enough intense signal
      DeleteWorkspace('events_histo')
      continue
    line='%05d %4.0f %6f.0 %7d %2d\n'%(det_id,peak_center,heigth,n_events,ibank)
    sys.stdout.write(line)
    if heigth > min_heigth:
      lines+=line
    for suffix in ('_single','_histo','_histo_cloned','_histo_peaks','_histo_NormalisedCovarianceMatrix','_histo_Parameters','_histo_Workspace'):
      DeleteWorkspace('events'+suffix)
  pFile.write(lines)
  DeleteWorkspace('events')
  pFile.close()
