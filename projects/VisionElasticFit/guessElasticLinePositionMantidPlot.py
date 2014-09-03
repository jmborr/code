'''
This script will try to find the elastic line and the adjacent peak
for a certain number of pixels.
Edit the below variables, then run the script
'''

nexus_file='/projects/development/VisionElasticFit/test2/VIS_1076_event.nxs'  # file containing the run It must be an EVENT nexus file
tofRange=(2000,5000) #range where we expect the elastic line to show up
dt=5.0  # bin spacing (in microseconds)
Max=5000  # maximum time of flight where we expect to find the elastic line
outParmFile='/tmp/parameters.dat' # file to output summary fitting data
nbank=20  # number of inelastic banks

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


def findElasticLinePosition(spectrum):
  """Fit the elastic line with a Gaussian, and return the center position  """

  sp = CloneWorkspace(spectrum,OutputWorkspace=spectrum.getName()+'_cloned')
  peak = FindHighestPeak(spectrum,tofRange=tofRange,dx=dt) #find the strongest peak in the TOF range

  #Fit the line to a Gaussian. Initial guess from output from FindHighestPeak
  centre = peak['centre']; sigma = peak['width']/2.35 #From FWHM to standard deviation
  startX = centre-1.5*sigma;  endX = centre+1.5*sigma #set the fitting boundaries
  funcStr = 'name=Gaussian, Height=%f, PeakCentre=%f, Sigma=%f'%(peak['height'], centre, sigma)
  gfit = Fit(funcStr, spectrum.getName(), StartX=startX, EndX=endX, CreateOutput='1')

  #Retrieve the optimized fitting parameters
  fittedParams = { 'Height':0, 'PeakCentre':0, 'Sigma':0}
  for iRow in range( gfit[3].rowCount() - 1):
    row = gfit[3].row( iRow )
    fittedParams[ row['Name'] ] = row['Value']
  
  return fittedParams['PeakCentre'], fittedParams['Height']


pFile = open(outParmFile,'w')
pFile.write('#1.detectorID 2.TOF 3.Nevents 4.bank\n')
#for ibank in (1,2,4,5,7,8,10,11,13,14,16,17,19,20):
for ibank in (4,5,7,8,16,17,19,20):
  Load(Filename=nexus_file,OutputWorkspace='events',BankName=r'bank%d'%ibank) # load only inelastic banks
  ws=mantid.getMatrixWorkspace('events')
  min_nevents = 20000 # assessed from script collectNevents.py
  min_heigth = 125 # required peak heigth
  lines='' # buffer per bank analyzed
  for ix in range( ws.getNumberHistograms() ):
    sp = ExtractSingleSpectrum(InputWorkspace='events',WorkspaceIndex=ix,OutputWorkspace='events_single')
    n_events=sp.readY(0)[0]
    if n_events < min_nevents:
      DeleteWorkspace('events_single') # require spectrum contains a minimum number of events
      continue
    det_id=sp.getDetector(0).getID() # pixel ID  
    spRebinned = Rebin( sp, OutputWorkspace='events_histo', params=(0,dt,Max), PreserveEvents=0)
    peak_center,heigth=findElasticLinePosition(spRebinned)
    line='%05d %4.0f %6f.0 %7d %2d\n'%(det_id,peak_center,heigth,n_events,ibank)
    sys.stdout.write(line)
    if heigth > min_heigth:
	lines+=line
    for suffix in ('_single','_histo','_histo_cloned','_histo_peaks','_histo_NormalisedCovarianceMatrix','_histo_Parameters','_histo_Workspace'):
      DeleteWorkspace('events'+suffix)
  pFile.write(lines)
  DeleteWorkspace('events')
pFile.close()
