'''
This script will try to find the elastic line and the adjacent peak
for a certain number of pixels.
Edit the below variables, then run the script
'''

wsName='VIS_1085_event'  # workspace name of the bank. It must be an EVENT workspace
centerPixel=567   # workspace index of the pixel with maximum intensity
pixelRange=0   # fit for pixelRange pixels to the left and right of centerPixel, in addition to centerPixel
dt=5                 # bin spacing (in microseconds)
Max=16700       # maximum time of flight
outParmFile='/tmp/parameters.dat' #file to output summary fitting data
outFitFilePrefix='/tmp/pixel_XXX.dat' #contains the data and fitted spectrum for a given pixel
############################

from math import exp, erfc, sqrt
from copy import deepcopy


def back2backExp(x,params):
  I=params['I'];  a=params['A'];  b=params['B'];  x0=params['X0'];  s=params['S']
  s2 = s*s
  diff = x-x0
  y1 = 0.0
  try:
    y1 = exp( a/2*(a*s2+2*diff) ) * erfc( (a*s2+diff) / sqrt(2*s2) )
  except:
    y1 = 0.0
  y2 = 0.0
  try:
    y2 = exp( b/2*(b*s2-2*diff) ) * erfc( (b*s2-diff) / sqrt(2*s2) )
  except: 
    y2 = 0.0
    
  return I * ( y1 + y2 )


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
  print '\nFindPeaks returns=',line
  return line


def substractFittedFunction(sp, function, params,background=None, startX=None, endX=None):
  """Set sp.dataY to zero within the function domain
  if we pass a background, substract the background in the (startX, endX) range
  """
  X = sp.dataX(0);  Y = sp.dataY(0)
  for i in range( len(X)-1 ):
    x = (X[i]+X[i+1])/2
    Z = function(x,params)
    if background :
      if startX < X[i] and X[i] < endX:
        Z+= background['A0'] + background['A1'] * X[i]
    Y[i] = max( Y[i] - Z, 0.0)


def fitDataToBackToBackExponential(spectrum,tofRange):
  """Fit spectrum to a BackToBackExponential peak
  spectrum: single spectrum histogram data
  tofRange: (TOFmin,TOFmax) range to find the peak
  Steps for this fitting recipe:
  1. Use algorithm FindPeaks to find the strongest peak
  2. Fit the peak to a Gaussian. Use the previous fit as starting guess.
  3. Fit the peak to a BackToBackExponential. Use the previous fit as starting guess.
  Return the BackToBackExponential optimized parameters
  """
  peak = FindHighestPeak(spectrum,tofRange=tofRange,dx=5.0) #find the strongest peak in the TOF range

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
  print '\nGaussian fit in fitDataToBackToBackExponential', fittedParams
  
  #Fit the line to a BackToBackExponential plus linear background. Initial guess from the previous Gaussian fit
  dX0 = 125
  I = fittedParams['Height'];  X0 = fittedParams['PeakCentre'];  S = fittedParams['Sigma']
  funcStr = 'name=BackToBackExponential,I=%f,A=0.01,B=0.01,X0=%f,S=%f,constraints=(A>0,B>0,X0>%f,X0<%f,S>%f)'%(I,X0,S,X0-dX0,X0+dX0,S/2.0)
  funcStr += ';name=LinearBackground,A0=%f,A1=%f,constraints=(A0>0,A1>0)'%(10.0,0.0)
  gfit = Fit(funcStr, spectrum.getName(), StartX=startX, EndX=endX, CreateOutput='1')

  #Retrieve the optimized fitting parameters
  fittedParams = {'I':0, 'A':0, 'B':0, 'X0':0, 'S':0}
  background = {'A0':0, 'A1':0}
  for iRow in range( gfit[3].rowCount() - 1):
    row = gfit[3].row( iRow )
    (fid,par) = row['Name'].split('.')
    if fid=='f0': fittedParams[ par ] = row['Value']
    if fid == 'f1': background[ par ] = row['Value']
  #Make sure slopes are positive
  for parms in ('A','B'):
    if fittedParams[parms] < 0: fittedParams[parms] = 0.001
    
  return fittedParams,background,(startX,endX)


def fitDataToLinearBackground( spectrum, tofRange ):
  """Fit spectra to a linear background
  """
  funcStr = 'name=LinearBackground'
  gfit = Fit(funcStr, spectrum.getName(), StartX=tofRange[0], EndX=tofRange[1], CreateOutput='1')

  #Retrieve the optimized fitting parameters
  fittedParams = {'A0':0, 'A1':0}
  for iRow in range( gfit[3].rowCount() - 1):
    row = gfit[3].row( iRow )
    fittedParams[ row['Name'] ] = row['Value']

def fitSpectrum(spectrum,fitMethod=None,isToTheRight=True):
  """Fit a single spectrum, histogram data
  isToTheRight: position of the adjacent peak with respect to the elastic line
  """

  if fitMethod == 'two BackToBackExponentials':
    
    sp = CloneWorkspace(spectrum,OutputWorkspace=spectrum.getName()+'_cloned')
    
    #Estimate the elastic peak
    elasticPeak,background,(sE,eE) = fitDataToBackToBackExponential( sp, (2000,5000) )
    print '\nElasticPeak=',elasticPeak
    funcStr = 'name=BackToBackExponential,I=%f,A=%f,B=%f,X0=%f,S=%f,constraints=(A>0,B>0,S>%f)'%(elasticPeak['I'],elasticPeak['A'],elasticPeak['B'],elasticPeak['X0'],elasticPeak['S'],elasticPeak['S']/2.0)
    
    #Estimate the adjacent peak
    substractFittedFunction( sp, back2backExp, elasticPeak,background=background,startX=sE, endX=eE) #substract peak    
    CloneWorkspace(sp,OutputWorkspace=spectrum.getName()+'_removedElastic')
    tofRange = ( elasticPeak['X0'] + elasticPeak['S'],elasticPeak['X0'] + 2000 )
    if not isToTheRight:
      tofRange = ( elasticPeak['X0'] - 2000 , elasticPeak['X0'] - elasticPeak['S'])
    adjacentPeak,background,(sA,eA) = fitDataToBackToBackExponential( sp, tofRange)
    print '\nAdjacentPeak=',adjacentPeak
    
    #Impose constrains in the putative position of the center of the adjacent peak when doing the global fit
    minX0 = adjacentPeak['X0'] - adjacentPeak['S']
    maxX0 = adjacentPeak['X0'] + adjacentPeak['S']
    funcStr += ';name=BackToBackExponential,I=%f,A=%f,B=%f,X0=%f,S=%f,constraints=(A>0,B>0,%f<X0<%f,S>%f)'%(adjacentPeak['I'],adjacentPeak['A'],adjacentPeak['B'],adjacentPeak['X0'],adjacentPeak['S'],minX0,maxX0,adjacentPeak['S']/2.0)

    # Estimate the background
    substractFittedFunction( sp, back2backExp, adjacentPeak ) #substract peak
    CloneWorkspace(sp,OutputWorkspace=spectrum.getName()+'_removedAdjacent')
    startX=min(sE,sA); endX=max(eE,eA)
    startX = max(0,startX-500); endX+=500
    print (startX,endX),sp.getName()
    DeleteWorkspace( sp )
    sp = CloneWorkspace(spectrum,OutputWorkspace=spectrum.getName()+'_cloned')
    substractFittedFunction( sp, back2backExp, elasticPeak )
    substractFittedFunction( sp, back2backExp, adjacentPeak ) 
    background,(sB,eB) = fitDataToLinearBackground( sp, (startX,endX) )
    print '\nBackground=',background
    funcStr += ';name=LinearBackground,A0=%f,A1=%f'%(background['A0'],background['A1'])
    
    #Do the global fit taking into account all the previous estimates
    print '\nFunction: ',funcStr
    print '\nStartX = ',startX,'  EndX = ',endX
    gfit = Fit(funcStr, spectrum.getName(), StartX=startX, EndX=endX, CreateOutput='1')
    
    #Extract fitting parameters
    parameters={'f0':{'I':None,'A':None,'B':None,'X0':None,'S':None},
      'f1':{'I':None,'A':None,'B':None,'X0':None,'S':None},
      'f2':{'A0':None,'A1':None}
      }
    for iRow in range( gfit[3].rowCount() - 1):
      row = gfit[3].row( iRow )
      functionID,parName=row['Name'].split('.')
      parameters[functionID][parName] = row['Value']
    
    workspaces={'X values':gfit[4].dataX(0),'Data':gfit[4].dataY(0),'Calc':gfit[4].dataY(1),'Diff':gfit[4].dataY(2)}

    #clean up
    suffixes=('_cloned','_cloned_peaks','_cloned_NormalisedCovarianceMatrix','_cloned_Parameters','_cloned_Workspace','_NormalisedCovarianceMatrix','_Workspace','_Parameters','_removedElastic','_removedAdjacent')
    #for suffix in suffixes: DeleteWorkspace(spectrum.getName()+suffix)
    
    return parameters,workspaces

pFile = open(outParmFile,'w')
pFile.write('#1.Pixel 2.El.I 3.El.X0 4.El.A 5.El.B 6.El.S 7.Il.I 8.Il.X0 9.Il.A 10.Il.B 11.Il.S 12.B.A0 13.B.A1\n')
for ix in range(centerPixel-pixelRange, centerPixel+pixelRange+1):
  isToTheRight = False
  if ix < centerPixel: isToTheRight = True #isToTheRight refers to position of the adjacent peak with respect to the elastic line

  sp = ExtractSingleSpectrum(InputWorkspace=wsName,WorkspaceIndex=ix,OutputWorkspace=str(ix))
  spRebinned = Rebin( sp, OutputWorkspace=sp.getName()+'_'+str(dt), params=(0,dt,Max), PreserveEvents=0)

  parameters,workspaces = fitSpectrum(spRebinned, fitMethod='two BackToBackExponentials',isToTheRight=isToTheRight)
  f0=parameters['f0'];  f1=parameters['f1'];  f2=parameters['f2']
  pFile.write('%03d %f %f %f %f %f %f %f %f %f %f %f %f\n'%(ix,f0['I'],f0['X0'],f0['A'],f0['B'],f0['S'], f1['I'],f1['X0'],f1['A'],f1['B'],f1['S'] ,f2['A0'],f2['A1']))
  pFile.flush()

  outFile = outParmFile.replace('XXX','%03d'%ix)
  wFile = open(outFile,'w')
  wFile.write('# 1.TOF  2.Data 3.Calc 4.Diff\n')
  for iTOF in range( 1, len(workspaces['X values'])-1 ):
    tof = ( workspaces['X values'][iTOF] + workspaces['X values'][1+iTOF] ) / 2.0
    wFile.write('%f %f %f %f\n'%(tof,workspaces['Data'][iTOF],workspaces['Calc'][iTOF],workspaces['Diff'][iTOF] ) )
  wFile.close()