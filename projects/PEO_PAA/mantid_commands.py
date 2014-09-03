# MonojoyHundredThousandTimeSteps I(Q,t) have 200,000 time points, way
# too many. We calculate S(Q,E) from I(Q,t) in the [-1ns,1ns]
cwd='/projects/research/PEO_PAA/MonojoyHundredThousandTimeSteps'
for D in '02 04 06 08 10'.split():
  LoadSassena(Filename='%s/charge%sfqt.h5'%(cwd,D), OutputWorkspace='charge'+D, TimeUnit=0.1)
  Rebin(InputWorkspace='charge%s_fqt.Re'%D, OutputWorkspace='charge%s_fqt.Re'%D, Params=[-2000,0.1,2000])
  SassenaFFT(InputWorkspace='charge'+D, FFTonlyRealPart=1, DetailedBalance=1, Temp=300)
  ExtractSingleSpectrum(InputWorkspace='charge%s_sqw'%D, OutputWorkspace='w'+D, WorkspaceIndex=0)
  SmoothData(InputWorkspace='w'+D,OutputWorkspace='w'+D, NPoints=10)
  Scale(InputWorkspace='w'+D,OutputWorkspace='w'+D, Factor=0.000000023, Operation='Multiply')
