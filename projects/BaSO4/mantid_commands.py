expdir='/projects/research/BaSO4/expdata'
simdir='/projects/research/BaSO4/andrew_simulation'
LoadDaveGrp(Filename='%s/resolution.dat'%expdir,OutputWorkspace='res',XAxisUnits='DeltaE',YAxisUnits='Empty',IsMicroEV=1)
e_range=(-0.1,0.0004,0.5)
Rebin(InputWorkspace='res',OutputWorkspace='res',Params=e_range)

for temp in '230 240 260 282 300'.split():
    LoadDaveGrp(Filename='%s/q%s.dat'%(expdir,temp),OutputWorkspace='exp%s'%temp,XAxisUnits='DeltaE',YAxisUnits='Empty',IsMicroEV=1)
    LoadSassena(Filename='%s/fqt_inc_T%s_water.h5'%(simdir,temp), OutputWorkspace='sim%s'%temp, TimeUnit=0.1)
    SassenaFFT(InputWorkspace='sim%s'%temp,FFTonlyRealPart=1,DetailedBalance=1,Temp=float(temp))
    SumSpectra(InputWorkspace='sim%s_sqw'%temp,OutputWorkspace='sumQ%s'%temp)
    Scale(InputWorkspace='sumQ%s'%temp,OutputWorkspace='sumQ%s'%temp, Factor=1e-9, Operation='Multiply')
    Rebin(InputWorkspace='sumQ%s'%temp,OutputWorkspace='sumQ%s'%temp,Params=e_range)
    fitstr='name=TabulatedFunction,Workspace=res,Scaling=1;(composite=Convolution;name=TabulatedFunction,Workspace=res,Scaling=1,ties=(Scaling=1);name=TabulatedFunction,Workspace=sumQ%s,Scaling=1);name=LinearBackground,A0=1e-05,A1=1e-05'%temp  
    SaveNexus(InputWorkspace='exp%s_Workspace'%temp,Filename='%s/fit2exp_T%s_water_Workspace.nxs'%(simdir,temp),Title='Structure factors for the fit of fqt_inc_T%s_water.h5 to experiment'%(temp))
    SaveNexus(InputWorkspace='exp%s_Parameters'%temp,Filename='%s/fit2exp_T%s_water_Parameters.nxs'%(simdir,temp),Title='Parameters for the fit of fqt_inc_T%s_water.h5 to experiment'%(temp))
    SaveNexus(InputWorkspace='exp%s_NormalisedCovarianceMatrix'%temp,Filename='%s/fit2exp_T%s_water_NormalisedCovarianceMatrix.nxs'%(simdir,temp),Title='Normalised covariance matrix for the fit of fqt_inc_T%s_water.h5 to experiment'%(temp))


    
