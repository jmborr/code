rootd='/projects/research/POSS/me8t8/charmm/nmolec_1'
for T in (100,125,150,175,200,225,250,300,350):
	LoadSassena(Filename='%s/T_%d/fqt_inc.h5'%(rootd,T),TimeUnit=0.1,OutputWorkspace='inc%d'%T)
	SassenaFFT(InputWorkspace='inc%d'%T,FFTonlyRealPart=1,DetailedBalance=1,Temp=T)
	
rootd='/projects/research/POSS/me8t8/expdata'
for T in (100,150,175,200,225,250,300,350):
	print 'T',T
	#LoadDaveGrp(Filename='%s/exp%dK.dat'%(rootd,T),XAxisUnits='DeltaE',YAxisUnits='Empty',IsMicroEV=1,OutputWorkspace='exp%dK'%T)
	#Rebin(InputWorkspace='exp%dK' %T,Params=[-0.14,0.0004,0.14],OutputWorkspace='exp%dK_2' %T)
	SaveNexus(InputWorkspace='exp%dK_2' %T,Filename='/projects/research/POSS/me8t8/expdata/exp%dK.nxs'%T)