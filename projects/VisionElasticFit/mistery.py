# Initial Guess:
# name=BackToBackExponential,I=271.487858,A=0.01,B=0.01,X0=5024.888620,S=393.321667,constraints=(A>0,B>0,S>196.660833)
# startX=4775.085702480463
# endX=5527.747977019568

# After the fit in MantidPlot:
# FittedParams={'I': 226.886431, 'A': 0.008677, 'X0': 4880.107234, 'B': 0.002034, 'S': 196.504648}

# After the fit with the Python script (note B is negative even though we constrained to be positive)
# FittedParams={'I': 631.0404666884342, 'A': 0.004825081924682037, 'X0': 5199.0811591762495, 'B': -0.03677555258585037, 'S': 350.69664884391074}

import sys
sys.path.append('/home/jmborr/devel/mantidproject/mantid/Code/debug/bin')
from mantid.simpleapi import *

sp = LoadAscii(Filename='/tmp/mistery.dat',Unit='TOF',Separator='CSV')

#funcStr = 'name=BackToBackExponential,I=4256.245644,A=0.010059,B=0.010062,X0=4076.714432,S=116.902398,constraints=(A>0,B>0,S>58.451199);name=BackToBackExponential,I=372.302715,A=0.001591,B=0.007044,X0=4905.004347,S=455.125745,constraints=(A>0,B>0,X0>4780.004347,X0<5030.004347,S>227.562872);name=LinearBackground,A0=-76.378293,A1=0.021016'

#funcStr = 'name=BackToBackExponential,I=4256,A=0.010059,B=0.010062,X0=4076,S=117,constraints=(A>0,B>0,S>58);name=BackToBackExponential,I=3725,A=0.001591,B=0.007044,X0=4905,S=455,constraints=(A>0,B>0,X0>4780,S>227);name=LinearBackground,A0=-76,A1=0.021016'

funcStr = 'name=BackToBackExponential,I=4256.245644,A=0.010059,B=0.010062,X0=4076.714432,S=116.902398,constraints=(A>0,B>0,S>58.451199);name=BackToBackExponential,I=372.302715,A=0.001591,B=0.007044,X0=4905.004347,S=455.125745,constraints=(A>0,B>0,4780.004347<X0<5030.004347,S>227.562872);name=LinearBackground,A0=-76.378293,A1=0.021016'


startX = 3411.607
endX   = 5672.455

print 'INITIAL :',funcStr
gfit = Fit(funcStr, sp, WorkspaceIndex=0, StartX=startX, EndX=endX, CreateOutput='1')

#Retrieve the optimized fitting parameters
parameters={'f0':{'I':None,'A':None,'B':None,'X0':None,'S':None},
            'f1':{'I':None,'A':None,'B':None,'X0':None,'S':None},
            'f2':{'A0':None,'A1':None}
            }
for iRow in range( gfit[3].rowCount() - 1):
    row = gfit[3].row( iRow )
    functionID,parName=row['Name'].split('.')
    parameters[functionID][parName] = row['Value']

print 'FINAL :',parameters
