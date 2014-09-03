#!/usr/bin/python

import sys
from spicker.spickerTMresultsManager import spickOut
from utilities.small_utilities import chomp

#listf=sys.argv[1] #first argument
#root=sys.argv[2]
listf='/gpfs1/scratch/jose/checkSpickerTM/input/list'
root='/gpfs1/scratch/jose/checkSpickerTM/spickerTMmodified'

print 'header length densest    biggest       bestTM        bestRMSD'
for header in chomp(open(listf,'r').readlines()):
    dir=root+'/'+header[1]+'/'+header  #;print 'dir=',dir
    s=spickOut(dir=dir)
    if s.readError:  print s.readError
    else:
        idd=s.getDensestID(top=5) #cluster if of densest cluster
        idt=s.rankIDsByTMtoNat(list=s.ids[0:5])[0] #cluster if with best TM score
        idr=s.rankIDsByRmsdtoNat(list=s.ids[0:5])[0]
        line=header+' %3d %2d %5.3lf %4.1lf %2d %5.3lf %4.1lf %2d %5.3lf %4.1lf %2d %5.3lf %4.1lf'%(s.length,idd,s.tms[idd],s.rmsds[idd],1,s.tms[1],s.rmsds[1],idt,s.tms[idt],s.rmsds[idt],idr,s.tms[idr],s.rmsds[idr])
        print line
