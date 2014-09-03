#!/usr/bin/python

import sys
from spicker.spickerTMresultsManager import spickOut
from utilities.small_utilities import chomp

#listf=sys.argv[1] #first argument
#root=sys.argv[2]
listf='/gpfs1/scratch/jose/checkSpickerTM/input/list'
root='/gpfs1/scratch/jose/checkSpickerTM/spickerTMmodified2000'

print '#header length densest   bestTM        bestRMSD'
for header in chomp(open(listf,'r').readlines()):
    ranks={}
    dir=root+'/'+header[1]+'/'+header  #;print 'dir=',dir
    nat='/gpfs1/scratch/jose/checkSpickerTM/input/'+header[1]+'/'+header+'/CA'
    s=spickOut(dir=dir,nat=nat)
    if s.readError:  print s.readError
    else:
        ranked=s.rankIDsByDens()  #cluster ID's ranked by decreasing density
        for i in range(0,5): ranks[ranked[i]]=1+i #returns rank for a given cluster ID
        idd=ranked[0]             #cluster ID of densest cluster
        for id in ranked[0:5]: s.tms[id]=s.TMtoNat(id) #recalculate TM-score
        idt=s.rankIDsByTMtoNat(list=ranked[0:5])[0] #cluster if with best TM score
        rt=ranks[idt] #rank of cluster with best TM score
        idr=s.rankIDsByRmsdtoNat(list=ranked[0:5])[0] #cluster ID with best RMSD
        rr=ranks[idr] #rank of cluster with best RMSD
        line=header+' %3d  1 %5.3lf %4.1lf %2d %5.3lf %4.1lf %2d %5.3lf %4.1lf'%(s.length,s.tms[idd],s.rmsds[idd],rt,s.tms[idt],s.rmsds[idt],rr,s.tms[idr],s.rmsds[idr])
        print line
