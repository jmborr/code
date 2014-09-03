#!/usr/bin/python


import sys,os,re
from utilities.small_utilities import chomp
from spicker.spickerResultsManager import spickOut
from math import sqrt
from inputArgs.inputArgs import deglobb

TMco=0.40
Zco=-3.0
dL=10


files=('/gpfs1/active/jose/projects/nonhomol_pdbJul05/all_alpha/summaryMoreABelow200Filt.all','/gpfs1/active/jose/projects/nonhomol_pdbJul05/all_beta/summaryMoreBBelow200Filt.all','/gpfs1/active/jose/projects/nonhomol_pdbJul05/alpha_beta/summaryMoreAplusBBelow200Filt.all')

#For a particular NSS, let's suppose there are N counts and M are above the TM (or Z) cutoff. Then an estimation of the error of the M positives is sqrt(M*(N-M))
for file in files:
    all=[0] * 100
    aboveTM=[0] * 100
    aboveZ=[0] * 100
    allL=[0] * 100
    aboveL=[0] * 100
    print '##### '+file
    print '#NSS, all[NSS], N(TM | NSS), Pf(TM | NSS), error, N(Z | NSS), Pf(Z | NSS), error, L, allL[L], N(TM | L), Pf(TM | L), error'
    for line in open(file,'r').readlines():
        if line[0]=='#': continue
        items=line.split()
        L=int(items[1])
        if L<200:
            NSS=int(items[2])
            Z=float(items[16])
            tm=float(items[20])
            all[NSS]+=1
            allL[int(L/dL)]+=1
            if tm >= TMco:
                aboveTM[NSS]+=1
                aboveL[int(L/dL)]+=1
            if Z<Zco:
                aboveZ[NSS]+=1

    for i in range(0,22):
        if all[i]==0: all[i]=1   #pseudo-count to avoid dividing by zero
        if allL[i]==0: allL[i]=1
        line='%2d %3d'%(i, all[i])
        x=(1.*aboveTM[i])/all[i]
        line+=' %3d %5.1f %5.1f'%( aboveTM[i], 100*x, 100*sqrt(x*(1-x)/all[i]) )
        x=(1.*aboveZ[i])/all[i]
        line+=' %3d %5.1f %5.1f'%( aboveZ[i], 100*x, 100*sqrt(x*(1-x)/all[i]) )
        line+=' %4d %3d'%( i*dL, allL[i] )
        x=(1.*aboveL[i])/allL[i]
        line+=' %3d  %5.1f %5.1f'%( aboveL[i], 100*x, 100*sqrt(x*(1-x)/allL[i]) )
        print line

    
