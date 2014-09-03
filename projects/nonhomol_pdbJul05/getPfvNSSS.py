#!/usr/bin/python


import sys,os,re
from utilities.small_utilities import chomp
from spicker.spickerResultsManager import spickOut
from math import sqrt
from inputArgs.inputArgs import deglobb

TMco=0.42
Zco=-3.0
dL=10

all=[0] * 100
aboveTM=[0] * 100
aboveZ=[0] * 100
allL=[0] * 100
aboveL=[0] * 100
files=('/gpfs1/active/jose/projects/nonhomol_pdbJul05/all_alpha/summaryMoreA.all','/gpfs1/active/jose/projects/nonhomol_pdbJul05/all_beta/summaryMoreB.all','/gpfs1/active/jose/projects/nonhomol_pdbJul05/alpha_beta/summaryAplusBMore.all')

#For a particular NSS, let's suppose there are N counts and M are above the TM (or Z) cutoff. Then an estimation of the error of the M positives is sqrt(M*(N-M))
for file in files:
    print '##### '+file
    print '#NSS,  all[NSS],  Pf(TM | NSS), erro,  Pf(Z | NSS), error, L,  allL[L],  Pf(TM | L), error'
    for line in open(file,'r').readlines():
        if line[0]=='#': continue
        items=line.split()
        L=int(items[1])
        if L<200:
            NSS=int(items[2])
            Z=float(items[16])
            tm=float(items[19])
            all[NSS]+=1
            allL[int(L/dL)]+=1
            if tm >= TMco:
                aboveTM[NSS]+=1
                aboveL[int(L/dL)]+=1
            if Z<Zco:
                aboveZ[NSS]+=1

    for i in range(0,22):
        if all[i]==0: all[i]=1
        if allL[i]==0: allL[i]=1
        print '%2d %3d %5.1f %5.1f %5.1f %5.1f %4d %3d %5.1f %5.1f'%\
              ( i, all[i], (100.0*aboveTM[i])/all[i], 100*sqrt(aboveTM[i]*(all[i]-aboveTM[i])), (100.0*aboveZ[i])/all[i], 100*sqrt(aboveZ[i]*(all[i]-aboveZ[i])), i*dL, allL[i], (100.0*aboveL[i])/allL[i] ), 100*sqrt(aboveL[i]*(allL[i]-aboveL[i]))

    
