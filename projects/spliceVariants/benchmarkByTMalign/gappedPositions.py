#!/usr/bin/python

import os,sys
from utilities.small_utilities import chomp,junkName,Bye
from TMalign.TMalignManageResults import TMalignOut
from seq.alignYangManageResults import gappedSeq

nmax=300 #sequence with maximum number of residues to consider

algsf='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign/alignmentsByTM.dat'
pin=open(algsf,'r')
algsidx='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign/alignmentsByTM.idx'

pidx=open(algsidx,'r') ; line=pidx.readline()
#junkN=2 ; N=0
while line :
    idx=int(line.split(':')[0]) #;print 'idx=',idx#index in alignmentsByTM.dat
    pin.seek(idx)
    headerline=pin.readline()  #info on the two sequences being read
    algOut=TMalignOut(pin)     #object to store alignment
    gs={}
    for i in (1,2): gs[i]=gappedSeq(algOut.seq[i]) #create gapped sequence objects
    L={}                       #lengths of each sequence    
    for i in (1,2): L[i]=len(gs[i].gs.replace('-',''))
    if L[1]>nmax or L[2]>nmax:
        line=pidx.readline() #next index pinpointing next alignment
        continue             #go to next iteration
    sys.stdout.write(headerline)        
    for i in (1,2):
        seq=gs[i]
        gaps=seq.extractGaps() #;print gaps
        poss=seq.gapPositionsUngappedSeq() #;print poss #positions of the beginning of the gaps
        for j in range(len(poss)):
            pos=1+poss[j] #shift index by one to abandon C-array convention
            lgap=len(gaps[j])
            if L[i]-pos<pos: pos=L[i]-pos #position counting from the end        
            sys.stdout.write('%3d %2d %3d\n'%(pos,(100*pos)/L[i],lgap))
    line=pidx.readline() #next index pinpointing next alignment
#    if N>junkN:Bye('')
#    else: N+=1
pin.close()
pidx.close()
sys.exit(0)
