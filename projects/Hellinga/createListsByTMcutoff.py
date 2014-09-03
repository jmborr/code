#!/usr/bin/python
'''

create a series of lists, each one containing a subset of all
structures in /gpfs1/scratch/jose/Hellinga/out/structures. For every
subset, there are no two structures with a TM score bigger than some
TM-score cut-off. The name of the list denotes the cut-off.

'''
import os,sys,re
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp
from jobs.job import job
from tmScore.tm_score import TMscore

scd='/gpfs1/scratch/jose/Hellinga'

#strs=chomp(open(scd+'/structures/list','r').readlines()) #list of structures
strs=chomp(open(scd+'/structures/list','r').readlines()) #list of structures

#initialize subset, tms, and list with 0.85 TM cutoff
tms={}

#within parenthesis formula of Betancourt & Skolnick, Biopolymers 59, 305-309 (2001)
Rg_max=1.5*( 3.08*150**0.3333 ) #; print Rg_max
maxtm=0.35
subset=[]
nn=0
m=0
for pdb1 in strs:
    n=0
    m+=1
    #print pdb1,m
    fullpdb1=scd+'/structures/'+pdb1 #full path
    #print 'fullpdb1=',fullpdb1
    #filter by Rg
    Rg=float(os.popen('Rg.x -a '+fullpdb1).readline().strip()) ; #print Rg
    if Rg > Rg_max: continue
    homolog=0
    tmsTemp={pdb1:{}}
    for pdb2 in subset:
        fullpdb2=scd+'/structures/'+pdb2 #full path
        #print 'fullpdb2=',fullpdb2
        
        tm=TMscore(fullpdb1,fullpdb2)
        n+=1
        #print pdb1,pdb2,tm
        if tm>maxtm: #we found a structural homologous in subset
            homolog=1
            break
        tmsTemp[pdb1][pdb2]=tm
        tmsTemp[pdb2]={}
        tmsTemp[pdb2][pdb1]=tm
    nn+=n
    #print n
    if not homolog:
        subset.append(pdb1)
        for key in tmsTemp.keys():
            if key not in tms.keys():
                tms[key]=tmsTemp[key]
            else:
                tms[key].update(tmsTemp[key]) #merge the two dictionaries

#print nn
open('subset'+'%4.2lf'%(maxtm)+'.list','w').write('\n'.join(subset))

#create lists for more relaxed TM-score cut-offs
maxtm=0.80
while maxtm>=0.30:
    prevsubset=subset
    subset=[] #this will not set prevsubset to [] too
    
    homolog=0
    for pdb1 in prevsubset:
        fullpdb1=scd+'/structures/'+pdb1
        for pdb2 in subset:
            fullpdb2=scd+'/structures/'+pdb2
            if tms[pdb1][pdb2]>maxtm:
                homolog=1
                break
        if not homolog: subset.append(pdb1)
                
    #output subset of structures
    open('subset'+'%4.2lf'%(maxtm)+'.list','w').write('\n'.join(subset))
    
    #next iteration
    maxtm=maxtm-0.05
    
sys.exit(0)
