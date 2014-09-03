#!/usr/bin/python

import os,sys
from utilities.small_utilities import Bye,chomp

all='/gpfs1/scratch/jose/nonhomol_pdbJul05/analysis/summaryMoreBelow200.all'
root='/gpfs1/scratch/jose/createInputsTasser/'

n={}
for i in range(201): n[i]=0.0
m={}
for i in range(201): m[i]=0.0

for line in open(all,'r').readlines():
    if line[0]=='#': continue #comment line
    items=line.split()    
    header=items[0] #;print header
    tm=items[19]
    pdbf=root+'/'+header[1]+'/'+header+'/abinitio/best.pdb'
    cs=chomp(os.popen('/library/jose/code/python/tertiary/contactMap.py -a '+pdbf).readlines())
    for c in cs:
        i,j=c.split(); l=abs(int(i)-int(j)) #;print i,j,l
        if l>0: n[l]+=1
    pdbf=root+'/'+header[1]+'/'+header+'/abinitio/CA'
    cs=chomp(os.popen('/library/jose/code/python/tertiary/contactMap.py -a '+pdbf).readlines())
    for c in cs:
        i,j=c.split(); l=abs(int(i)-int(j))
        if l>0: m[l]+=1

for i in range(201):
    print '%3d %6d %6d'%(i,n[i],m[i])

sys.exit(0)
