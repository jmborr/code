#!/usr/bin/python

import os,sys,re
from utilities.small_utilities import chomp,junkName,Bye
from TMalign.TMalignManageResults import TMalignOut
from seq.alignYangManageResults import gappedSeq
from math import sqrt

SCroot='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
templd='/gpfs1/archive/adrian/dat/pdb_aug06.apr06' #template directory
pdbd=templd+'/CA'
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'

junk=junkName()

n_overlap=[0.0]*300   #number of insertions of a given length
av_overlap=[0.0]*300  #average overlap for insertions of a given length
sig_overlap=[0.0]*300 #standar deviation of the overlap for insertions of a given length

M=0;
for line in open(scd+'/findCommonTemplate.dat','r').readlines():
    if line[0]=='#': continue #line is a comment
    (short,long,templ)=line.split()
    spdb=pdbd+'/'+short+'.pdb' ; lpdb=pdbd+'/'+long+'.pdb' ; tpdb=pdbd+'/'+templ+'.pdb'
    pdbs={short:spdb,long:lpdb}
    #gather info on the alignment between short and long 
    cmd='TMalign '+spdb+' '+lpdb+' > '+junk  ;  os.system(cmd)
    alg=TMalignOut(junk)
    constitutive=alg.nAligned               #number of aligned residues
    alg.initIsAlignedArray() #initialize alg.isAlined dictionary
    arrayOfAligned={short:alg.isAligned[1],long:alg.isAligned[2]}
    p=re.compile('0+')
    for header in (short,long):
        onezero=''
        for i in range( len(arrayOfAligned[header]) ):
            onezero+=`arrayOfAligned[header][i]`
        pdb=pdbs[header]
        cmd='TMalign '+pdb+' '+tpdb+'  -L '+`constitutive`+' > '+junk  ;  os.system(cmd)
        alg2=TMalignOut(junk)
        alg2.initIsAlignedArray()        
        iterator= p.finditer(onezero)
        if iterator:
            for match in iterator:
                span=match.span()
                l=span[1]-span[0]
                n=0
                for i in range(span[0],span[1]):
                    if alg2.isAligned[1][i]==1: n+=1 #template covers this insertion site
                overlap=(1.0*n)/l
                n_overlap[l]+=1
                av_overlap[l]+=overlap
                sig_overlap[l]+=overlap*overlap
    #if M>100:break #this for debuggin purposes
    #M+=1

print '#insertion-length  number-of-insertions average-overlap-with-template standar-deviation'
for i in range( len(n_overlap) ):
    if n_overlap[i]>0:
        av_overlap[i]=av_overlap[i]/n_overlap[i]
        sig_overlap[i]=sqrt(sig_overlap[i]/n_overlap[i] - av_overlap[i]*av_overlap[i])
        print '%3d %3d %5.3lf %6.4lf'%(i,n_overlap[i],av_overlap[i],sig_overlap[i])
    
os.system('/bin/rm '+junk)
sys.exit(0)
