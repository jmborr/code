#!/usr/bin/python
'''

This script will cathegorize structures accoding to their long
strands. It takes into account the fact of simbolic links to
out_filtered. It can be equally applied to beta and to alpha/beta
structures

'''
import os,sys,re,copy
from inputArgs.inputArgs import inpHand,deglobb,extractHeader
from utilities.small_utilities import Bye,chomp,junkName
from secondary.seqdat2ss import translateSec
from random import randint
from copy import copy

ih=inpHand('Usage: checkLongStrands.py [options]\nExample:./checkLongStrands.py -a /gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/out -b /gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/out/b.list -c /gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/structures -d /gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/structures/b.filtered.list -e /gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/filtered/real2symb.list',
           ' -a _RA_rootd root directory for seq.dat files',
           ' -b _RA_listf list of headers',
           ' -c _RA_rootdstr structures directory',
           ' -d _RA_filtl list of filtered structures',
           ' -e _RA_real2symbf file relating each filtered structure to its symbolic name'
           )
ih.parse(locals(),sys.argv)

#read filtered structures
fs=chomp(open(filtl,'r').readlines())
nfilt=len(fs)

#stablish relationship between symbolic link and real file
real2symb={}
for line in open(real2symbf,'r').readlines():
    items=line.split()
    real2symb[items[0]]=items[1]

max=42
seqdathits=[]
for i in range(max): seqdathits.append([])
pat=re.compile('E+')
headers=chomp(open(listf,'r').readlines())
mult={}     #structures associated to each header
filtmult={} #filtered structures associated to each header
nstruct=0
for header in headers:
    mult[header]=chomp(os.popen('ls -1 '+rootdstr+'/'+header+'* 2>/dev/null').readlines())
    nstruct+=len(mult[header])
    filtmult[header]=[]
    for i in range( len(mult[header]) ):
        f=os.path.basename(mult[header][i])
        mult[header][i]=f
        if f in fs: #f is a filtered structure
            filtmult[header].append(f)
        
    ss=chomp(os.popen('seqdat2ss.py -a '+rootd+'/'+header+'/seq.dat').readline())
    it=pat.finditer(ss)
    for m in it:
        (b,e)=m.span() #beginning and ending+1 of the strand
        l=e-b #length of strand
        seqdathits[l].append(header)

#remove duplicates. Now seqdathits[l] will contain headers with one strand of length l
for i in range(max): seqdathits[l]={}.fromkeys( seqdathits[l] ).keys()

#cummulative and remove duplicates. seqdathits[l] will contain headers
#with one strand equal or longer than l.
l=max-2
while l:
    seqdathits[l]=seqdathits[l]+seqdathits[l+1][:]
    seqdathits[l]={}.fromkeys( seqdathits[l] ).keys() #remove duplicates
    l-=1
#seqdatNOhits[l] will contain headers with all strands smaller than l
seqdatNOhits=[]
for l in range(max): seqdatNOhits.append([])
for l in range(max):
    seqdatNOhits[l]=copy(headers) #hard copy, not a reference to headers
    for header in seqdathits[l]: seqdatNOhits[l].remove(header)

#find number of structures and filtered structures containing one
#strand of length equal or bigger than some length
structhits=[0]*max
filthits=[0]*max
buf='#strand-length number-unfiltedred-struct-with-strands-smaller-than number-filtered-str...\n'
for l in range(8,max):
    #print l,len(seqdathits[l]),seqdathits[l]
    for header in seqdathits[l]:
        structhits[l]+= len(mult[header])
        for f in mult[header]:
            if f in fs: filthits[l]+=1
    buf+='%2d %4d %4d\n'%(l,nstruct-structhits[l],nfilt-filthits[l])
open('hist.strandsSmaller.dat','w').write(buf)

#print a list of symbolic links to filtered structures with all strands smaller than some length
for l in range(9,max):
    outp=open('allStrandsShorterThan_'+`l`+'.list','w')
    buf=''
    for header in seqdatNOhits[l]:#cycle over all header having all strands smaller than l
        for f in filtmult[header]: buf+=real2symb[f]+'\n'
    outp.write(buf)
    outp.close()

#print a list of filtered structures having at least one strand of
#some length, but no strands bigger than this particular length
for l in range(9,max):
    outp=open('oneStrand_'+`l`+'_noStrandsLonger.list','w')
    buf=''
    for header in seqdathits[l]: #header containing one strand of length l or bigger
        if header not in seqdathits[l+1]:
            for f in filtmult[header]: buf+=real2symb[f]+'\n'
                
    outp.write(buf)
    outp.close()
    
sys.exit(0)


        
