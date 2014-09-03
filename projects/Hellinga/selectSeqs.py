#!/usr/bin/python

import os,sys,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import Bye,chomp
from math import sqrt

#parse input options
inp=inpHand('Usage: selectSeqs.py [options]\ngiven rankSeqsBySSsimilarity.dat, select the sequences based on the same Z-score cutoff for all properties\n',
            ' -a _RA_rankssf rankSeqsBySSsimilarity.dat',
            ' -b _A_filtf file for output filtered list of sequences (def:STDOUT)',
            ' -c __Zco Z-score cutoff (def:0.0)',
            ' -d __nfilt output at least nfilt sequences (def:500)'
            )
inp.parse(locals(),sys.argv)

if not filtf:filtp=sys.stdout
else: filtp=open(filtf,'w')
if not Zco:Zco=0.0
if not nfilt: nfilt=500

#read rankSeqsBySSsimilarity.dat file
propl=('score', 'ss_seqid', 'seqid', 'confidence', 'ss_confidence')
maths={'averages':{} ,
       'std-devs':{} }

all={}
for line in chomp(open(rankssf,'r').readlines()):
    if line[0]=='#':continue #line is a comment line
    items=line.split()
    if len(items)!=1+len(propl): Bye('ERROR: number of properties does not agree')
    header=items[0]
    all[header]={'line':line}
    for i in range(len(items)-1): all[header][ propl[i] ]=float(items[i+1])
nall=len(all)

#calculate averages and standard deviation and Z-score
for math in maths.keys(): #initialize each math property for each prop to zero
    for prop in propl:  maths[math][prop]=0.0
    
for header in all.keys():
    for prop in propl:
        x=all[header][prop]
        maths['averages'][prop]+=x
        maths['std-devs'][prop]+=x*x

for prop in propl:
    maths['averages'][prop]=float(maths['averages'][prop])/nall
    x=maths['averages'][prop]
    maths['std-devs'][prop]=sqrt(float(maths['std-devs'][prop])/nall-x*x)

#calculate Z-scores. normalize so that the more negative the z-score,
#the better no matter what property we're dealing with
signs={'score':1, 'ss_seqid':-1, 'seqid':-1, 'confidence':-1,'ss_confidence':-1}
zscores={}
for header in all.keys():
    zscores[header]={}
    for prop in propl:
        av=maths['averages'][prop]
        std=maths['std-devs'][prop]
        zscores[header][prop]=signs[prop]*(all[header][prop]-av)/std

#select only those below cutoff Z-score
buf=''
while(1):
    nout=0
    buf=''
    for header in all.keys():
        output=True
        for prop in propl:
            if zscores[header][prop]>Zco:
                output=False
                break
        if output:
            nout+=1
            buf+=all[header]['line']+'\n'
    if nout<nfilt and nout<nall: Zco+=0.01
    else:
        print Zco
        break
#
filtp.write(buf)
filtp.close()
    
