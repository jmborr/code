#!/usr/bin/python

'''******************** Pseudocode *****************
1 open native seq.dat -> extract secondary assigments
2 open rep1Rg.tra for reading, rep1RgSS.tra for writing
3 read snapshot by snapshot of "rep1Rg.tra"
  3.1 obtain secondary assignments from snapshot and compare to native
  3.2 if assignments sufficient similar, then output
**************************************************'''

import sys,os
from utilities.small_utilities import chomp,Bye,junkName
from inputArgs.inputArgs import inpHand
from secondary.pdb2ss import smooth,outputSec0
from secondary.seqdat2ss import translateSec

ih=inpHand(
    'Usage: filterByssSimilarity.py [options] ',
    '  -a _AR_rep rep1Rg.tra file',
    '  -b _AR_natseqdatf native seq.dat file',
    '  -c __repSS  filtered file (def=dirname(rep1Rg.tra file)/rep1RgSS.tra)',
    '  -d __minsim minimal similarity (def:0.75)'
    )
ih.parse(locals(),sys.argv)

#initialize optional input
if not repSS: repSS=os.path.dirname(rep)+'/rep1RgSS.tra'
if not minsim: minsim=0.75
else: minsim=float(minsim)

natss=translateSec(natseqdatf) #native secondary assignments
L=len(natss) #get sequence length

#open rep1Rg.tra, and read snapshot by snapshot
p=open(rep,'r')
pout=open(repSS,'w')
header=p.readline()
while(header!=''):
    axyzs=[]#initialize snapshot coordinates, as string type 
    xs=[] ; ys=[] ; zs=[] #initialize snapshot coordinates, as float type
    for i in range(0,L): axyzs.append(p.readline())
    for i in range(0,L):
        x,y,z=axyzs[i].split()
        xs.append(float(x)) ; ys.append(float(y)) ; zs.append(float(z))
    ss=smooth(outputSec0(xs,ys,zs)).replace('T','C') #snapshot secondary assignments
    sim=0.0
    for i in range(L):
        if ss[i]==natss[i]: sim+=1.0
    sim=sim/L
    if sim>minsim:
        header=chomp(header)+' ss=%4.2lf\n'%(sim)
        pout.write(header)
        pout.writelines(axyzs)
    header=p.readline()
    
p.close()      
pout.close()
