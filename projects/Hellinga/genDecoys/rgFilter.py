#!/usr/bin/python

'''******************** Pseudocode *****************
1 open native file -> determine protein length and radius of gyration
2 set RgM=1.5*native_Rg
3 open "rep1.tra" and "rep1Rg.tra"
4 read snapshot by snapshot of "rep1.tra"
  4.1 create junkRANDOM.pdb with the snapshot
  4.2 run Rg.x junkRANDOM.pdb. If the Rg < RgM, store snapshot in "rep1Rg.tra"
**************************************************'''

import sys,os
from utilities.small_utilities import chomp,Bye,junkName
from random import randint
from inputArgs.inputArgs import inpHand
from math import sqrt
###############################################################
def rg(xyzs):
    N=len(xyzs)
    #first calculate center of mass
    ccc=[0.0,0.0,0.0]
    for xyz in xyzs:
        for i in range(3): ccc[i]+=xyz[i]
    for i in range(3): ccc[i]=ccc[i]/N
    #Now calculate radius of gyration
    v=0.0
    for xyz in xyzs:#for every particle
        for i in range(3):#x,y,z components
            d=xyz[i]-ccc[i]
            v+=d*d
    return sqrt(v/N)
###############################################################
ih=inpHand(
    'Usage: rgFilter.py ',
    '  -a _AR_rep rep1.tra file',
    '  -b _AR_nat native pdb file',
    '  -c __repR  filtered file (def=dirname(rep1.tra file)/rep1RgF.tra)'
    )
ih.parse(locals(),sys.argv)

#initialize optional input
if not repR: repR=os.path.dirname(rep)+'/rep1RgF.tra'

#create a junk file
junk=junkName()

#create native file with only CA atoms
natCA=junkName()
os.system('grep -e " CA " '+nat+' > '+natCA)
          
#get sequence length
L=int(os.popen('wc -l '+natCA).readline().split()[0])

#get radius of gyration cut-off
RgM=1.5 * float( chomp(os.popen('Rg.x -a '+natCA).readline()) )

#get amino and numbering info
prefixes=[]
p=open(natCA,'r')
for i in range(0,L):
    prefixes.append( p.readline()[0:31] )
p.close()

#open rep1.tra, rep1Rg.tra, and read trayectory by trayectory
p=open(rep,'r')
pout=open(repR,'w')
header=p.readline()
while(header!=''):
    axyzs=[]#initialize snapshot coordinates, as string type 
    xyzs=[] #initialize snapshot coordinates, as float type
    for i in range(0,L): axyzs.append(p.readline())
    for i in range(0,L):
        xyz=axyzs[i].split()
        for i in range(3):xyz[i]=float(xyz[i])
        xyzs.append(xyz)
    #calculate radius of gyration of the snapshot
    Rg=rg(xyzs)
    #store snapshot in filtered trayectory if compact
    if Rg<RgM:
        header=chomp(header)+' Rg=%5.1f\n'%(Rg) #attach Rg value
        pout.write(header)
        pout.writelines(axyzs)
    #read next snapshot-header
    header=p.readline()
p.close()
pout.close()

#cleanup
os.system('/bin/rm '+junk)
os.system('/bin/rm '+natCA)

sys.exit(0)
