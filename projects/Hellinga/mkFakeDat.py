#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp
from random import random,randint

root='/gpfs1/scratch/jose/Hellinga'
N=150 #sequence length
Z=400  #number of structures
nXinY=['nHinH','nEinE','nHinHE','nEinHE']
#load distributions and normalize
alph=['H','E','H_H','E_E','H_E','E_H','Nt','Ct']
ssx={'H':2,'E':4,'C':1}  #numeric code for secondary structures
dis={} #dictionary of normalized distributions
rx={} #dictionary of histogram of number of placed elements of secondary structure
for x in nXinY: rx[x]=[0]*264 #initialize histogram

for x in alph:
    pf=open(root+'/'+x+'.dis','r')
    dis[x]=chomp(pf.readlines())
    sum=0
    for i in range(len(dis[x])):
        dis[x][i]=float(dis[x][i].split()[1])
        sum+=dis[x][i]
    for i in range(len(dis[x])): dis[x][i]=dis[x][i]/sum #normalize the distribution
    pf.close()

#calculate mean helix length and mean strand length
avL={'H':0.0,'E':0.0}
for x in avL.keys():
    for i in range(len(dis[x])): avL[x]+=i*dis[x][i]
#for alpha/beta proteins, the ratio of helix to strand content is
#approximately equal (27%). Since strands are shorter, we have to put
#more number of strands than helixes
pHE=avL['E']/(avL['H']+avL['E'])

#extend to 300 residues long
for x in alph:
    for i in range(len(dis[x]),N): dis[x].append(0.0)


tofile={'H':'a','E':'b','X':'c'}
for x in ['H','E','X']: #H:alpha protein E:beta X:alpha/beta
    for i in range(Z):
        nplaced={'H':0,'E':0} #number of placed helixes
        m=0    #assigned residues
        mm=N   #remaining residues to assign
        header=tofile[x]+'%04d'%(i)
        print header
        d=root+'/'+tofile[x]+'/'+header
        os.system('mkdir -p '+d)
        p=open(d+'/seq.dat','w')
        
        #pick one N-terminal length
        while(1):
            n=randint(0,N-1)
            if random()<dis['Nt'][n]: #accept move
                for i in range(n):
                    m+=1  ;  mm-=1
                    p.write('%5d   VAL%5d   10\n'%(m,ssx['C']))
                break

        #select secondary structure lengths and interior loops
        nexty='' #variable holding the secondary structure type to be placed
        if x=='X':#select helix or strand for alpha/beta class with same probability
            if random()<pHE: nexty='H'
            else           : nexty='E'
        else: nexty=x
        y=nexty
        while(mm>3):
            #select secondary structure
            while(1):
                n=randint(0,mm)  #pick one helix or strand length
                if random()<dis[y][n]: #accept move
                    nplaced[y]+=1
                    for i in range(n):
                        m+=1  ;  mm-=1
                        p.write('%5d   VAL%5d   10\n'%(m,ssx[y]))
                    break
            #select interior loop, or maybe the loop will be the C-terminal
            while(1):
                #the loop may be the C-terminal
                if random()<dis['Ct'][mm]:
                    for i in range(mm):
                        m+=1
                        p.write('%5d   VAL%5d   10\n'%(m,ssx['C']))
                        mm=0
                    break
                if mm==0: break
                #the loop is an interior loop
                n=randint(0,mm)
                #determine which secondary structure will follow
                if x=='X':#select helix or strand for alpha/beta class with same probability
                    if random()<pHE: nexty='H'
                    else           : nexty='E'
                else: nexty=x
                if random()<dis[y+'_'+nexty][n]: #pull the correct loop length distribution
                    for i in range(n):
                        m+=1  ;  mm-=1
                        p.write('%5d   VAL%5d   10\n'%(m,ssx['C']))
                    break
                y=nexty #update the secondary structure type
                
        #assign remnant residues as loop
        for i in range(0,mm):
            m+=1
            p.write('%5d   VAL%5d   10\n'%(m,ssx['C']))
        p.close()

        #update rx
        if x=='H':
            rx['nHinH'][ nplaced['H'] ]+=1
        elif x=='E':
            rx['nEinE'][ nplaced['E'] ]+=1
        else:
            rx['nHinHE'][ nplaced['H'] ]+=1
            rx['nEinHE'][ nplaced['E'] ]+=1

#output distribution of number of placed elements of secondary structure
for x in nXinY:
    pf=open(x+'.dis','w')
    for y in range(N):  pf.write('%3d %4d\n'%(y,rx[x][y]))
    pf.close()

sys.exit(0)
