#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp
from random import random,randint

root='/gpfs1/scratch/jose/foldsOrigin'
N=300 #sequence length
Z=200  #number of structures
#load distributions and normalize
alph=['H','E','_','L']         #H:helix E:strand -:(interior)loop L:terminal
ssx={'H':2,'E':4,'_':1,'L':1}  #numeric code for secondary structures
dis={} #dictionary of normalized distributions
for x in alph:
    pf=open(root+'/'+x+'.dis','r')
    dis[x]=chomp(pf.readlines())
    sum=0
    for i in range(len(dis[x])):
        dis[x][i]=float(dis[x][i].split()[1])
        sum+=dis[x][i]
    for i in range(len(dis[x])): dis[x][i]=dis[x][i]/sum #normalize the distribution
    pf.close()
#extend to 300 residues long
for x in alph:
    for i in range(len(dis[x]),N): dis[x].append(0.0)

#generate one alpha file

tofile={'H':'a','E':'b','X':'c'}
for x in ['H','E','X']: #H:alpha protein E:beta X:alpha/beta
    for i in range(Z):
        m=0    #assigned residues
        mm=N   #remaining residues to assign
        header=tofile[x]+'%04d'%(i)
        print header
        d=root+'/'+tofile[x]+'/'+header
        os.system('mkdir -p '+d)
        p=open(d+'/seq.dat','w')
        while(1):
            n=randint(0,N-1)  #pick one N-terminal length
            if random()<dis['L'][n]: #accept move
                for i in range(n):
                    m+=1  ;  mm-=1
                    p.write('%5d   ALA%5d   10\n'%(m,ssx['L']))
                break
        while(mm>3): #less than three residues is automatically assigned as loop
            #print mm
            while(1):
                if x=='X': #select helix or strand for alpha/beta class
                    if random()<0.5: y='H'
                    else           : y='E'
                else: y=x
                n=randint(0,mm)  #pick one helix or strand length
                if random()<dis[y][n]: #accept move
                    for i in range(n):
                        m+=1  ;  mm-=1
                        p.write('%5d   ALA%5d   10\n'%(m,ssx[y]))
                    break
            while(1):
                if random()<dis['L'][mm]: #the loop may be the C-terminal
                    for i in range(mm):
                        m+=1
                        p.write('%5d   ALA%5d   10\n'%(m,ssx['L']))
                        mm=0
                    break
                if mm==0: break
                n=randint(0,mm)   #pick one loop lenght
                random()<dis['_'][n]
                if random()<dis['_'][n]: #the loop may be the C-terminal
                    for i in range(n):
                        m+=1  ;  mm-=1
                        p.write('%5d   ALA%5d   10\n'%(m,ssx['_']))
                    break
                
        for i in range(0,mm):
            m+=1
            p.write('%5d   ALA%5d   10\n'%(m,ssx['_']))
        p.close()    

sys.exit(0)
