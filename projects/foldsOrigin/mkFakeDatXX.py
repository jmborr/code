#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp
from random import random,randint

root='/gpfs1/scratch/jose/foldsOrigin'
N=300  #sequence length
Z=100  #number of seq.dat files to output

#load secondary structure length distributions and normalize
alph=['H','E','L','_']         #H:helix E:strand -:(interior)loop L:terminal
ssx={'H':2,'E':4,'_':1,'L':1}  #numeric code for secondary structures
dis={} #dictionary to store the normalized distributions
for x in alph:
    pf=open(root+'/'+x+'.dis','r') #open appropriate file containing histogram of lengths
    dis[x]=chomp(pf.readlines())
    sum=0
    for i in range(len(dis[x])):
        dis[x][i]=float(dis[x][i].split()[1])
        sum+=dis[x][i]
    for i in range(len(dis[x])): dis[x][i]=dis[x][i]/sum #normalize the distribution
    pf.close()
#extend the distributions to 300 residues long
for x in alph:
    for i in range(len(dis[x]),N): dis[x].append(0.0)
    
#special for Tim Barrel: do not accept strands with less than four
#residues, helixes with less than 8 res. and interior loops with less
#than four residues ***NOTE***: since the minimal helix length is eight
#residues, we will stop our assigning secondary structure if there are
#less than 8  unassigned residues
dis['E'][1]=0;dis['E'][2]=0.0;dis['E'][3]=0.0
dis['H'][4]=0;dis['H'][5]=0.0;dis['H'][6]=0.0;dis['H'][7]=0.0;
dis['_'][1]=0;dis['_'][2]=0.0;dis['_'][3]=0.0

for i in range(Z): #generate Z seq.dat-like files with N residues and nSec secondary structures
    print '********************'
    nSec=0 #number of introduced secondary structure elements
    attempts=0 #number of times we attempt to obtain a file with 8 strands and 8 helixes
    header='tb%04d'%(i) #header for the target
    d=root+'/input/tb/'+header #directory that will hold the seq.dat-like file
    os.system('mkdir -p '+d)   #create the directory

    #TimBarrel special: throw away seq.dat files not containing 8 strands plus 8 helices
    while nSec!=16:
        nSec=0 #reset nSec
        m=0    #number of assigned residues
        mm=N   #remaining number of residues to assign
        p=open(d+'/seq.dat','w') #open the seq.dat-like file
        
        #print '#pick one N-terminal length'
        while(1):
            n=randint(0,N-1)
            if random()<dis['L'][n]: #accept move
                for i in range(n):
                    m+=1  ;  mm-=1
                    p.write('%5d   ALA%5d   10\n'%(m,ssx['L']))
                break

        #Assign secondary structure and loops to the rest of the
        #sequence. If less than eigth residues remain unasigned (see above), then
        #automatically assign them as loop
        
        timeforbeta=1 #allways start secondary structure assignment with a strand
        while(mm>7):
            #pick the secondary structure element
            #alternate between beta and alpha
            if timeforbeta==1:
                y='E'
                timeforbeta=0
            else:
                y='H'
                timeforbeta=1
            #pick length of the secondary structure element
            while(1):
                n=randint(0,mm)
                if random()<dis[y][n]: #accept move
                    nSec+=1 #we introduce another secondary structure element
                    for i in range(n):
                        m+=1  ;  mm-=1
                        p.write('%5d   ALA%5d   10\n'%(m,ssx[y]))
                    break
            
            pick subsequent loop
            while(1):
                #the loop may be the C-terminal
                if random()<dis['L'][mm]:
                    for i in range(mm):
                        m+=1
                        p.write('%5d   ALA%5d   10\n'%(m,ssx['L']))
                        mm=0
                    break
                if mm==0: break
                #pick one interior loop length
                n=randint(0,mm)
                random()<dis['_'][n]
                if random()<dis['_'][n]:
                    for i in range(n):
                        m+=1  ;  mm-=1
                        p.write('%5d   ALA%5d   10\n'%(m,ssx['_']))
                    break

        #Assing the remaining residues (if any) as loop
        for i in range(0,mm):
            m+=1
            p.write('%5d   ALA%5d   10\n'%(m,ssx['_']))

        p.close()
        attempts+=1
        print header,`attempts`

sys.exit(0)
