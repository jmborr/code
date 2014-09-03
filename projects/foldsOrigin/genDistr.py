#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp

#f contains fasta-like entries for each header
f='/gpfs1/scratch/jose/foldsOrigin/ss.txt'

alph=['H','E','-']  #H:helix E:strand -:(interior)loop
px={}
hx={}
N=264 #maximum assumed length for a contiguous stretch
for x in alph:
    px[x]=re.compile(x+'+') #regular expression to detect contigous stretches
    hx[x]=[0.0]*N           #initialize histogram
pN=re.compile('^(-+)[^-]')  #regexp for N-terminal
pC=re.compile('[^-](-+)$')  #regexp for C-terminal
hx['L']=[0.0]*N 

pf=open(f,'r')
l=pf.readline()
while(l):
    if l[0]=='>':
        print l
        ll=''
        l=chomp(pf.readline())
        if l=='': l=chomp(pf.readline())#there's a (stupid) blank line in between entries
        while not l[0]=='>':
            ll=ll+l #glue all lines of the fasta entry
            l=chomp(pf.readline())
            if l=='': l=chomp(pf.readline())
            if not l: break
        for x in alph: #same routine for H,E,-
            a=px[x].findall(ll) #a is a list with all contiguos stretches of the same type
            if x=='-':
                if len(a)>2: a=a[1:-1] #remove termini for loops
                else: a=[]             #protein has no definite secondary structure
            for b in a:
                hx[x][len(b)]+=1
        #special treatment for termini
        if pN.search(ll): hx['L'][len( pN.search(ll).group(1) )]+=1
        if pC.search(ll): hx['L'][len( pC.search(ll).group(1) )]+=1
        
pf.close()

#output the four histograms
for x in alph: #same routine for H,E,-
    if x=='-': pf=open('_.dis','w') #use '_' instead of '-' in filename
    else: pf=open(x+'.dis','w')
    for y in range(N):  pf.write('%3d %4d\n'%(y,hx[x][y]))
    pf.close()
#special treatment for termini    
pf=open('L.dis','w')
for y in range(N): pf.write('%3d %4d\n'%(y,hx['L'][y]))
pf.close()

sys.exit(0)
