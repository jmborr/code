#!/usr/bin/python
'''

Create length distributions for helix, strands, and loops. Loops have
different distribution, wether they are at the termini, or bound by
different types of secondary structures

'''
import sys,os,re
from utilities.small_utilities import chomp

#f contains fasta-like entries for each header
f='/gpfs1/scratch/jose/foldsOrigin/ss.txt'

alph=['H','E','-']  #H:helix E:strand -:(interior)loop
#below are different chunks: helix, strand, N-termini-loop, C-termini-loop, loop between helixes, loop between strands, , loop-preceeded-of-helix-and-followed-by-strand, loop-preceeded-of-strand-and-followed-by-helix
alph2=['HHHH+','EE+','^-+','-+$','HHHH-+HHHH','EE-+EE','HHHH-+EE','EE-+HHHH']
px={}
hx={}
N=264 #maximum assumed length for a contiguous stretch
for x in alph2:
    px[x]=re.compile(x) #regular expression to detect contigous stretches
for x in alph2:
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
        for x in alph2: #same routine for H,E,-
            iterator = px[x].finditer(ll)
            for match in iterator: hx[x][len(match.group())]+=1
        
pf.close()

#output the  histograms
alph3={'HHHH+':'H', 'EE+':'E', '^-+':'Nt', '-+$':'Ct', 'HHHH-+HHHH':'H_H', 'EE-+EE':'E_E', 'HHHH-+EE':'H_E', 'EE-+HHHH':'E_H'}
for x in alph2: #same routine for H,E,-
    pf=open(alph3[x]+'.dis','w')
    for y in range(N):  pf.write('%3d %4d\n'%(y,hx[x][y]))
    pf.close()

sys.exit(0)
