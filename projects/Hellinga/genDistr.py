#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp

#f contains fasta-like entries for each header
f='/gpfs1/scratch/jose/foldsOrigin/ss.txt'
cad='/gpfs1/archive/adrian/dat/pdb_may06.apr06.aux/CA'
rgx='/gpfs1/active/jose/code/cpp/pdb_utils/Rg.x'

alph=['H','E','-']  #H:helix E:strand -:(interior)loop
px={}
hx={}
Nx={'H':0.0,'E':0.0,'-':0.0}
Ntarget=150
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
                header=l[1:6]

        #filter according to protein length and radius of gyration
        caf=cad+'/'+header+'.pdb' #CA file
        L=int(os.popen('grep " CA " '+caf+' | wc -l').readline().strip())
        if L<Ntarget*0.5 or Ntarget>150: continue #protein length not in target length
        rg,rg0=os.popen(rgx+' -a '+caf+' -b ').readline().split()
        if float(rg) > 1.25 * float(rg0): continue #protein is not compact

        for x in alph: #same routine for H,E,-
            a=px[x].findall(ll) #a is a list with all contiguos stretches of the same residue
            Nx[x]=len(''.join(a))
            if x=='-':
                if len(a)>2: a=a[1:-1] #remove termini for loops
                else: a=[]             #protein has no definite secondary structure
            for b in a:
                hx[x][len(b)]+=1
                
        #output secondary structure content
        L=len(ll)
        if L>0:
            for x in Nx.keys(): Nx[x]=(1.0*Nx[x])/L
            content='%4.2lf %4.2lf %4.2lf'%(Nx['H'],Nx['E'],1-Nx['H']-Nx['E'])
            if   Nx['E']<0.1 and Nx['H']>0.19: content='a '+content
            elif Nx['H']<0.1 and Nx['E']>0.19: content='b '+content
            elif Nx['H']>0.1 and Nx['E']>0.1: content='c '+content
            print content            
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
