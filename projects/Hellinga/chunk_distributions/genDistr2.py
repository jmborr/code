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
cad='/gpfs1/archive/adrian/dat/pdb_may06.apr06.aux/CA'
rgx='/gpfs1/active/jose/code/cpp/pdb_utils/Rg.x'
Nx={'H':0.0,'E':0.0,'-':0.0}
nXinY=['nHinH','nEinE','nHinHE','nEinHE']
Ntarget=150

alph=['H','E','-']  #H:helix E:strand -:(interior)loop
#below are different chunks: helix, strand, N-termini-loop, C-termini-loop, loop between helixes, loop between strands, , loop-preceeded-of-helix-and-followed-by-strand, loop-preceeded-of-strand-and-followed-by-helix
alph2=['(H+)','(E+)','^(-+)','(-+)$','H(-+)H','E(-+)E','H(-+)E','E(-+)H']
px={}
hx={}
qx={}
rx={}
N=264 #maximum assumed length for a contiguous stretch
for x in alph2:
    px[x]=re.compile(x) #regular expression to detect contigous stretches
for x in alph2:
    hx[x]=[0.0]*N           #initialize histogram
for x in nXinY:
    rx[x]=[0]*N           #initialize histogram
for x in Nx.keys():
    qx[x]=re.compile(x)
    
pf=open(f,'r')
l=pf.readline()
while(l):
    if l[0]=='>':
        print l
        header=l[1:6]
        
        #store secondary structure assignment in ll
        ll=''
        l=chomp(pf.readline())
        if l=='': l=chomp(pf.readline())#there's a (stupid) blank line in between entries
        while not l[0]=='>':
            ll=ll+l #glue all lines of the fasta entry
            l=chomp(pf.readline())
            if l=='': l=chomp(pf.readline())
            if not l: break

        #filter according to protein length and radius of gyration
        caf=cad+'/'+header+'.pdb' #CA file
        L=int(os.popen('grep " CA " '+caf+' | wc -l').readline().strip())
        if L<Ntarget*0.665 or L>Ntarget*1.333: continue #protein length not in target length
        rg,rg0=os.popen(rgx+' -a '+caf+' -b ').readline().split()
        rg=float(rg) ; rg0=float(rg0)
        if rg > 1.25 * rg0: continue #protein is not compact

        for x in alph2: #same routine for H,E,-
            iterator = px[x].finditer(ll)
            for match in iterator: hx[x][len(match.group(1))]+=1

        #output secondary structure content
        LL=len(ll)
        if LL>0:
            for x in Nx.keys(): Nx[x]=len(''.join(qx[x].findall(ll))) #number of H's, E's, -'s
            for x in Nx.keys(): Nx[x]=(1.0*Nx[x])/LL
            content='%4.2lf %4.2lf %4.2lf'%(Nx['H'],Nx['E'],1-Nx['H']-Nx['E'])
            if   Nx['E']<0.1 and Nx['H']>0.19: type='H'
            elif Nx['H']<0.1 and Nx['E']>0.19: type='E'
            elif Nx['H']>0.1 and Nx['E']>0.1:  type='HB'
            print type+' '+content+' %3d %4.1lf %4.1lf'%(L,rg,rg0)+' '+header

        #store number of helices and strands per sequence
        #print ll
        if LL>125 and LL<175:
            if type=='H':
                rx['nHinH'][len(px['(H+)'].findall(ll))]+=1 #number of helixes in this protein
            elif type=='E':
                rx['nEinE'][len(px['(E+)'].findall(ll))]+=1 #number of strands in this protein
            else:
                rx['nHinHE'][len(px['(H+)'].findall(ll))]+=1
                rx['nEinHE'][len(px['(E+)'].findall(ll))]+=1

pf.close()

#output the  histograms
alph3={'(H+)':'H', '(E+)':'E', '^(-+)':'Nt', '(-+)$':'Ct', 'H(-+)H':'H_H', 'E(-+)E':'E_E', 'H(-+)E':'H_E', 'E(-+)H':'E_H'}
for x in alph2: #same routine for H,E,-
    pf=open(alph3[x]+'.dis','w')
    for y in range(N):  pf.write('%3d %4d\n'%(y,hx[x][y]))
    pf.close()

for x in nXinY:
    pf=open(x+'.dis','w')
    for y in range(N):  pf.write('%3d %4d\n'%(y,rx[x][y]))
    pf.close()
    
sys.exit(0)
