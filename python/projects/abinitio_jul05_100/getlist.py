#!/usr/bin/python
import os,sys,re
from utilities.codedir import codedir,libseq
from Bio.PDB.PDBParser import PDBParser

N=100 #limit of too long sequences
#directory with list of homologous for each file
nonHomD='/net/dell/01/users/skolnick/pdbhomoljul05'

#retrieve all headers
headers=map( (lambda x:x[:-1]),
             os.popen('ls -1 '+nonHomD+' | grep -e ".homol$" | cut -d \'.\' -f 1').readlines()
             )

#filter those with less than 100 amino acids
for header in headers:
    P=[]
    chId=header[-1]
    if chId=='_': chId=' '
    pdbf='/library/pdb/pdb'+header[0:4]+'.ent'
    if not os.path.exists(pdbf): ih.abort(' !!!Non existent PDB file '+pdbf)
    n=PDBParser(PERMISSIVE=1).get_structure(header[0:4],pdbf)[0][chId].__len__()
    print header,n
    if n<N: P.append(header)

#make a list with homologous headers to the one passed
def nonhom(hd):
    return map( (lambda x:x[:-1]),
                os.popen('cat '+nonHomD+'/'+hd+'.homol | cut -d \' \' -f 1').readlines()
             )
    
#filter those homologous
i=0
while i<len(P):
    P2=P[0:i+1][:] #copy constructor
    Q=nonhom(P[i])
    for p in P[i+1:]:
        if not p in Q:
            P2.append(p)
    P=P2[:]
    i=i+1

for p in P: print p
