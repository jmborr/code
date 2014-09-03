#!/usr/bin/python

"""For every reference sequence, find residues that participate at least
in one variant fragment. This will define the \"enclosing\" variant
fragment. The remaining residues will be the \"minimal\" constitutive
region."""

import os,sys
from copy import copy
from jobs.job import Bye
from seq.alignYangManageResults import alignOut
from seq.letters import valid_A,background_probs

refs={}
seqs={}
pidx=open('/gpfs1/scratch/jose/spliceVariants/preparing_database/var0_var2.alg.idx','r')
pin=open('/gpfs1/scratch/jose/spliceVariants/preparing_database/var0_var2.alg','r')

line=pidx.readline()
while line:    
    pos,ref,var=line.split()    
    pin.seek(int(pos))
    alg=alignOut(pin).alg #read one alignment
    ll=alg.isAlignedLines()[1]
    if ref not in seqs.keys(): seqs[ref]=alg.gs[1].s #ungapped sequence of reference
    if ref not in refs.keys(): refs[ref]=ll
    else:
        l=refs[ref]
        for i in range( len(ll) ):
            if ll[i]=='0': l=l[:i]+'0'+l[i+1:]
    line=pidx.readline()

initcomp={}
for X in valid_A: initcomp[X]=0
comp={'1':copy(initcomp),'0':copy(initcomp)} #compositions for constitutive ('1') and variant ('0')
n={'1':0,'0':0}
for ref in refs.keys():
    seq=seqs[ref] #reference sequence
    l=refs[ref]   #string of '1' and '0' if residue is constitutive or variant
    for i in range(len(l)):
        comp[ l[i] ][ seq[i] ]+=1
        n[ l[i] ]+=1

#normalize compositions
for type in comp.keys():
    for X in valid_A:
        comp[type][X]=float(comp[type][X])/n[type]

#print compositions
print '#AA PDB constitutive variant'
for x in valid_A:
    print x,'%5.3lf %5.3lf %5.3lf'%(background_probs[x],comp['1'][x],comp['0'][x])
        

sys.exit(0)

