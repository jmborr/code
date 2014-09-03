#!/usr/bin/python

"""For every reference sequence, find residues that participate at
least in one variant fragment smaller than some particular size. This
will define the \"enclosing\" variant fragment. The remaining residues
will be the \"minimal\" constitutive region."""

import os,sys,re
from copy import copy
from jobs.job import Bye
from seq.alignYangManageResults import alignOut
from seq.letters import valid_A,background_probs
from inputArgs.inputArgs import inpHand

ih=inpHand('Usage: secondaryIdentityOverConstitutive.py [options]',
           '  -a __maxL maximum length of variant region (def: 20)',
           )
ih.parse(locals(),sys.argv)

if not maxL: maxL=20
else: maxL=int(maxL)
minL=4

refs={}
seqs={}
pidx=open('/gpfs1/scratch/jose/spliceVariants/preparing_database/var0_var2.alg.idx','r')
pin=open('/gpfs1/scratch/jose/spliceVariants/preparing_database/var0_var2.alg','r')
rex=re.compile('0+') #iterator over stretches of zeroes
line=pidx.readline()
while line:    
    pos,ref,var=line.split()    
    pin.seek(int(pos))
    alg=alignOut(pin).alg #read one alignment
    ll=alg.isAlignedLines()[1]
    ll=ll.replace('010','000')
    ll=ll.replace('101','111')
    ll=ll.replace('001100','000000')
    if ref not in seqs.keys(): seqs[ref]=alg.gs[1].s #ungapped sequence of reference
    if ref not in refs.keys(): refs[ref]='1'*len(ll) #initialize whole sequence as constitutive
    matches=rex.finditer(ll)
    for match in matches:
        (a,b)=match.span()
        if b-a<maxL and b-a>minL:
            l=refs[ref]
            refs[ref]=l[:a]+match.group()+l[b:]
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

