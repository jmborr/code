#!/usr/bin/python
'''
for all pairs of alignments in refsAboveC1, create superposition of
the structures with TM score for the aligned residues according to the
global sequence alignment passed in file alfg
'''
import os,sys
from utilities.small_utilities import chomp,Bye

root='/gpfs1/scratch/jose/spliceVariants/wrong/refsAboveC1'
root2='/gpfs1/scratch/jose/spliceVariants/wrong/tasser'

listf=root+'/pairs.list'
#listf=root+'/toy.list'
pairs=chomp(open(listf,'r').readlines())

for pair in pairs:
    p1,p2=pair.split() ; p=p1[:-1] ; s1=p1[-1]; s2=p2[-1]
    algf=root+'/'+p+'/'+p+'_'+s1+'_'+s2+'.alignment.dat' #alignment file
    pdb1=root2+'/'+p[1]+'/'+p1+'/combo1.pdb' #assume combo1 is the top cluster
    pdb2=root2+'/'+p[1]+'/'+p2+'/combo1.pdb' #assume combo1 is the top cluster

    if not os.path.exists(algf) or not os.path.exists(pdb1) or not os.path.exists(pdb2): continue

    sup=root+'/'+p+'/'+p+'_'+s1+'_'+s2+'.alignment.sup.pdb' #output superpostion file
    line=os.popen('tm_score.py -a '+pdb1+' -b '+pdb2+' -c '+sup+' -i '+algf+' -d yes').readline()
    sys.stdout.write(pair+' '+line)
sys.exit(0)
