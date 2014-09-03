#!/usr/bin/python

import os,sys
from utilities.small_utilities import chomp,Bye

def rgfun(L):
    return 1.20*2.2*L**0.38
    
rootd='/gpfs1/scratch/jose/Hellinga'
    
for header in chomp(open(rootd+'/pdb_aug06.apr06.CA/list','r').readlines()):
    pdbf=rootd+'/pdb_aug06.apr06.CA/'+header+'.pdb'
    cmd='grep "ATOM" '+pdbf+' | grep " CA "|wc -l'
    L=int( os.popen(cmd,'r').readline().split()[0] )
    rgmax=rgfun(L)
    cmd='Rg.x -a '+pdbf
    rg=float(os.popen(cmd,'r').readline())
    if L>=100 and L<200 and rg<rgmax: print header+' %3d %5.2lf %5.2lf'%(L,rg,rgmax)

sys.exit(0)
