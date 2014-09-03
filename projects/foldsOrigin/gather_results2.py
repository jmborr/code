#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp

file='/gpfs1/scratch/jose/foldsOrigin/600x10inPDB/a/tm_cov_lengthHit_lx.list'

pt=open(file,'r')

pt.readline() #first line is a comment
print '#header tm #aligned-residues coverage(%) PDB-hit-length extended-alignment-length'

for line in chomp(pt.readlines()):
    items=line.split() #split, separator is blank space
    header=items[0]
    tm=float(items[1])
    l=int(items[3])        #number of aligned residues
    l0=int(items[5])       #length of the PDB hit
    cov=int( (100*l)/300 ) #coverage of the alignment. Remember template length is 300
    lx=int(float(items[8]))       #extended alignment length
    print '%-9s %5.3lf %3d %3d %4d %4d'%(header,tm,l,cov,l0,lx)
