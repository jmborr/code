#!/usr/bin/python

import os,sys

n=0
rootd='/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/out'
cmd='for h in `cat '+rootd+'/list`;do echo -n "$h " && grep -v " 1 " $h/seq.dat|wc -l;done > '+rootd+'/'+junk
os.system(cmd)
lines=open(rootd+'/junk','r').readlines()
for line in lines:
    id,m=line.split()
    if int(m)<65:
        cmd='/bin/mv '+rootd+'/'+id+' '+rootd+'/d%04d'%(n) ;print cmd
        n+=1
        os.system(cmd)
sys.exit(0)
