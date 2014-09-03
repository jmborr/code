#!/usr/bin/python

import os,sys
from utilities.small_utilities import chomp,Bye

"""
Find percentage of pairs (reference,isoform) having the same first template

We use file '/gpfs1/scratch/jose/spliceVariants/out/best_prospector_templ.dat'
"""

#initialize templs, which stores templates for each variant in a dictionary
templs={}
for line in os.popen('cat /gpfs1/scratch/jose/spliceVariants/out/best_prospector_templ.dat|grep -v "#"','r').readlines():
    variant,templ,Z=line.split()
    root=variant[0:6]
    isoform=variant[6:8]
    if root not in templs.keys(): templs[root]={}
    if isoform not in templs[root].keys(): templs[root][isoform]=templ

npairs=0
nequal=0
buf=''
for root in templs.keys():
    dict=templs[root]
    isoforms=dict.keys() ; isoforms.sort()
    if '00' not in isoforms: continue
    templ0=dict['00']
    for isoform in isoforms[1:]: #the first one is '00', the reference
        buf+=root+'00 '+root+isoform+'\n'
        npairs+=1
        if templ0==dict[isoform]: nequal+=1

open('/gpfs1/scratch/jose/spliceVariants/out/ref_var_pairs.dat','w').write(buf)
print float(nequal)/npairs

sys.exit(0)

