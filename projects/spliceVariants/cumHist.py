#!/usr/bin/python

""" calculate the cumulative histogram of (reference,isoform) pairs
depending on TM(templ_ref,temp_iso) and min(Z(templ_ref),Z(temp_iso))"""

import os,sys
from jobs.job import Bye
from math import log10

cumhist=[] #cumhist[Z][tm]
for i in range(71): cumhist.append([1]*101)

cmd='cat /gpfs1/scratch/jose/spliceVariants/out/ref_var_tms2.dat|grep -v -e "#"'
for line in os.popen(cmd).readlines():
    it=line.split()
    tm=it[6]
    if tm!='inf':
        tm=1+int(float(tm)*100) #;print tm
    else:
        continue
    Z=1+int(float(it[7]))
    for z in range(Z):
        hist=cumhist[z]
        for t in range(tm): hist[t]+=1

for z in range(0,71):
    hist=cumhist[z]
    for t in range(101):
        print '%2d %4.2f %d %f'%(z,0.01*t,hist[t],log10(hist[t]))
        

sys.exit(0)
