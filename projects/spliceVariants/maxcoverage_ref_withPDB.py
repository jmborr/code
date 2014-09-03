#!/usr/bin/python

import os,sys
from utilities.codedir import scratchdir

pin=open( os.path.join(scratchdir,'spliceVariants/preparing_database/coverage_of_refs_withPDB.dat'),'r')
l=pin.readline()
while l and l[0:2]=='AC':
    buf=l.strip()
    l=pin.readline()
    maxcov=0.0
    while l and l[0:2]=='DR':
        cov=float(l[-4:]) #;print cov;sys.exit(1)
        if cov > maxcov: maxcov=cov
        l=pin.readline()
    print buf,'%4.2lf'%maxcov

sys.exit(0)
