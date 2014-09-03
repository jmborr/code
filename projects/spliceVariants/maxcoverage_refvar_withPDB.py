#!/usr/bin/python

import os,sys,re
from utilities.codedir import scratchdir

pat=re.compile('id=(\d\.\d+)')
pin=open( os.path.join(scratchdir,'spliceVariants/preparing_database/coverage_of_refvars_withPDB.dat'),'r')
l=pin.readline()
while l and l[0:2]=='AC':
    buf=l.strip()
    l=pin.readline()
    maxcov=0.0
    id=0.0
    while l and l[0:2]=='DR':
        cov=float(l[-4:]) #;print cov;sys.exit(1)
        if cov > maxcov:
            maxcov=cov
            id=float(pat.search(l).group(1))
        l=pin.readline()
    print buf,'%4.2lf %4.2lf'%(id,maxcov)

sys.exit(0)
