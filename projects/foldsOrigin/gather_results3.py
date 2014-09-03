#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp

re1=re.compile('Chain 1.*Size=\s*(\d+)')  #retrieve PDB-hit's length
re2=re.compile('Aligned length=\s*(\d+)') #retrieve coverage
re3=re.compile('RMSD=\s*(\d+\.\d+)')      #retrieve RMSD
re4=re.compile('TM-score=(0.\d+)')        #retrieve TM-score

print '#header PDB-hit-length aligned-length coverage=100*l/300 RMSD-aligned TM-score'

root='/gpfs1/scratch/jose/foldsOrigin/600x10inPDB/c'
headers=chomp(open(root+'/tms.list','r').readlines())
for header in headers:
    tmfile=root+'/'+header+'.tm'
    all=''.join(open(tmfile,'r').readlines()) #file in one long string
    L=int(re1.search(all).group(1))
    l=int(re2.search(all).group(1))
    cov=(100*l)/300
    rmsd=float(re3.search(all).group(1))
    tm=float(re4.search(all).group(1))
    print '%-26s %4d %3d %3d %5.2lf %5.3lf'%(header,L,l,cov,rmsd,tm)
