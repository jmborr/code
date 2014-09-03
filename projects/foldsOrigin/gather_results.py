#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp

re1=re.compile('(\d+) \(TM-score') #retrieve header's length
re2=re.compile('Aligned length=\s*(\d+)') #retrieve coverage
re3=re.compile('RMSD=\s*(\d+\.\d+)') #retrieve RMSD
re4=re.compile('TM-score=(0.\d+)') #retrieve TM-score

print '#header L=length l=aligned-length coverage=100*l/L RMSD-aligned TM-score'

root='/gpfs1/scratch/jose/foldsOrigin/pdbIn600x10'
headers=chomp(open(root+'/success.list','r').readlines())
for header in headers:
    tmfile=root+'/'+header[1]+'/'+header+'/'+header+'.tm'
    all=''.join(open(tmfile,'r').readlines()) #file in one long string
    L=int(re1.search(all).group(1))
    l=int(re2.search(all).group(1))
    cov=(100*l)/L
    rmsd=float(re3.search(all).group(1))
    tm=float(re4.search(all).group(1))
    print header+' %3d %3d %3d %5.2lf %5.3lf'%(L,l,cov,rmsd,tm)
