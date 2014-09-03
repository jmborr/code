#!/usr/bin/python

import os,sys,re
from utilities.small_utilities import chomp,Bye

    
rootd='/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/out'

for header in chomp(open(rootd+'/a.list','r').readlines()):
    seq=''
    seqdatf=rootd+'/'+header+'/seq.dat'
    lines=open(seqdatf,'r').readlines()
    for line in lines:
        seq+=line.split()[2]
    print len(re.compile('2+').findall(seq))

sys.exit(0)
