#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp,junkName,Bye

from seq.alignYangManageResults import alignOut

root='/gpfs1/scratch/jose/G8'

headers=chomp(open(root+'/input/list','r').readlines())
N=len(headers)
avseqid=0.0  ;  n=0
for i in range(0,N-1):
    hi=headers[i] #; print hi
    fastai=root+'/input/fastas/'+hi #; print fastai
    for j in range(i+1,N):
        hj=headers[j] #; print '  '+hj
        fastaj=root+'/input/fastas/'+hj
        si=float(alignOut(os.popen('$HOME/bin/align '+fastai+' '+fastaj+' 0').readlines()).seqid)
        print si
        avseqid+=si  ;  n+=1
print avseqid/n
sys.exit(0)
