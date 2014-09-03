#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq
from jobs.job import job

ih=inpHand('Usage: sleepy_list.py',
           ' -l _AR_list (list of five-letter codes)',
           ' -u __submode submission mode of each job (sub,bsub,rsub) def=bsub',
           ' -h __help this script runs (possibly) prospector, tasser and spicker on a list of five-letter codes'
           )
ih.parse(locals(),sys.argv)
#if not submode: submode='bsub'
if not submode: submode='info'

#instantiate a sleepy job
j=job(exe='sleepy.py', exed='/gpfs1/active/jose/code/python/combo_jobs', args='' )

#bsub one job at a time, change j.args for each job
for header in open(list,'r').readlines():
    j.args=''
    #getattr returns reference to object method (it is also an attribute of
    #jobs j) associated to string submode. Then () obligues execution
    getattr(j,submode)() 

