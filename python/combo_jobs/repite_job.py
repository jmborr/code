#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq
from jobs.job import job

ih=inpHand('Usage: repite_job.py',
           ' -x _R_exec executable',
           ' -d _AR_execd directory where to find the executable',
           ' -a __arguments arguments to pass the executable (enclose in double quotes if necessary',
           ' -n _R_nr number of repetitions of the job',
           ' -u __submode submission mode of each job (sub,bsub,rsub) def=bsub',
           ' -v __xflags bsub may need these extra flags. Enclose them within double quotes',
           ' -h __help this script runs (possibly) prospector, tasser and spicker on a list of five-letter codes'
           )
ih.parse(locals(),sys.argv)

if not submode: submode='rsub'
j=job(exe=exec, exed=execd, args=arguments)

for i in range(nr):
    getattr(j,submode)(extraflags=xflags)


