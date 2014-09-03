#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq
from utilities.small_utilities import chomp
from jobs.job import job

ih=inpHand('Usage: create_input_for_tasser_list.py',
           ' -l _AR_list (list of five-letter codes)',
           ' -c _R_createargs arguments to pass to create_input_for_tasser.py, enclosed whithin double quotes. Do not pass the -i flag (the header flag) and the -x flag (identifier flag). These flags will be set automatically',
           ' -u __submode submission mode of each job (sub,bsub,rsub,qsub) def=bsub',
           ' -q __queue submission queue(default=standard)',
           ' -b __memlimit memory reserved, in megabyes (def=850)',
           ' -v __xflags some possible extra flags to the submission command. Enclose them within double quotes',
           ' -x _R_identifier one word for job-name purposes and for logfile name. Will append xxxxx',
           ' -h __help this script generates input files for future tasser runs, for a list of five-letter pdb codes'
           )
ih.parse(locals(),sys.argv)

if not submode: submode='bqsub'
if not queue: queue='batch'
if not memlimit or memlimit>1000: memlimit=850

#instantiate a template job
j=job(exe='create_input_for_tasser.py', exed='/gpfs1/active/jose/code/python/combo_jobs', args='')

#check the list of headers
if not os.path.exists(list):
    stderr.write('ERROR: no '+os.basename(list)+' found in '+os.dirname(list))
    stdout.write('program exiting')
    sys.exit(1)

#bsub one job at a time, change j.args for each job
for header in open(list,'r').readlines():
    header=chomp(header)
    j.args=createargs+' -i '+header+' -x '+identifier+header
    print j.args
    sys.exit(1)
    #getattr returns reference to object method (it is also an attribute of
    #jobs j) associated to string submode. Then, () obligues execution
    getattr(j,submode)(queue=queue,mem_limit=memlimit,extraflags=xflags)

sys.exit(0)
