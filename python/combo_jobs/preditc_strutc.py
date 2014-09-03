#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq
from jobs.job import gen_lf,pastry

ih=inpHand('Usage: preditc_struct.py',
           ' -a __iidp run ID, if submitting more than one run per target (def=None)',
           ' -m __move list of files to move  (in double quotes, globbing allowed)',
           ' -o _A_outd directory to save output (default is current dir)',
           ' -p _A_prospector prospector job (in double quotes)',
           ' -t _A_tasser tasser job    (in double quotes)',
           ' -x _R_identifier identifier string to generate the logfile',
           ' -y __logrequest either _ASK_, or _LOG_ (def=_LOG_). _ASK_ and _LOG_ cause to automatically generate a logfile name with the identifier (flag -x), but _ASK_ will only output the logfile to STDOUT and exit. logfile will reside in the output directory',
           ' -z _A_spicker spicker job  (in double quotes)'
           )
ih.parse(locals(),sys.argv)
if not outd: outd=os.getcwd();
print 'identifier=',identifier
logfile=gen_lf(logrequest,outd,identifier) #generate a suitable logfile

#create output directory before running the job. If job is killed, at
#least the logfile will be stored in the directory
pastry('/bin/mkdir -p '+outd)

if prospector:
    print '...predict_struct.py executing prospector job...'
    pastry(prosp)
if tasser:
    print '...predict_struct.py executing tasser job...'
    #add run ID if we provided one
    if iidp: tasser=tasser+' -j '+iidp
    pastry(tasser)
if spicker:
    print '...predict_struct.py executing spicker job...'
    pastry(spicker)

pastry('/bin/mv -f '+move+' '+outd)
sys.exit(0)
