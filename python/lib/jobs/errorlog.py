#!/usr/bin/python
import os,re,sys
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp #remove '\n' of strings and lists

inpHand('Usage: error.py',
        ' -a _AR_rootD root directory. Subdirectories contain the log files',
        ' -b _R_q select which logs to output (one of COMPLETED,ERROR,ALL)',
        ' -c _A_outf output file (default=standar output)',
        ' -h __help this script will output the log files under a root directory. Will output all logfile, or only those reporting error, or those that report a completed job'
        ).parse(locals(),sys.argv)

if q not in ('COMPLETED','ERROR','ALL'):
    sys.stderr.write('ERROR: the selection must be one of "COMPLETED","ERROR","ALL"\n')
    sys.exit(1)

out=sys.stdout
if outf: out=open(outf,"w")

os.chdir(rootD)
logs=os.popen('find . -name "*.log"').readlines() #print logs
chomp(logs)

for log in logs:
    if not os.system('grep "exit mode = 0" '+log+
                     ' >>/dev/null 2>>/dev/null' ):#log reports a completed job
        if q=='COMPLETED' or q=='ALL':
            out.write(log+'\n')
    else:#logfile reports an unsuccessful job
        if q=='ERROR' or q=='ALL':
            out.write(log+'\n')

print 'errolog.py finished'
sys.exit(0)
