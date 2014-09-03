#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand
from jobs.job import gen_lf,pastry

ih=inpHand('Usage: generic_job.py',
           ' -a __inpcmd input command line (in double quotes)',
           ' -b _R_execcmd executable command line, run from the temp dir (in double quotes)',
           ' -c __outcmd output command line  (in double quotes)',
           ' -d _RA_outd output directory where the logfile will go',
           ' -x _R_identifier identifier string to generate the logfile',
           ' -y  __logrequest either _ASK_, or _LOG_. _ASK_ and _LOG_ cause to automatically generate a logfile name with the identifier (flag -x), but _ASK_ will only output the logfile to STDOUT and exit(0). logfile will reside in the output directory',)
ih.parse(locals(),sys.argv)
logfile=gen_lf(logrequest,outd,identifier) #generate a suitable logfile

#create first output directory where logfile go
if not os.path.exists(outd): pastry('mkdir -p '+outd)
#bring input
if inpcmd: pastry(inpcmd)
#execute job
pastry(execcmd)
#move output
if outcmd: pastry(outcmd)

sys.exit(0)
