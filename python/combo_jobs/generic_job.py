#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand
from jobs.job import pastry

#assumed outcmd='/bin/mv file1 file2 ... outdir'
def getOutDir(outcmd):
    return os.path.dirname(outcmd.split()[-1])

ih=inpHand('Usage: generic_job.py',
           ' -a __inpcmd input command line (in double quotes)',
           ' -b _R_execcmd executable command line, run from the temp dir (in double quotes)',
           ' -c __outcmd output command line  (in double quotes)',
           )
ih.parse(locals(),sys.argv)

if outcmd:
    outdir=getOutDir(outcmd)
    if not os.path.exists(outdir): pastry('mkdir -p '+outdir)

#bring input
if inpcmd: pastry(inpcmd)
#execute job
pastry(execcmd)
#move output
if outcmd: pastry(outcmd)

sys.exit(0)
