#!/usr/bin/python

"""Regularly checks size of a file. Outputs and kills a job if size
exceeded some limit"""

import os,sys
from time import sleep
from utilities.small_utilities import Bye
from inputArgs.inputArgs import inpHand,addAbsPath


#global variables
interval0=10
inpHand('Usage: contrologsize.py [options]',
        ' -a _R_jobid jobid',
        ' -b _RA_logf  log file',
        ' -c __maxsize maximum log file size, in Kb (def:100)',
        ' -d __interval check log file every certain number of seconds (def:'+`interval0`+')',
        ).parse(locals(),sys.argv)

#parse input arguments
if not interval: interval=interval0
else: interval=min(interval0, int(interval) )
if not maxsize: maxsize=100000
else: maxsize=min(100000,int(maxsize)*1000)

while True:
    sleep(interval)
    size=int( os.popen('du -b '+logf).readline().split()[0] )
    
    if size>maxsize:
        sys.stderr.write('ERROR in contrologsize: exceeded allowed log file size\n')
        os.system('kill -9 '+jobid)
        sys.exit(1)

sys.exit(0)
