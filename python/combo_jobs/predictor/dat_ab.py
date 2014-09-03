#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand

from seq.fastaManager import importFastaEntry

inpHand('Usage: dat_ab.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_outdir output directory',
        ).parse(locals(),sys.argv)

pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

pastry('mkrmsinp.py -i '+header+' -o . -f '+fastaf)
pastry('/bin/tar cf '+header+'.dat.tar rmsinp')
pastry('/bin/mv '+header+'.dat.tar '+outdir)

sys.exit(0)
