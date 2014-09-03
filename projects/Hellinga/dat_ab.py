#!/usr/bin/python

import os,sys
from jobs.job import pastry
from utilities.small_utilities import Bye
from inputArgs.inputArgs import inpHand

from seq.fastaManager import importFastaEntry

inpHand('Usage: dat_ab.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_outdir output directory',
        '  -c _RA_seqdatf seq.dat file',
        '  -d _A_chaindatf chain.dat file (def:none)'
        ).parse(locals(),sys.argv)

if chaindatf and not os.path.exists(chaindatf):
    Bye('ERROR from dat_ab.py: file '+chaindatf+' does not exists')
    
pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

#create rmsinp file
pastry('/library/jose/code/python/tasser/tasser1.0/mkrmsinp.py -i '+header+' -o . -f '+fastaf)

#create chain.dat file
pastry('/bin/cp /gpfs1/active/jose/code/projects/Hellinga/genChainDat.py /gpfs1/active/jose/code/projects/Hellinga/betahelix.pdb .')

if not chaindatf:
    pastry('./genChainDat.py -a '+seqdatf+' -b betahelix.pdb -c chain.dat -d init.pdb -f no')
    pastry('/bin/tar jcf '+header+'.dat_ab.tar.bz2 rmsinp chain.dat init.pdb')
else:
    pastry('/bin/cp '+chaindatf+' ./chain.dat')
    pastry('/bin/tar jcf '+header+'.dat_ab.tar.bz2 rmsinp chain.dat')

pastry('/bin/mv '+header+'.dat_ab.tar.bz2 '+outdir)

sys.exit(0)
