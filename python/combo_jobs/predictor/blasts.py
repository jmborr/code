#!/usr/bin/python

import os,sys
from time import sleep
from seq.blastManager import blastRun
from seq.fastaManager import importFastaEntry
from tasser.mksolve import gen_solve
from inputArgs.inputArgs import inpHand
from jobs.job import pastry

inpHand('Usage: blasts.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_outdir output directory'
        ).parse(locals(),sys.argv)

#run blast
pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))
blastOuts=blastRun(fastaf,header=header)
if blastOuts['errors']['blastpgp']: sys.exit(1) #check blast job did run
pastry('tar jcf '+header+'.out.blast.tar.bz2 *.chk *.fasta *.blast')
pastry('/bin/mv '+header+'.out.blast.tar.bz2 '+outdir)

#run blast for the exposed surface predictor
gen_solve(fastaf) #create exp.dat
pastry('/bin/mv exp.dat '+outdir)

sleep(10) ; sys.exit(0)
