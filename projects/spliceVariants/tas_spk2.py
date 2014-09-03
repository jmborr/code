#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir

from seq.fastaManager import importFastaEntry

inpHand('Usage: tas_spk.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_intar compressed tar file containing all neccessary input for tasser, except for common',
        '  -f _RA_outdir output directory',
        '  -g _A_cas cas version (def:codedir+/f77/tasser/tasser.1.0/original/cas',
        '  -i _A_spicker spicker version (def:codedir+/f77/spicker/spicker0.0/spicker)'
        ).parse(locals(),sys.argv)

pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

root=codedir+'/f77/tasser/tasser.1.0'
os.system('/bin/cp '+root+'/common.tbz2 . && tar jxf common.tbz2')

pastry('/bin/tar jxf '+intar)

if not cas: cas=codedir+'/f77/tasser/tasser.1.0/original/cas'
pastry('/bin/cp '+cas+' .')
pastry('./'+os.path.basename(cas))
pastry('echo `ls -1 rep*.tra|wc -l` > tra.in && ls -1 rep*.tra >> tra.in')
pastry('tar jcf '+header+'.in.spicker.tar.bz2 rep*tra tra.in seq.dat chain.dat rmsinp exp.dat')
pastry('/bin/mv '+header+'.in.spicker.tar.bz2 '+outdir)

if not spicker: spicker=codedir+'/f77/spicker/spicker0.0/spicker'
pastry('/bin/cp '+spicker+' .')
pastry('./'+os.path.basename(spicker))
if not os.system('ls loop*.dat &>/dev/null'):
    pastry('tar jcf '+header+'.out.spicker.tar.bz2 combo* closc* str.txt rst*.dat loop*.dat')
else:
    pastry('tar jcf '+header+'.out.spicker.tar.bz2 combo* closc* str.txt rst*.dat')
pastry('/bin/mv '+header+'.out.spicker.tar.bz2 '+outdir)

sys.exit(0)
