#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir

from seq.fastaManager import importFastaEntry

inpHand('Usage: tas_spk.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_dattar tarred file from dat_ab.py containing rmsinp',
        '  -c _RA_pairtar tarred file from pair.py containing pair1.dat, pair3.dat',
        '  -d _RA_inprosp tarred file from inprosp_ab.py containing seq.dat',
        '  -e _RA_expdat exp.dat file',
        '  -f _RA_outdir output directory'
        ).parse(locals(),sys.argv)

pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

root=codedir+'/f77/tasser/tasser.1.0'
os.system('/bin/cp '+root+'/abinitio/in.dd .')
pastry('/bin/tar jxf '+root+'/common.tbz2')
pastry('/bin/tar xf '+dattar)
pastry('/bin/tar xf '+pairtar)
pastry('/bin/tar xf '+inprosp)
pastry('/bin/cp seq/'+header+'.SEQ seq.dat')
pastry('/bin/cp '+expdat+' .')
pastry('tar cf '+header+'.in.tasser.tar exp.dat pa*dat seq.dat in.dd rmsinp')
pastry('/bin/mv '+header+'.in.tasser.tar '+outdir)

pastry('/bin/cp '+codedir+'/f77/tasser/tasser.1.0/abinitio/cas_abinitio .')
pastry('./cas_abinitio >/dev/null')
ntra=os.popen('ls -1 ./rep*tra|wc -l').readline().strip()
pastry('echo '+ntra+' > tra.in && ls -1 rep*.tra >> tra.in')
pastry('tar cf '+header+'.in.spicker.tar rep*tra tra.in')
pastry('tar cf '+header+'.in.spicker.tar rep*tra tra.in')
pastry('/bin/mv '+header+'.in.spicker.tar '+outdir)

spicker=codedir+'/f77/spicker/spicker0.0/spicker'
pastry('/bin/cp '+spicker+' .')
pastry('./spicker')
pastry('tar cf '+header+'.out.spicker.tar combo* closc* str.txt rst.dat')
pastry('/bin/mv '+header+'.out.spicker.tar '+outdir)

sys.exit(0)
