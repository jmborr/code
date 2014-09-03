#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir

from seq.fastaManager import importFastaEntry

inpHand('Usage: tas_spk.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_dattar contains rmsinp, comb* dist* chain.dat',
        '  -c _RA_pairtar contains pair1.dat, pair3.dat',
        '  -d _RA_inprosp contains seq.dat',
        '  -e _RA_expdat exp.dat file',
        '  -f _RA_outdir output directory',
	'  -g _A_inddf in.dd file (def:none)'
        ).parse(locals(),sys.argv)

pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

root=codedir+'/f77/tasser/tasser.1.0'

if inddf and os.path.exists(inddf):
    os.system('/bin/cp '+inddf+' . ')
else:
    os.system('/bin/cp '+root+'/original/in.dd . ')
    
os.system('/bin/cp '+root+'/common.tbz2 . && tar jxf common.tbz2')

pastry('/bin/tar xf '+dattar)
pastry('/bin/tar xf '+pairtar)
pastry('/bin/tar xf '+inprosp) ;pastry('/bin/cp seq/'+header+'.SEQ seq.dat')
pastry('/bin/cp '+expdat+' .')
pastry('tar cf '+header+'.in.tasser.tar comb*dat chain.dat dist*dat exp.dat pa*dat seq.dat in.dd rmsinp')
pastry('/bin/mv '+header+'.in.tasser.tar '+outdir)

cas=codedir+'/f77/tasser/tasser.1.0/original/cas'
pastry('/bin/cp '+cas+' .')
pastry('./cas')

pastry('echo `ls -1 rep*.tra|wc -l` > tra.in && ls -1 rep*.tra >> tra.in')
pastry('tar cf '+header+'.in.spicker.tar rep*tra tra.in')
pastry('/bin/mv '+header+'.in.spicker.tar '+outdir)

spicker=codedir+'/f77/spicker/spicker0.0/spicker'
pastry('/bin/cp '+spicker+' .')
pastry('./spicker')
pastry('tar cf '+header+'.out.spicker.tar combo* closc* str.txt rst.dat')
pastry('/bin/mv '+header+'.out.spicker.tar '+outdir)

sys.exit(0)
