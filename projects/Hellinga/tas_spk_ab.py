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
        '  -f _RA_outdir output directory',
        '  -g _A_tasserx tasser executable (def:codedir/f77/tasser/tasser.1.0/abinitio/cas_abinitio)',
        '  -i _A_inddf in.dd file (def:codedir/projects/Hellinga/brun/in.dd)',
        '  -j _A_CA native file (def:none)'
        ).parse(locals(),sys.argv)

if not tasserx: tasserx=codedir+'/f77/tasser/tasser.1.0/abinitio/cas_abinitio'
if not inddf: inddf=codedir+'/projects/Hellinga/brun/in.dd'

pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

root=codedir+'/f77/tasser/tasser.1.0'
pastry('/bin/cp '+inddf+' .')
pastry('/bin/tar jxf '+root+'/common.tbz2')
pastry('/bin/tar jxf '+dattar)
pastry('/bin/tar jxf '+pairtar)
pastry('/bin/tar jxf '+inprosp)
pastry('/bin/cp seq/'+header+'.SEQ seq.dat')
pastry('/bin/cp '+expdat+' .')
pastry('tar jcf '+header+'.in.ab_tasser.tar.bz2 exp.dat pa*dat seq.dat in.dd rmsinp')
pastry('/bin/mv '+header+'.in.ab_tasser.tar.bz2 '+outdir)

pastry('/bin/cp '+tasserx+' .')
tasserx=os.path.basename(tasserx)
pastry('./'+tasserx+' >/dev/null')
ntra=os.popen('ls -1 ./rep*tra|wc -l').readline().strip()
pastry('echo '+ntra+' > tra.in && ls -1 rep*.tra >> tra.in')
pastry('tar jcf '+header+'.in.ab_spicker.tar.bz2 rep*tra tra.in')
pastry('/bin/mv '+header+'.in.ab_spicker.tar.bz2 '+outdir)

spicker=codedir+'/f77/spicker/spicker0.0/spicker'
pastry('/bin/cp '+spicker+' .')
if CA: pastry('/bin/cp '+CA+' .')
pastry('./spicker')
pastry('tar jcf '+header+'.out.ab_spicker.tar.bz2 combo* closc* str.txt rst.dat')
pastry('/bin/mv '+header+'.out.ab_spicker.tar.bz2 '+outdir)
sys.exit(0)
