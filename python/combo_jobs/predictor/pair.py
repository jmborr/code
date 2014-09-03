#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand

from seq.fastaManager import importFastaEntry
from seq.msa import gen_msa
from tasser.mkpair import gen_pair

inpHand('Usage: pair.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_blastout tarred file with previously run blast',
        '  -c _RA_inprosp tarred file with input for prospector',
        '  -d __co sequence cut-off (def:0.3)',
        '  -e _RA_outdir output directory',
        '  -f __dotar do we tar the results? (default:yes)'
        ).parse(locals(),sys.argv)

if co:co=float(co)
else: co=0.3
if dotar not in ('no','n','N','NO'): dotar=1

pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

taropts='xf'
if 'bz2' in inprosp: taropts='jxf'
pastry('/bin/tar '+taropts+' '+inprosp)

taropts='xf'
if 'bz2' in blastout: taropts='jxf'
pastry('/bin/tar '+taropts+' '+blastout) #below we assume certain naming convention for the files
blastOuts={'outf':'', 'blast':header+'.blast', 'chk':header+'.chk', 'fasta':fastaf}

msas=gen_msa(seq,header,Eco=10.0  ,maxId=0.90,minId=0.10,red=0.95,blastOuts=blastOuts)
gen_pair(msas, 'seq/'+header+'.SEQ', '.', co=co)

if dotar==1:
    pastry('/bin/tar jcf '+header+'.pair.tar.bz2 pair1.dat pair3.dat')
    pastry('/bin/mv '+header+'.pair.tar.bz2 '+outdir)
else:
    pastry('/bin/mv pair1.dat pair3.dat '+outdir)
    
sys.exit(0)
