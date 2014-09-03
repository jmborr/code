#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand

from seq.fastaManager import importFastaEntry

inpHand('Usage: dat.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_blastout tarred file with previously run blast',
        '  -c _RA_inprosp tarred file with input for prospector',
        '  -d _RA_prosptar tarred file with output for prospector',
        '  -e _RA_outdir output directory',
        '  -f __co sequence cut-off (def:0.3)'
        ).parse(locals(),sys.argv)

if co:co=float(co)
else: co=0.3

pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

pastry('mkrmsinp.py -i '+header+' -o . -f '+fastaf)

taropts='xf'
if 'bz2' in inprosp: taropts='jxf'
pastry('/bin/tar '+taropts+' '+inprosp)

taropts='xf'
if 'bz2' in prosptar: taropts='jxf'
pastry('/bin/tar '+taropts+' '+prosptar)

taropts='xf'
if 'bz2' in blastout: taropts='jxf'
pastry('/bin/tar '+taropts+' '+blastout) #below we assume certain naming convention for the files
blastOuts={'outf':'', 'blast':header+'.blast', 'chk':header+'.chk', 'fasta':fastaf}

pastry('mkdat.py -a '+header+' -b '+fastaf+' -c seq/'+header+'.SEQ -e pdbbpdborienrev/'+header+'rap3orienrev5s.pdb -g pdbbzpot4a/'+header+'.predictedrap3orienrev -l '+blastOuts['blast']+' -m '+`co`)

pastry('/bin/tar jcf '+header+'.dat.tar.bz2 rmsinp comb*dat chain.dat dist*dat par.dat')
pastry('/bin/mv '+header+'.dat.tar.bz2  '+outdir)

sys.exit(0)
