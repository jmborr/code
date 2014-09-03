#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.codedir import codeNodeLib

from seq.fastaManager import importFastaEntry

inpHand('Usage: prosp.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_inprosp tarred file with input for prospector',
        '  -c _RA_outdir output directory',
        '  -d __co remove homologs above sequence cutoff for prospector  (def:1.0)'
        ).parse(locals(),sys.argv)

if not co: co='1.0'

pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)

header,seq=importFastaEntry(open(fastaf,'r'))

if header!=fastaf: os.system('/bin/cp '+fastaf+' '+header)

taropts='xf'
if 'bz2' in inprosp: taropts='jxf'
pastry('/bin/tar '+taropts+' '+inprosp) #; pastry('ls')

prospector=codeNodeLib+'/python/prospector/prospector_genomesmay06/modified/prospector.py'
os.system(prospector+' -a '+header+' -b . -g '+co+' -i '+fastaf)
pastry('tar jcf '+header+'.out.prospector.tar.bz2 pdbb*')
pastry('/bin/mv '+header+'.out.prospector.tar.bz2  '+outdir)

sys.exit(0)
