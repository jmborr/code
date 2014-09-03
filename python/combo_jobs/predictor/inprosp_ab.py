#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand

from seq.fastaManager import importFastaEntry
from seq.msa import gen_msa,output_msa
from seq.mtx import gen_mtx
from secondary.psipred.runpsipred import run_psipred
from secondary.mkseq import output_sec

inpHand('Usage: inprosp_ab.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_blastout tarred file with previously run blast',
        '  -c _RA_outdir output directory'
        ).parse(locals(),sys.argv)

msas={}
pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

pastry('/bin/tar xf '+blastout) #below we assume certain naming convention for the files
blastOuts={'outf':'', 'blast':header+'.blast', 'chk':header+'.chk', 'pssm':header+'.pssm', 'fasta':fastaf}

pastry('/bin/mkdir -p profile seq')

gen_mtx(seq,'profile/'+header+'.mtx',header=header,blastOuts=blastOuts) #profile

prediction=run_psipred(fastaf,mtxf='profile/'+header+'.mtx')#run psipred with generated profile
output_sec(seq,'seq/'+header+'.SEQ',prediction=prediction)

pastry('tar cf '+header+'.in.prospector.tar profile seq')
pastry('/bin/mv '+header+'.in.prospector.tar '+outdir)

sys.exit(0)
