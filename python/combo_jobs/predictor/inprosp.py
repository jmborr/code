#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand

from utilities.small_utilities import abort
from seq.fastaManager import importFastaEntry
from seq.msa import gen_msa,output_msa
from seq.mtx import gen_mtx
from secondary.psipred.runpsipred import run_psipred
from secondary.mkseq import output_sec

inpHand('Usage: inprosp.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _RA_blastout tarred file with previously run blast',
        '  -c _RA_outdir output directory'
        ).parse(locals(),sys.argv)

msas={}
if not os.path.exists('./'+os.path.basename(fastaf)): pastry('/bin/cp '+fastaf+' .')
fastaf=os.path.basename(fastaf)
header,seq=importFastaEntry(open(fastaf,'r'))

taropts='xf'
if '.bz2' in blastout or 'tbz2' in blastout: taropts='jxf'
if '.gz' in blastout or 'tgz' in blastout: taropts='zxf'
pastry('/bin/tar '+taropts+' '+blastout) #below we assume certain naming convention for the files
blastOuts={'outf':'', 'blast':header+'.blast', 'chk':header+'.chk', 'fasta':fastaf}

pastry('/bin/mkdir -p 3590 e10 profile seq')

msas['3590']=gen_msa(seq,header,Eco=0.001,maxId=0.90,minId=0.35,red=0.95,blastOuts=blastOuts)
if not output_msa(msas['3590'],'3590/'+header+'.aln'):
    abort('ERROR from inprosp.py: no alignments for 3590')

msas['e10'] =gen_msa(seq,header,Eco=0.001,maxId=0.90,minId=0.10,red=0.95,blastOuts=blastOuts)
if not output_msa(msas['e10'],'e10/'+header+'.aln'):
    abort('ERROR from inprosp.py: no alignments for e10')

gen_mtx(seq,'profile/'+header+'.mtx',header=header,blastOuts=blastOuts) #profile

prediction=run_psipred(fastaf,mtxf='profile/'+header+'.mtx')
output_sec(seq,'seq/'+header+'.SEQ',prediction=prediction)

pastry('tar jcf '+header+'.in.prospector.tar.bz2 3590 e10 profile seq')
pastry('/bin/mv '+header+'.in.prospector.tar.bz2 '+outdir)

sys.exit(0)
