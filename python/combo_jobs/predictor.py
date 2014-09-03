#!/usr/bin/python

import os,sys
from time import sleep

from utilities.codedir import codedir
from utilities.small_utilities import Bye
from inputArgs.inputArgs import inpHand
from jobs.job import pastry

from seq.fastaManager import importFastaEntry
from seq.blastManager import blastRun
from seq.msa import gen_msa,output_msa
from seq.mtx import gen_mtx

from secondary.psipred.runpsipred import run_psipred
from secondary.mkseq import output_sec

from tasser.mkpair import gen_pair
from tasser.mksolve import gen_solve

#parse input options
inpHand('Usage: predictor.py [options]',
        '  -a _A_fastaf sequence file in fasta format',
        '  -b __co sequence cutoff (def:0.3)'
        ).parse(locals(),sys.argv)

if not co:co=0.3

print '\ngather header and sequence from fasta file\n';sleep(1)
header,seq=importFastaEntry(open(fastaf,'r')) #;print header,seq;sys.exit(1)

print '\ndo blast run with default options\n';sleep(1)
print 'blastRun(fastaf,header=header)';sleep(1)
blastOuts=blastRun(fastaf,header=header) #;print blastOuts;sys.exit(1)
cmd='tar cf '+header+'.out.blast.tar *.chk *.fasta *.blast' ;print cmd ;sleep(1) ;pastry(cmd)

print '\ncreate input files for prospector\n';sleep(1)
msas={}
os.system('/bin/mkdir -p 3590 e10 profile seq')
print 'gen_msa(seq,header,Eco=0.001,maxId=0.90,minId=0.35,red=0.95,blastOuts=blastOuts)';sleep(1)
msas['3590']=gen_msa(seq,header,Eco=0.001,maxId=0.90,minId=0.35,red=0.95,blastOuts=blastOuts) #mult.seq.alg.
output_msa(msas['3590'],'3590/'+header+'.aln') #close multiple alignemt #output the mult.seq.alg. to a file
print 'gen_msa(seq,header,Eco=0.001,maxId=0.90,minId=0.10,red=0.95,blastOuts=blastOuts)';sleep(1)
msas['e10'] =gen_msa(seq,header,Eco=0.001,maxId=0.90,minId=0.10,red=0.95,blastOuts=blastOuts)
output_msa(msas['e10'],'e10/'+header+'.aln')   #distant multiple alignment
print 'gen_mtx(seq,\'profile/\'+header+\'.mtx\',header=header,blastOuts=blastOuts)';sleep(1)
gen_mtx(seq,'profile/'+header+'.mtx',header=header,blastOuts=blastOuts) #profile
print 'run_psipred(fastaf,mtxf=\'profile/\'+header+\'.mtx\')';sleep(1)
prediction=run_psipred(fastaf,mtxf='profile/'+header+'.mtx')#run psipred with generated profile
output_sec(seq,'seq/'+header+'.SEQ',prediction=prediction) #format psipred output
cmd='tar cf '+header+'.in.prospector.tar 3590 e10 profile seq' ;print cmd ;sleep(1) ;pastry(cmd)

print '\nrun prospector job (we NEED /tmp/jose/templatedir !!!)\n';sleep(1)
prospector=codedir+'/python/prospector/prospector_genomesmay06/modified/prospector.py'
cmd=prospector+' -a '+header+' -b . ' ;print cmd ;sleep(1) ;os.system(cmd)
cmd='tar cf '+header+'.out.prospector.tar pdbb*' ;print cmd ;sleep(1) ;pastry(cmd)

print '\ncreate input for tasser\n'; sleep(1)
cmd='mkrmsinp.py -i '+header+' -o . -f '+fastaf; print cmd ;sleep(1) ;os.system(cmd) #rmsinp file
cmd='mkdat.py -a '+header+' -b '+fastaf+' -c seq/'+header+'.SEQ -e pdbbpdborienrev/'+header+'rap3orienrev5s.pdb -g pdbbzpot4a/'+header+'.predictedrap3orienrev -l '+blastOuts['blast']+' -m '+`co`
print cmd ;sleep(1) ;pastry(cmd) #comb*.dat dist*.dat par.dat chain.dat
print 'gen_solve(fastaf)';sleep(1)
gen_solve(fastaf) #create exp.dat
print 'gen_msa(seq,header,Eco=1.0  ,maxId=0.90,minId=0.10,red=0.95,blastOuts=blastOuts)';sleep(1)
msas['pair']=gen_msa(seq,header,Eco=1.0  ,maxId=0.90,minId=0.10,red=0.95,blastOuts=blastOuts)
print 'gen_pair(msas[\'pair\'],\'seq.dat\',\'.\',co=co)';sleep(1)
gen_pair(msas['pair'],'seq.dat','.',co=co) #pair1.dat and pair3.dat
root=codedir+'/f77/tasser/tasser.1.0'
os.system('/bin/cp '+root+'/original/in.dd '+root+'/common.tbz2 . && tar jxf common.tbz2')
cmd='tar cf '+header+'.in.tasser.tar comb*dat chain.dat dist*dat exp.dat pa*dat seq.dat in.dd rmsinp'; print cmd ;sleep(1) ;pastry(cmd)

print '\nrun tasser & spicker\n';sleep(1)
cas=codedir+'/f77/tasser/tasser.1.0/original/cas'
spicker=+'/f77/spicker/spicker0.0/spicker'
pastry('/bin/cp '+cas+' '+spicker+' .')
pastry('./cas')
pastry('echo 16 > tra.in && ls -1 rep*.tra >> tra.in')
pastry('./spicker')
cmd='tar cf '+header+'.in.spicker.tar rep*tra tra.in' ;print cmd ;sleep(1) ;pastry(cmd)
cmd='tar cf '+header+'.out.spicker.tar combo* closc* str.txt rst.dat' ;print cmd ;sleep(1)
pastry(cmd)

print '\nrefine built-in side-chains for combo*.pdb\n';sleep(1)
cmd='for c in `ls combo*.pdb`; do dukka.x $c && pulchra best.pdb && /bin/mv best.pdb.rebuilt ${c}.rebuilt; done'; print cmd ;sleep(1) ;pastry(cmd)
cmd='tar cf '+header+'.out.dukka_pulchra.tar combo*rebuilt' ;print cmd ;sleep(1) ;pastry(cmd)

sys.exit(0)
