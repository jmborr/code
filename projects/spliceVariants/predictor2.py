#!/usr/bin/python

import os,sys

from time import sleep
from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import job,pastry,add_dir_for_backup,genJobList
from seq.fastaManager import importFastaEntry

root1='/gpfs1/active/jose/code/projects/spliceVariants'
root2='/library/jose/code/python/combo_jobs/predictor'

#parse input options
inp=inpHand('Usage: predictor.py [options]',
            ' -a _RA_fastaf0 (globbed) sequence file in fasta format',
            ' -b _RA_outdir0 (globbed) output directory',
            ' -c _A_listf file with list of headers',
            ' -d __header single header (in place of -a option)',
            ' -f __wt walltime for tasser+spicker together, in days (def:7)',
            ' -i __prefix an optional prefix for the jobname (def:None)',
            )
inp.parse(locals(),sys.argv)

#Parse input
if not wt:wt='7'
if not prefix: prefix=''

if not listf and not header:inp.abort('At least provide either -a or -b flag')
if not listf:
    list=[header,]
else: list=chomp(open(listf,'r').readlines())

joblist=genJobList()

for header in list:
    fastaf=deglobb(fastaf0,header)
    outdir=deglobb(outdir0,header)

    header,seq=importFastaEntry(open(fastaf,'r'))

    jobID={} #job id of the diferent submitted jobs


    #run tasser and spicker
    dependencies=[]
    intar=outdir+'/'+header+'.in.tasser.tar'
    jobID['tasser']=job(exe='tas_spk2.py',exed=root1,args='-a '+fastaf+' -b '+intar+' -f '+outdir).qsub(prefix+'j6'+header,outdir,wallt=wt,mem_limit='900',depend=dependencies,joblist=joblist)
    if jobID['tasser']=='':Bye('tasser job could not be submitted')
    spktar=outdir+'/'+header+'.out.spicker.tar'
#    add_dir_for_backup(outdir,depth=2)

    #run dukka and pulchra (assume spicker outputs 10 clusters, submit one job per combo)
    if jobID['tasser']=='COMPLETED': dependencies=None
    else: dependencies=[jobID['tasser'],]
    job(exe='dukka_pulchra.py',exed=root2,args='-a '+fastaf+' -b '+spktar+' -c '+outdir).qsub(prefix+'j7'+header,outdir,wallt=0.249,mem_limit='40',depend=dependencies,joblist=joblist)
#    add_dir_for_backup(outdir,depth=2)

sys.exit(0)
