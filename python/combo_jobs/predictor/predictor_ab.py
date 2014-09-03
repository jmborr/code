#!/usr/bin/python

import os,sys

from time import sleep
from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import job,pastry,genJobList
from seq.fastaManager import importFastaEntry

root='/library/jose/code/python/combo_jobs/predictor'

#parse input options
inp=inpHand('Usage: predictor_ab.py [options]',
            ' -a _RA_fastaf0 (globbed) sequence file in fasta format',
            ' -b _RA_outdir0 (globbed) output directory',
            ' -c _A_listf file with list of headers',
            ' -d __header single header (in place of -a option)',        
            ' -e __co sequence cutoff for mkpair(def:0.3)',
            ' -f __wt walltime for tasser+spicker together, in days (def:7)'
            )
inp.parse(locals(),sys.argv)

#Parse input
if not co:co='0.3'
if not wt:wt='7'
joblist=genJobList()

if not listf and not header:inp.abort('At least provide either -a or -b flag')
if not listf:
    list=[header,]
else: list=chomp(open(listf,'r').readlines())

for header in list:
    fastaf=deglobb(fastaf0,header)
    outdir=deglobb(outdir0,header)

    header,seq=importFastaEntry(open(fastaf,'r'))

    jobID={} #job id of the diferent submitted jobs

    #blast jobs
    jobID['blasts']=job(exe='blasts.py',exed=root,args='-a '+fastaf+' -b '+outdir).qsub('j1'+header,outdir,wallt='0.10',mem_limit='1333',extraflags='-Wx=qos:blast',joblist=joblist)
    if jobID['blasts']=='':Bye('blast job could not be submitted')
    blastout=outdir+'/'+header+'.out.blast.tar'
    expdat=outdir+'/exp.dat'

    #generate prospector input files (here we only need to create header.SEQ)
    if jobID['blasts']=='COMPLETED': dependencies=[] #previous job was completed earlier
    else: dependencies=[jobID['blasts'],]
    jobID['inprosp']=job(exe='inprosp_ab.py',exed=root,args='-a '+fastaf+' -b '+blastout+' -c '+outdir).qsub('j2'+header,outdir,wallt='0.25',mem_limit='40',depend=dependencies,joblist=joblist)
    if jobID['inprosp']=='':Bye('inprosp job could not be submitted')
    inprosp=outdir+'/'+header+'.in.prospector.tar'

    #generate pair files
    if jobID['inprosp']=='COMPLETED': dependencies=[]
    else: dependencies=[jobID['inprosp'],]    
    jobID['pair']=job(exe='pair.py',exed=root,args='-a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+co+' -e '+outdir).qsub('j3'+header,outdir,wallt='0.50',mem_limit='400',depend=dependencies,joblist=joblist)
    if jobID['pair']=='':Bye('pair job could not be submitted')
    pairtar=outdir+'/'+header+'.pair.tar'

    #run mkrmsinp and fake mkdat (this job has no dependencies)
    jobID['dat']=job(exe='dat_ab.py',exed=root,args='-a '+fastaf+' -b '+outdir).qsub('j5'+header,outdir,wallt='0.10',mem_limit='40',depend=dependencies,joblist=joblist)
    if jobID['dat']=='':Bye('dat_ab job could not be submitted')
    dattar=outdir+'/'+header+'.dat.tar'

    #run tasser and spicker
    dependencies=[]
    if jobID['pair']!='COMPLETED': dependencies.append(jobID['pair'])
    if jobID['dat'] !='COMPLETED': dependencies.append(jobID['dat'])
    jobID['tasser']=job(exe='tas_spk_ab.py',exed=root,args='-a '+fastaf+' -b '+dattar+' -c '+pairtar+' -d '+inprosp+' -e '+expdat+' -f '+outdir).qsub('j6'+header,outdir,wallt=wt,mem_limit='900',depend=dependencies,joblist=joblist)
    if jobID['tasser']=='':Bye('tasser job could not be submitted')
    spktar=outdir+'/'+header+'.out.spicker.tar'

    #run dukka and pulchra (assume spicker outputs 10 clusters, submit one job per combo)
    if jobID['tasser']=='COMPLETED': dependencies=None
    else: dependencies=[jobID['tasser'],]
    job(exe='dukka_pulchra.py',exed=root,args='-a '+fastaf+' -b '+spktar+' -c '+outdir).qsub('j7'+header,outdir,wallt=0.25,mem_limit='40',depend=dependencies,joblist=joblist)

sys.exit(0)
