#!/usr/bin/python

import os,sys

from time import sleep
from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import job,pastry,add_dir_for_backup,genJobList
from seq.fastaManager import importFastaEntry

root='/library/jose/code/python/combo_jobs/predictor'

#parse input options
inp=inpHand('Usage: predictor.py [options]',
            ' -a _RA_fastaf0 (globbed) sequence file in fasta format',
            ' -b _RA_outdir0 (globbed) output directory',
            ' -c _A_listf file with list of headers',
            ' -d __header single header (in place of -a option)',        
            ' -e __co sequence cutoff for mkdat and mkpair(def:0.35)',            
            ' -f __wt walltime for tasser+spicker together, in days (def:6.99)',
            ' -g __co2 remove homologs above sequence cutoff for prospector  (def:0.35)'
            )
inp.parse(locals(),sys.argv)

#Parse input
if not co:co='0.35'
if not co2:co2='0.35'
if not wt:wt='6.99'

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

    #blast jobs (including exp.dat)
    #jobID['blasts']=job(exe='blasts.py',exed=root,args='-a '+fastaf+' -b '+outdir).qsub('j1'+header,outdir,wallt='0.162',mem_limit='900',extraflags='',joblist=joblist)
    #if jobID['blasts']=='':Bye('blast job could not be submitted')
    blastout=outdir+'/'+header+'.out.blast.tar'
    expdat=outdir+'/exp.dat'
#    add_dir_for_backup(outdir,depth=2)    

    #generate prospector input files
    #if jobID['blasts']=='COMPLETED': dependencies=[] #previous job was completed earlier
    #else: dependencies=[jobID['blasts'],]
    #jobID['inprosp']=job(exe='inprosp.py',exed=root,args='-a '+fastaf+' -b '+blastout+' -c '+outdir).qsub('j2'+header,outdir,wallt='0.249',mem_limit='40',depend=dependencies,joblist=joblist)
    #if jobID['inprosp']=='':Bye('inprosp job could not be submitted')
    inprosp=outdir+'/'+header+'.in.prospector.tar'
#    add_dir_for_backup(outdir,depth=2)

    #generate pair files
    if jobID['inprosp']=='COMPLETED': dependencies=[]
    else: dependencies=[jobID['inprosp'],]    
    jobID['pair']=job(exe='pair.py',exed=root,args='-a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+co+' -e '+outdir).qsub('j3'+header,outdir,wallt='0.999',mem_limit='400',depend=dependencies,joblist=joblist)
    if jobID['pair']=='':Bye('pair job could not be submitted')
    pairtar=outdir+'/'+header+'.pair.tar'
#    add_dir_for_backup(outdir,depth=2)

    #run prospector
    jobID['prosp']=job(exe='prosp.py',exed=root,args='-a '+fastaf+' -b '+inprosp+' -c '+outdir+' -d '+co2,shared='/tmp/jose/templatedir:codedir/f77/prospector/prospector_genomesmay06/templatedir.tbz2').qsub('j4'+header,outdir,wallt='0.249',mem_limit='900',depend=dependencies,joblist=joblist)
    if jobID['prosp']=='':Bye('prosp job could not be submitted')
    prosptar=outdir+'/'+header+'.out.prospector.tar'
#    add_dir_for_backup(outdir,depth=2)

    #run mkrmsinp and mkdat
    if jobID['prosp']=='COMPLETED': dependencies=[]
    else: dependencies=[jobID['prosp']]
    jobID['dat']=job(exe='dat.py',exed=root,args='-a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+prosptar+' -e '+outdir +' -f '+co).qsub('j5'+header,outdir,wallt='0.99',mem_limit='100',depend=dependencies,joblist=joblist)
    if jobID['dat']=='':Bye('dat job could not be submitted')
    dattar=outdir+'/'+header+'.dat.tar'
#    add_dir_for_backup(outdir,depth=2)

    #run tasser and spicker
    dependencies=[]
    if jobID['pair']!='COMPLETED': dependencies.append(jobID['pair'])
    if jobID['dat'] !='COMPLETED': dependencies.append(jobID['dat'])
    jobID['tasser']=job(exe='tas_spk.py',exed=root,args='-a '+fastaf+' -b '+dattar+' -c '+pairtar+' -d '+inprosp+' -e '+expdat+' -f '+outdir).qsub('j6'+header,outdir,wallt=wt,mem_limit='900',depend=dependencies,joblist=joblist)
    if jobID['tasser']=='':Bye('tasser job could not be submitted')
    spktar=outdir+'/'+header+'.out.spicker.tar'
#    add_dir_for_backup(outdir,depth=2)

    #run dukka and pulchra (assume spicker outputs 10 clusters, submit one job per combo)
    if jobID['tasser']=='COMPLETED': dependencies=None
    else: dependencies=[jobID['tasser'],]
    job(exe='dukka_pulchra.py',exed=root,args='-a '+fastaf+' -b '+spktar+' -c '+outdir).qsub('j7'+header,outdir,wallt=0.249,mem_limit='100',depend=dependencies,joblist=joblist)
#    add_dir_for_backup(outdir,depth=2)

sys.exit(0)
