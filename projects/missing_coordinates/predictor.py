#!/usr/bin/python

'''
'''

import os,sys
from time import sleep
from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import job,pastry,genJobList
from seq.fastaManager import importFastaEntry
from tasser.utilities import calc_wall_time

root='/library/jose/code/python/combo_jobs/predictor'
rootH='/gpfs1/active/jose/code/projects/missing_coordinates'
#parse input options
inp=inpHand('Usage: predictor.py [options]\n\n',
            ' -a _RA_fastaf0 (globbed) sequence file in fasta format',
            ' -b _RA_outdir0 (globbed) output directory',
            ' -c _A_listf file with list of headers',
            ' -d __header single header (in place of list of headers)',        
            ' -i __extraflags extra flags (def:none)',
            ' -l __co2 remove homologs above sequence cutoff for prospector  (def:1.0)'
            )
inp.parse(locals(),sys.argv)

#Parse input
if not co2:co2='1.0'
if float(co2)>1: inp.abort('ERROR: cutoff for prospector should be less than one\n')
#if not extraflags: extraflags='-Wx=qos:critical'
joblist=genJobList()

if not listf and not header:inp.abort('At least provide a header or a list of headers')
if not listf:
    list=[header,]
else: list=chomp(open(listf,'r').readlines())

for header in list:
    
    fastaf=deglobb(fastaf0,header)
    outdir=deglobb(outdir0,header)
    header,seq=importFastaEntry(open(fastaf,'r'))

    jobID={} #job id of the diferent submitted jobs

    #blast jobs
    blastout=outdir+'/'+header+'.out.blast.tar.bz2'
    expdat=outdir+'/exp.dat'
    if not os.path.exists(blastout) or not os.path.exists(expdat):
        jobID['blasts']=job(exe='blasts.py',exed=root,args='-a '+fastaf+' -b '+outdir).qsub('j1'+header,outdir,wallt='0.10',mem_limit='900',extraflags='-Wx=qos:blast',joblist=joblist)
        if jobID['blasts']=='':Bye('blast job could not be submitted')
    else:
        jobID['blasts']='COMPLETED'

    #generate prospector input files (here we only need to create header.SEQ)
    if jobID['blasts']=='COMPLETED': dependencies=[] #previous job was completed earlier
    else: dependencies=[jobID['blasts'],]
    args='-a '+fastaf+' -b '+blastout+' -c '+outdir
    jobID['inprosp']=job(exe='inprosp.py',exed=root,args=args).qsub('j2'+header,outdir,wallt='0.99',mem_limit='40',depend=dependencies,joblist=joblist,extraflags=extraflags)
    if jobID['inprosp']=='':Bye('inprosp job could not be submitted')
    inprosp=outdir+'/'+header+'.in.prospector.tar.bz2'
    
    #generate prospector run
    if jobID['inprosp']=='COMPLETED': dependencies=[]
    else: dependencies=[jobID['inprosp'],]

    prosptar=outdir+'/'+header+'.out.prospector.tar.bz2'
    jobID['prosp']=job(exe='prosp.py',exed=root,args='-a '+fastaf+' -b '+inprosp+' -c '+outdir+' -d '+co2,shared='/tmp/jose/templatedir:codedir/f77/prospector/prospector_genomesmay06/templatedir.tbz2').qsub('j4'+header,outdir,wallt='0.99',mem_limit='900',depend=dependencies,joblist=joblist)
    if jobID['prosp']=='':Bye('prosp job could not be submitted')

sys.exit(0)
