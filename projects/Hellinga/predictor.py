#!/usr/bin/python

'''
allows to pass precomputed seq.dat file
'''

import os,sys
from time import sleep
from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import job,pastry,genJobList
from seq.fastaManager import importFastaEntry

root='/library/jose/code/python/combo_jobs/predictor'
rootH='/gpfs1/active/jose/code/projects/Hellinga'
#parse input options
inp=inpHand('Usage: predictor_ab.py [options]\ncalculate ab initio and regular structure prediction\n',
            ' -a _RA_fastaf0 (globbed) sequence file in fasta format',
            ' -b _RA_outdir0 (globbed) output directory',
            ' -c _A_listf file with list of headers',
            ' -d __header single header (in place of list of headers)',        
            ' -e __co sequence cutoff for mkpair(def:1.0)',
            ' -f __wt walltime for tasser+spicker together, in days (def:7)',
            ' -g _AR_seqdatf0 pass (globbed) seq.dat file so that we don\'t calculate it',
            ' -i __extraflags extra flags (def:none)',
            ' -j _A_inddf in.dd file (def:none)',
            ' -k _A_chaindatf0 globbed chain.dat file (def:none, dat_ab.py will generate)',
            ' -l __co2 remove homologs above sequence cutoff for prospector  (def:1.0)'
            )
inp.parse(locals(),sys.argv)

#Parse input
if not co:co='1.0'
if not wt:wt='7'
if not inddf: inddf=codedir+'/projects/Hellinga/brun/in.dd'
if not co2:co2='1.0'
else: co2=int(co2)
#if not extraflags: extraflags='-Wx=qos:critical'
joblist=genJobList()

if not listf and not header:inp.abort('At least provide a header or a list of headers')
if not listf:
    list=[header,]
else: list=chomp(open(listf,'r').readlines())

for header in list:
    
    chaindatf=''
    if chaindatf0: chaindatf=deglobb(chaindatf0,header)

    fastaf=deglobb(fastaf0,header)
    outdir=deglobb(outdir0,header)
    seqdatf=deglobb(seqdatf0,header)
    header,seq=importFastaEntry(open(fastaf,'r'))

    jobID={} #job id of the diferent submitted jobs

    #blast jobs
    jobID['blasts']=job(exe='blasts.py',exed=root,args='-a '+fastaf+' -b '+outdir).qsub('j1'+header,outdir,wallt='0.10',mem_limit='1800',extraflags='-Wx=qos:blast',joblist=joblist)
    if jobID['blasts']=='':Bye('blast job could not be submitted')
    blastout=outdir+'/'+header+'.out.blast.tar'
    expdat=outdir+'/exp.dat'

    #generate prospector input files (here we only need to create header.SEQ)
    if jobID['blasts']=='COMPLETED': dependencies=[] #previous job was completed earlier
    else: dependencies=[jobID['blasts'],]
    args='-a '+fastaf+' -b '+blastout+' -c '+outdir+' -d '+seqdatf
    jobID['inprosp']=job(exe='inprosp.py',exed=rootH,args=args).qsub('j2'+header,outdir,wallt='0.25',mem_limit='40',depend=dependencies,joblist=joblist,extraflags=extraflags)
    if jobID['inprosp']=='':Bye('inprosp job could not be submitted')
    inprosp=outdir+'/'+header+'.in.prospector.tar'

    #generate pair files and prospector run
    if jobID['inprosp']=='COMPLETED': dependencies=[]
    else: dependencies=[jobID['inprosp'],]    
    jobID['pair']=job(exe='pair.py',exed=root,args='-a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+co+' -e '+outdir).qsub('j3'+header,outdir,wallt='0.50',mem_limit='400',depend=dependencies,joblist=joblist,extraflags=extraflags)
    if jobID['pair']=='':Bye('pair job could not be submitted')
    pairtar=outdir+'/'+header+'.pair.tar'
    jobID['prosp']=job(exe='prosp.py',exed=root,args='-a '+fastaf+' -b '+inprosp+' -c '+outdir+' -d '+co2,shared='/tmp/jose/templatedir:codedir/f77/prospector/prospector_genomesmay06/templatedir.tbz2').qsub('j4'+header,outdir,wallt='0.99',mem_limit='1100',depend=dependencies)
    if jobID['prosp']=='':Bye('prosp job could not be submitted')
    prosptar=outdir+'/'+header+'.out.prospector.tar'

    #run mkrmsinp and fake mkdat (this job has no dependencies)
    args='-a '+fastaf+' -b '+outdir+' -c '+seqdatf
    if chaindatf: args+=' -d '+chaindatf
    jobID['dat_ab']=job(exe='dat_ab.py',exed=rootH,args=args).qsub('j5'+header,outdir,wallt='0.10',mem_limit='40',depend=dependencies,joblist=joblist,extraflags=extraflags)
    if jobID['dat_ab']=='':Bye('dat_ab job could not be submitted')
    dat_abtar=outdir+'/'+header+'.dat_ab.tar'

    #run mkrmsinp and  mkdat
    if jobID['prosp']=='COMPLETED': dependencies=[]
    else: dependencies=[jobID['prosp']]
    args='-a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+prosptar+' -e '+outdir+' -c '+seqdatf+' -f '+co2
    jobID['dat']=job(exe='dat.py',exed=rootH,args=args).qsub('j6'+header,outdir,wallt='0.99',mem_limit='800',depend=dependencies,joblist=joblist,extraflags=extraflags)
    if jobID['dat']=='':Bye('dat job could not be submitted')
    dattar=outdir+'/'+header+'.dat.tar'

    #run ab initio tasser and spicker
    dependencies=[]
    if jobID['pair']!='COMPLETED': dependencies.append(jobID['pair'])
    if jobID['dat_ab'] !='COMPLETED': dependencies.append(jobID['dat_ab'])
    args='-a '+fastaf+' -b '+dat_abtar+' -c '+pairtar+' -d '+inprosp+' -e '+expdat+' -f '+outdir+' -i '+inddf
    jobID['ab_tasser']=job(exe='tas_spk_ab.py',exed=rootH,args=args).qsub('j7'+header,outdir,wallt=wt,mem_limit='900',depend=dependencies,joblist=joblist,extraflags=extraflags)
    if jobID['ab_tasser']=='':Bye('tasser job could not be submitted')
    ab_spktar=outdir+'/'+header+'.out.ab_spicker.tar'

    #run tasser and spicker
    dependencies=[]
    if jobID['pair']!='COMPLETED': dependencies.append(jobID['pair'])
    if jobID['dat'] !='COMPLETED': dependencies.append(jobID['dat'])
    args='-a '+fastaf+' -b '+dattar+' -c '+pairtar+' -d '+inprosp+' -e '+expdat+' -f '+outdir+' -g '+inddf
    jobID['tasser']=job(exe='tas_spk.py',exed=rootH,args=args).qsub('j8'+header,outdir,wallt=wt,mem_limit='900',depend=dependencies,joblist=joblist,extraflags=extraflags)
    if jobID['tasser']=='':Bye('tasser job could not be submitted')
    spktar=outdir+'/'+header+'.out.spicker.tar'


    #run dukka and pulchra (assume spicker outputs 10 clusters, submit one job per combo)
    if jobID['ab_tasser']=='COMPLETED': dependencies=None
    else: dependencies=[jobID['ab_tasser'],]
    args='-a '+fastaf+' -b '+ab_spktar+' -c '+outdir+' -f '+header+'abdukka.tar'
    job(exe='dukka_pulchra.py',exed=root,args=args).qsub('j9'+header,outdir,wallt=0.25,mem_limit='40',depend=dependencies,joblist=joblist,extraflags=extraflags)

    #run dukka and pulchra (assume spicker outputs 10 clusters, submit one job per combo)
    if jobID['tasser']=='COMPLETED': dependencies=None
    else: dependencies=[jobID['tasser'],]
    args='-a '+fastaf+' -b '+spktar+' -c '+outdir+' -f '+header+'dukka.tar'
    job(exe='dukka_pulchra.py',exed=root,args=args).qsub('j10'+header,outdir,wallt=0.25,mem_limit='40',depend=dependencies,joblist=joblist,extraflags=extraflags)




sys.exit(0)
