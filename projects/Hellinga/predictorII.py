#!/usr/bin/python

'''
allows to pass precomputed blast and exp.dat results
Example:
 /gpfs1/active/jose/code/projects/Hellinga/predictorII.py
 -a $g/xxxxx/xxxxx -b $g/xxxxx -c $g/best.list
 -g $g/xxxxx/seq.dat -i "-Wx=qos:critical"
 -m $g/xxxxx/xxxxx.out.blast.tar -n $g/xxxxx/exp.dat
'''

import os,sys
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
            ' -j _A_inddf in.dd file (def:codedir/projects/Hellinga/brun/in.dd)',
            ' -k _A_chaindatf0 globbed chain.dat file (def:none, dat_ab.py will generate)',
            ' -l __co2 remove homologs above sequence cutoff for prospector  (def:1.0)',
            ' -m _RA_blastout0 (globbed) blastout tarred results from blast.py script',
            ' -n _RA_expdat0 (globbed) exp.dat file also result from blast.py scripte'
            )
inp.parse(locals(),sys.argv)

#Parse input
if not co:co='1.0'
if not wt:wt='7'
if not inddf: inddf=codedir+'/projects/Hellinga/brun/in.dd'
if not co2:co2='1.0'
else: co2=int(co2)
if not extraflags: extraflags='-Wx=qos:critical'

joblist=genJobList()

if not listf and not header:inp.abort('At least provide a header or a list of headers')
if not listf:
    list=[header,]
else: list=chomp(open(listf,'r').readlines())

for header in list:
    
    chaindatf=''
    if chaindatf0: chaindatf=deglobb(chaindatf0,header)

    fastaf=deglobb(fastaf0,header)
    seqdatf=deglobb(seqdatf0,header)
    blastout=deglobb(blastout0,header)
    expdatf=deglobb(expdat0,header)
    
    #check existence of input files
    notinput=False #flag there's not enough input to submit job
    for file in (fastaf,seqdatf,blastout,expdatf):
        if not os.path.exists(file):
            sys.stderr.write(file+' does not exists\n')
            notinput=True
            break
    if notinput: continue #go to another header
    
    outdir=deglobb(outdir0,header)
    jobID={} #job id of the diferent submitted jobs
    
    #generate prospector input files and pair files
    dependencies=[] #previous job was completed earlier
    args='-a '+fastaf+' -b '+outdir+' -c '+header+' -d '+co+' -e '+seqdatf+' -f '+blastout
    jobID['inprosp_pair']=job(exe='inprosp_pair.py',exed=rootH,args=args).qsub('k1'+header,outdir,wallt='0.50',mem_limit='400',depend=dependencies,joblist=joblist,extraflags='-Wx=qos:blast')
    if jobID['inprosp_pair']=='':Bye('inprosp_pair job could not be submitted')
    inprosp=outdir+'/'+header+'.in.prospector.tar'
    pairtar=outdir+'/'+header+'.pair.tar'

    #stablish dependencies for tasser and for ab initio tasser
    if jobID['inprosp_pair']=='COMPLETED': dependencies=[]
    else: dependencies=[jobID['inprosp_pair']]


    #do ab initio tasser run (mkdat+cas+spicker+dukka)
    args='-a '+fastaf+' -b '+inprosp+' -c '+outdir+' -d '+header+' -g '+expdatf+' -i '+pairtar+' -j '+inddf+' -k '+seqdatf
    if chaindatf: args+=' -l '+chaindatf
    jobID['ab_pro_dat_tas_skp_duk']=job(exe='ab_pro_dat_tas_skp_duk.py',exed=rootH,args=args).qsub('k3'+header,outdir,wallt='7',mem_limit='1020',depend=dependencies,joblist=joblist,extraflags=extraflags)
    if jobID['ab_pro_dat_tas_skp_duk']=='':Bye('ab initio tasser could not be submitted')

    dat_abtar=outdir+'/'+header+'.dat_ab.tar'
    ab_spktar=outdir+'/'+header+'.out.ab_spicker.tar'
    abduktar=outdir+'/'+header+'abdukka.tar'

    #do a normal tasser run (prosp+mkdat+cas+spicker+dukka)
    #args='-a '+fastaf+' -b '+inprosp+' -c '+outdir+' -d '+header+' -e '+co+' -f '+blastout+' -g '+expdatf+' -i '+pairtar+' -j '+inddf
    #jobID['pro_dat_tas_skp_duk']=job(exe='pro_dat_tas_skp_duk.py',exed=rootH,args=args,shared='/tmp/jose/templatedir:codedir/f77/prospector/prospector_genomesmay06/templatedir.tbz2').qsub('k2'+header,outdir,wallt='7',mem_limit='1020',depend=dependencies,joblist=joblist)
    #if jobID['pro_dat_tas_skp_duk']=='':Bye('normal tasser could not be submitted')
    
    #prosptar=outdir+'/'+header+'.out.prospector.tar'
    #dattar=outdir+'/'+header+'.dat.tar'
    #spktar=outdir+'/'+header+'.out.spicker.tar'
    #duktar=outdir+'/'+header+'dukka.tar'

sys.exit(0)
