#!/usr/bin/python

'''
allows to pass precomputed blast and exp.dat results
Example:
for g in a0216 a0512 a0653 a0857 a0984;do /gpfs1/active/jose/code/projects/Hellinga/predictorIII.py -a $g/xxxxx/xxxxx -b $g/xxxxx -g $g/xxxxx/seq.dat -m $g/xxxxx/xxxxx.out.blast.tar  -n $g/xxxxx/exp.dat -c $g/toy.list -i "-Wx=qos:critical" ;done
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
            ' -f __wt walltime for tasser+spicker together, in days (def:3)',
            ' -g _AR_seqdatf0 pass (globbed) seq.dat file so that we don\'t calculate it',
            ' -i __extraflags extra flags (def:none)',
            ' -j _A_inddf in.dd file (def:codedir/projects/Hellinga/brun/in.dd)',
            ' -k _A_chaindatf0 globbed chain.dat file (def:none, dat_ab.py will generate)',
            ' -l __co2 remove homologs above sequence cutoff for prospector  (def:1.0)',
            ' -m _RA_blastout0 (globbed) blastout tarred results from blast.py script',
            ' -n _RA_expdat0 (globbed) exp.dat file also result from blast.py scripte',
            ' -o __zzz lag between qsub calls, in seconds (def:1)',
            ' -p __prefix prefix to put to the jobnames (def:'')',
            )
inp.parse(locals(),sys.argv)

#Parse input
if not co:co='1.0'
if not wt:wt='0.99'
if not inddf: inddf=codedir+'/projects/Hellinga/brun/in.dd'
if not co2:co2='1.0'
else: co2=int(co2)
if not extraflags: extraflags='' #'-Wx=qos:critical'
if not zzz: zzz=0.2
else: zzz=float(zzz)
if not prefix: prefix=''

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
    dependencies=[] #previous job was completed earlier

    #do ab initio tasser run (mkdat+cas+spicker+dukka)
    args='-a '+fastaf+' -b '+blastout+' -c '+outdir+' -d '+header+' -g '+expdatf+' -i '+co+' -j '+inddf+' -k '+seqdatf
    if chaindatf: args+=' -l '+chaindatf
    jobID['ab_pair_pro_dat_tas_skp_duk']=job(exe='ab_pair_pro_dat_tas_skp_duk.py',exed=rootH,args=args).qsub(prefix+'k3'+header,outdir,wallt=wt,mem_limit='850',depend=dependencies,joblist=joblist,extraflags=extraflags,zzz=zzz)
    if jobID['ab_pair_pro_dat_tas_skp_duk']=='':
        sys.stderr.write('ab initio tasser for '+header+' could not be submitted\n')

    inprosp=outdir+'/'+header+'.in.prospector.tar'
    pairtar=outdir+'/'+header+'.pair.tar'
    dat_abtar=outdir+'/'+header+'.dat_ab.tar'
    in_ab_tassertar=outdir+'/'+header+'.in.ab_tasser.tar'
    inab_spickertar=outdir+'/'+header+'.in.ab_spicker.tar '
    ab_spktar=outdir+'/'+header+'.out.ab_spicker.tar'
    abduktar=outdir+'/'+header+'abdukka.tar'

    
sys.exit(0)
