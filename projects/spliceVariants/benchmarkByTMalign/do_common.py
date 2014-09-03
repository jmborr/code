#!/usr/bin/python

'''

For every target, get exp.dat seq.dat pair1.dat pair3.dat
rmsinp. These files are intrinsic to the target and independent of the
template we provide.

For some of the targets, Shashi already did the input, so we take it
from him. For the other targets, Adrian did seq.dat and we do the
rest.

'''
import os,sys
from time import sleep
from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp,junkName
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import job,pastry,add_dir_for_backup
from seq.fastaManager import importFastaEntry

currd=os.getcwd()
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
shashid='/gpfs1/scratch/shashi'
root='/library/jose/code/python/combo_jobs/predictor'
#benchd='/gpfs1/archive/adrian/dat/pdb_aug06.apr06'
benchd='/gpfs1/archive/adrian/dat/pdb_mar07'

#Some variables
co='0.3' #sequence cutoff for mkpair
fastaf0=benchd+'/input/xxxxx'
currd=os.getcwd()

#input list of shashi's already done inputs
slist=chomp(open(scd+'/shashi.list','r').readlines())
sall=''.join(slist)
junkd=junkName() ; pastry('/bin/mkdir '+junkd)

#for line in open(scd+'/findCommonTemplate.dat','r').readlines():
for header in chomp(open(scd+'/filtered.list','r').readlines()):

    #resolved some header-dependent variables
    fastaf=deglobb(fastaf0,header)
    if not os.path.exists(fastaf):
        fastaf=os.path.join(benchd,'aux',header)
    outdir=scd+'/out/'+header
    if not os.path.exists(header): os.system('/bin/mkdir -p '+outdir)

    jobID={} #job id of the diferent submitted jobs

    #check if we already did this header
    already_done=True
    for f in ('rmsinp','exp.dat','pair1.dat','pair3.dat','seq.dat'):
        if not os.path.exists( os.path.join(outdir,f) ):
            already_done=False
            break
    if already_done:
        print header,'already done'
        continue #go with next header
    
    #check if shashi has already done the pair files
    nextHeader=0
    for line in slist: #find which bench-x00 subdirectory is the input
        if header in line:
            print line
            sinp=shashid+line[1:] #tasser input file
            os.system('tar jxf '+sinp+' -C '+junkd)
            os.chdir(junkd)
            os.system('/bin/mv rmsinp exp.dat pair1.dat pair3.dat seq.dat '+outdir)
            os.chdir(currd)
            nextHeader=1
            break
    if nextHeader:
        print header,'bringing stuff from shashi'
        continue #go with next header
    
    #nothing in Shashi's inputs, thus do them
    #blast jobs (including exp.dat)
    jobID['blasts']=job(exe='blasts.py',exed=root,args='-a '+fastaf+' -b '+outdir).qsub('j1'+header,outdir,wallt='0.249',mem_limit='900',extraflags='-Wx=qos:personal')
    if jobID['blasts']=='':Bye('blast job could not be submitted')
    blastout=outdir+'/'+header+'.out.blast.tar.bz2'
    expdat=outdir+'/exp.dat'

    #tar mtx and seq files (for pair.py) and output seq.dat
    missing_pdb=False
    for f in (benchd+'/seq/'+header+'.SEQ', benchd+'/profile/'+header+'.mtx'):
        if not os.path.exists(f):
            missing_pdb=True
            break
    if missing_pdb:
        print header, 'no .SEQ and .mtx files'
        continue #go with next header
    os.chdir(outdir)
    os.system('/bin/cp '+benchd+'/seq/'+header+'.SEQ ./seq.dat')
    os.system('/bin/mkdir profile seq') #temporary directories
    os.system('/bin/cp '+benchd+'/profile/'+header+'.mtx profile')
    os.system('/bin/cp '+benchd+'/seq/'+header+'.SEQ seq')
    os.system('/bin/tar jcf '+header+'.in.prospector.tar.bz2 profile seq')
    os.system('/bin/rm -rf profile seq')
    inprosp=outdir+'/'+header+'.in.prospector.tar.bz2'
    os.chdir(currd)

    #generate pair files
    dependencies=[]
    if jobID['blasts']!='COMPLETED':dependencies.append(jobID['blasts'])
    jobID['pair']=job(exe='pair.py',exed=root,args='-a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+co+' -e '+outdir+' -f no').qsub('j3'+header,outdir,wallt='1.99',mem_limit='400',depend=dependencies)
    if jobID['pair']=='':Bye('pair job could not be submitted')    

    #generate mkrsinp
    os.system('mkrmsinp.py -i '+header+' -o '+outdir+' -f '+fastaf)
    
os.system('/bin/rm -rf '+junkd)
sys.exit(0)
