#!/usr/bin/python
'''
All inputs for tasser have already been done
'''
import os,sys

from time import sleep
from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import job,pastry,genJobList
from seq.fastaManager import importFastaEntry
from random import randint

def calc_wt(L): #estimate TASSER running time
    minwt=0.99    
    if L<100: wt=minwt
    else: wt=minwt+(L-100.0)/200*6
    return '%4.2lf'%(wt,)

joblist=genJobList()

currd=os.getcwd()
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
pcd='/gpfs1/active/jose/code/projects/spliceVariants/benchmarkByTMalign'
root='/library/jose/code/python/combo_jobs/predictor'
benchd='/gpfs1/archive/adrian/dat/pdb_mar07'

#Some options
wt='6.99' #walltime for tasser+spicker together
fastaf0=benchd+'/input/xxxxx'
CA0=benchd+'/CA/xxxxx.pdb'
#job dependencies
#blast->pair->tasser->dukka

#create list of header and output pairs
list=[]
for line in open(scd+'/findCommonTemplate.dat','r').readlines():
    if line[0]=='#': continue #line is just a comment
    (short,long,templ)=line.split()
    subdir=scd+'/out/'+short+'_'+long
    for header in (short,long):
        outdir=subdir+'/'+header
        list.append((header,outdir))

#go through every (header,outdir) pair
for i in range(len(list)):
#for i in range(1):
    jobname='%05d'%(i) #because multiple jobs have same header
    header,outdir=list[i]
    fastaf=deglobb(fastaf0,header)
    CA=deglobb(CA0,header)
    if not os.path.exists(CA): #pdb_mar07 does not contains this entry which was from aug06
        CA=os.path.join(benchd,'aux',header+'.pdb')
        fastaf=os.path.join(benchd,'aux',header)
    header,seq=importFastaEntry(open(fastaf,'r'))

    jobID={} #job id of the diferent submitted jobs
    
    #run tasser and spicker (note tas_spk.py located in directory pcd)

    inpdir1=scd+'/out/'+header #contains exp.dat seq.dat rmsinp pair1.dat pair3.dat
    inpdir2=outdir             #contains comb*dat dist*dat chain.dat par.dat

    wt=calc_wt(len(seq))
    jobID['tasser']=job(exe='tas_spk.py',exed=pcd,args='-a '+inpdir1+' -b '+inpdir2+' -c '+outdir+' -d '+CA+' -e '+header).qsub('xt'+jobname,outdir,wallt=wt,mem_limit='890',joblist=joblist)
    if jobID['tasser']=='':Bye('tasser job could not be submitted')
    spktar=outdir+'/'+header+'.out.spicker.tar.bz2'

    continue #insert this here if we only want to run tas_spk.py

    #run dukka and pulchra
    dependencies=[]
    if jobID['tasser']!='COMPLETED': dependencies.append(jobID['tasser'])
    job(exe='dukka_pulchra.py',exed=root,args=' -b '+spktar+' -c '+outdir).qsub(wallt=0.99,mem_limit='40',jobname='xd'+jobname,logname='xd'+jobname+'.log',outdir=outdir,depend=dependencies)
    
sys.exit(0)
