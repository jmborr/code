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


joblist=genJobList()

currd=os.getcwd()
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
pcd='/gpfs1/active/jose/code/projects/spliceVariants/benchmarkByTMalign'
root='/library/jose/code/python/combo_jobs/predictor'
benchd='/gpfs1/scratch/adrian/dat/pdb_mar07'

#Some options
wt='0.249' #less than siz hours
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
#for i in range(10):
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

    #run spicker    
    tasserf=outdir+'/'+header+'.out.tasser.tar.bz2'
    chaindatf=outdir+'/chain.dat'
    expf=scd+'/out/'+header+'/exp.dat'
    rmsinpf=scd+'/out/'+header+'/rmsinp'
    seqdatf=scd+'/out/'+header+'/seq.dat'
    spkout='closc* combo* loop*.dat rst*.dat' #output from spicker
    spktbz2=header+'.out.spicker.tar.bz2'      #compressed output from spicker
    spkx='/gpfs1/active/jose/code/projects/spliceVariants/benchmarkByTMalign/spk_rp'

    cmd='generic_job_list.py -c '+outdir+\
    ' -n spk'+\
    ' -b '+jobname+\
    ' -d '+wt+\
    ' -k 990'+\
    ' -f "/bin/cp '+tasserf+' '+chaindatf+' '+expf+' '+rmsinpf+' '+seqdatf+' '+spkx+' ."'+\
    ' -g "/bin/tar jxf '+tasserf+' && /bin/cp '+CA+' ./CA && ./spk_rp && tar jcf '+spktbz2+' '+spkout+'"'+\
    ' -i "/bin/mv '+spktbz2+' '+outdir+'"'+\
    ' -q no'
    print outdir#,cmd

    os.system(cmd)
    
sys.exit(0)
