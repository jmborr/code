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
scd='/gpfs1/scratch/jose/spliceVariants'
pcd='/gpfs1/active/jose/code/projects/spliceVariants'
root='/library/jose/code/python/combo_jobs/predictor'
#benchd='/gpfs1/archive/adrian/dat/pdb_mar07'

#Some options
wt='0.249' #less than siz hours
fastaf0=scd+'/input/xxxxx'

for header in chomp(open(scd+'/input/list','r').readlines()):
    outdir=os.path.join(scd,'out',header[0:6],header[6:])
    jobname=header #because multiple jobs have same header
    fastaf=deglobb(fastaf0,header)

    inpdir1=scd+'/out/'+header #contains exp.dat seq.dat rmsinp pair1.dat pair3.dat
    inpdir2=outdir             #contains comb*dat dist*dat chain.dat par.dat

    #run spicker
    intasserf=outdir+'/'+header+'.in.tasser.tar.bz2'
    tasserf=outdir+'/'+header+'.in.spicker.tar.bz2'
    expf=outdir+'/exp.dat'
    spkout='closc* combo* loop*.dat rst*.dat' #output from spicker
    spktbz2=header+'.out.spicker.tar.bz2'      #compressed output from spicker
    spkx='/gpfs1/active/jose/code/projects/spliceVariants/benchmarkByTMalign/spk_rp'

    cmd='generic_job_list.py -c '+outdir+\
    ' -n spk'+\
    ' -b '+jobname+\
    ' -d '+wt+\
    ' -k 990'+\
    ' -f "/bin/cp '+tasserf+' '+intasserf+' '+expf+' '+spkx+' ."'+\
    ' -g "/bin/tar jxf '+intasserf+' chain.dat rmsinp seq.dat && /bin/tar jxf '+tasserf+' && ./spk_rp && tar jcf '+spktbz2+' '+spkout+'"'+\
    ' -i "/bin/mv '+spktbz2+' '+outdir+'"'+\
    ' -q no'
    #print outdir,cmd

    os.system(cmd)
    
sys.exit(0)
