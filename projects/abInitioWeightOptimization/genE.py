#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp
from utilities.codedir import codedir
from inputArgs.inputArgs import inpHand
from jobs.job import pastry

#argument options
ih=inpHand(
    'Usage: genE.py ',
    '  -a _R_header header'
    ).parse(locals(),sys.argv)

#initialize some variables
dir='/gpfs1/scratch/jose/abInitioWeightOptimization/100to199/'+header[1]+'/'+header
natf=header+'.pdb'

#bring input to working directory
wd=os.getcwd() #working directory
os.system('/bin/cp /gpfs1/archive/adrian/dat/pdb_may06.apr06.aux/CA/'+natf+ ' .')
os.chdir(dir)
pastry('/bin/cp rep1.tra.bz2 exp.dat pair1.dat pair3.dat rmsinp seq.dat '+wd)
os.chdir(codedir+'/f77/tasser/tasser.1.0')
os.system('/bin/cp report_E/reportEab abinitio/in.dd '+wd)
os.system('/bin/tar jxf common.tbz2 -C '+wd)
os.chdir(codedir+'/projects/abInitioWeightOptimization')
os.system('/bin/cp rgFilter.py filterByRmsd.x tra2chain.py gterms.py '+wd)

#run necessary jobs
os.chdir(wd)
os.system('bunzip2 rep1.tra.bz2')
#Once we have the decoys in rep1.tra, we filter out all structures with
#radius of gyration bigger than 1.5*native-radius-of-gyration
pastry('./rgFilter.py -a rep1.tra -b '+natf) #-> rep1RgF.tra 
#Once we have the filtered trajectory rep1RgF.tra, we reduce redundancy
#by not allowing two structures with less than 1Angstrom RMSD. In
#addition, we don't keep more than 10,000 structures with TM-score
#to native in between [x,x+01], where x=0.4, 0.5, 0.6, 0.7, 0.8
pastry('./filterByRmsd.x') #-> rep1rmsdF.tra, tmHisto.dat
#rep1rmsdF.tra is the set of filtered structures, tmHisto.dat is the
#histogram of number of structures in the aforementioned TM-score bins.
pastry('./tra2chain.py -a rep1rmsdF.tra -b '+natf) #-> chain.dat
#Once we have the filtered trajectory rep1rmsdF.tra, we transform the
#format so that we cand feed it to a modified TASSER code that will
#output the energy terms for each decoy.
pastry('./reportEab') #-> energy_terms.dat
#Once we have the energy terms of each decoy, do different averages over all decoys
pastry('./gterms.py -a energy_terms.dat -b energy_averages.dat')

#bring output to dir
os.system('bzip2 rep1RgF.tra rep1rmsdF.tra chain.dat energy_terms.dat')
os.system('/bin/mv rep1RgF.tra.bz2 rep1rmsdF.tra.bz2 chain.dat.bz2 energy_terms.dat.bz2 tmHisto.dat energy_averages.dat '+natf+' '+dir)

sys.exit(0)
