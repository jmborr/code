#!/usr/bin/python

import os,sys
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir

from seq.fastaManager import importFastaEntry

inpHand('Usage: tas_spk.py [options]',
        '  -a _RA_inpdir1 directory containing exp.dat seq.dat rmsinp pair1.dat pair3.dat',
        '  -b _RA_inpdir2 directory containing comb*dat dist*dat chain.dat par.dat',
        '  -c _RA_outdir output directory',
        '  -d _RA_CA calpha file',
        '  -e _R_header header',
        '  -f _A_cas tasser executable',
        '  -g _A_spicker spicker executable',
        ).parse(locals(),sys.argv)

#Bring input
currd=os.getcwd()

os.chdir(inpdir1)
pastry('/bin/cp exp.dat seq.dat rmsinp pair1.dat pair3.dat '+currd)

os.chdir(inpdir2)
pastry('/bin/cp comb*dat dist*dat chain.dat par.dat '+currd)

os.chdir('/gpfs1/active/jose/code/f77/tasser/tasser.1.0')
pastry('/bin/cp tra.in original/in.dd '+currd)

os.chdir(currd)
pastry('tar jxf /gpfs1/active/jose/code/f77/tasser/tasser.1.0/common.tbz2 -C .')

pastry('grep "ATOM" '+CA+'|grep " CA " > CA')


#run cas
if not cas: cas=codedir+'/f77/tasser/tasser.1.0/learning/cas'
pastry('/bin/cp '+cas+' .')
pastry('./'+os.path.basename(cas))
pastry('tar jcf '+header+'.out.tasser.tar.bz2 rep*tra tra.in')
pastry('/bin/mv '+header+'.out.tasser.tar.bz2 '+outdir)

#run spicker (spicker version generating loop.dat)
if not spicker: spicker=codedir+'/projects/spliceVariants/benchmarkByTMalign/spk_rp'
pastry('/bin/cp '+spicker+' .')
pastry('./spk_rp')
pastry('tar jcf '+header+'.out.spicker.tar.bz2 combo* closc* str.txt rst*.dat loop*.dat')
pastry('/bin/mv '+header+'.out.spicker.tar.bz2 '+outdir)

sys.exit(0)
