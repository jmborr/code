#!/usr/bin/python

import os,sys
from utilities.small_utilities import Bye,chomp,junkName
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import pastry

inp=inpHand('Usage: pro_dat_tas_spk_duk.py [options]\ncalculate a normal tasser run\n',
            ' -a _RA_fastaf sequence file in fasta format',
            ' -b _RA_inprosp tarred file with input for prospector',
            ' -c _RA_outdir output directory',
            ' -d _R_header single header (in place of list of headers)',
            ' -e __co2 remove homologs above sequence cutoff for prospector and mkdat (def:1.0)',
            ' -f _RA_blastout tarred file with previously run blast',
            ' -g _RA_expdatf exp.dat file',
            ' -i _RA_pairtar contains pair1.dat, pair3.dat',
            ' -j _RA_inddf in.dd file'
            )
inp.parse(locals(),sys.argv)

currd=os.getcwd()
root='/library/jose/code/python/combo_jobs/predictor'
rootH='/gpfs1/active/jose/code/projects/Hellinga'

#prospector run
prosptar=outdir+'/'+header+'.out.prospector.tar'
if not os.path.exists(prosptar):
    sys.stdout.write('Running prosp.py')
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+root+'/prosp.py . ')
    cmd='./prosp.py -a '+fastaf+' -b '+inprosp+' -c '+outdir+' -d '+co2
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

#mkrmsinp and mkdat run
dattar=outdir+'/'+header+'.dat.tar'
if not os.path.exists(dattar):
    sys.stdout.write('Running dat.py')
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+rootH+'/dat.py . ')
    cmd='./dat.py -a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+prosptar+' -e '+outdir+' -f '+co2
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

#run tasser and spicker
spktar=outdir+'/'+header+'.out.spicker.tar'
if not os.path.exists(spktar):
    sys.stdout.write('Running tas_skp.py')
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+rootH+'/tas_spk.py . ')
    cmd='./tas_spk.py -a '+fastaf+' -b '+dattar+' -c '+pairtar+' -d '+inprosp+' -e '+expdatf+' -f '+outdir+' -g '+inddf
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

#run dukka and pulchra
duktar=outdir+'/'+header+'dukka.tar'
if not os.path.exists(duktar):
    sys.stdout.write('Running dukka_pulchra.py')
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+root+'/dukka_pulchra.py . ')
    cmd='./dukka_pulchra.py -b '+spktar+' -c '+outdir+' -d '+header+' -f '+header+'dukka.tar'
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

sys.exit(0)
