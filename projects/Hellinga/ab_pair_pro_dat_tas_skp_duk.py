#!/usr/bin/python

import os,sys
from utilities.small_utilities import Bye,chomp,junkName
from inputArgs.inputArgs import inpHand,deglobb
from utilities.codedir import codedir
from jobs.job import pastry

inp=inpHand('Usage: ab_pro_dat_tas_spk_duk.py [options]\ncalculate a normal tasser run\n',
            ' -a _RA_fastaf sequence file in fasta format',
            ' -b _RA_blastout blastout tarred results from blast.py script',
            ' -c _RA_outdir output directory',
            ' -d _R_header single header (in place of list of headers)',
            ' -g _RA_expdatf exp.dat file',
            ' -i __co sequence cutoff for mkpair(def:1.0)',
            ' -j _RA_inddf in.dd file',
            ' -k _RA_seqdatf seq.dat file',
            ' -l _A_chaindatf chain.dat file (def:none)',
            ' -m _A_CA native file(def:none)',
            )
inp.parse(locals(),sys.argv)

if not co: co='1.0'
currd=os.getcwd()
root=os.path.join(codedir,'python/combo_jobs/predictor')
rootH='/gpfs1/active/jose/code/projects/Hellinga'


#print '#run inprosp.py only if inprosp does not exists'
inprosp=outdir+'/'+header+'.in.prospector.tar.bz2'
if not os.path.exists(inprosp):
    wd=os.path.join(currd,junkName()) #create a subdirectory
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+rootH+'/inprosp.py . ')
    cmd='./inprosp.py -a '+fastaf+' -b '+blastout+' -c '+outdir+' -d '+seqdatf
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

#print '#run pair.py only if pairtar does not exists'
wd=os.path.join(currd,junkName()) #create a subdirectory
pairtar=outdir+'/'+header+'.pair.tar.bz2'
if not os.path.exists(pairtar):
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+root+'/pair.py . ')
    cmd='./pair.py -a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+co+' -e '+outdir
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)


#mkrmsinp and mkdat run
dat_abtar=outdir+'/'+header+'.dat_ab.tar.bz2'
if not os.path.exists(dat_abtar):
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+rootH+'/dat_ab.py . ')
    cmd='./dat_ab.py -a '+fastaf+' -b '+outdir+' -c '+seqdatf
    if chaindatf: cmd+=' -d '+chaindatf
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

#run tasser and spicker
ab_spktar=outdir+'/'+header+'.out.ab_spicker.tar.bz2'
if not os.path.exists(ab_spktar):
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+rootH+'/tas_spk_ab.py . ')
    args='-a '+fastaf+' -b '+dat_abtar+' -c '+pairtar+' -d '+inprosp+' -e '+expdatf+' -f '+outdir+' -i '+inddf
    if CA: args+=' -j '+CA
    cmd='./tas_spk_ab.py '+args
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

#run dukka and pulchra
abduktar=outdir+'/'+header+'abdukka.tar.bz2'
if not os.path.exists(abduktar):
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+root+'/dukka_pulchra.py . ')
    cmd='./dukka_pulchra.py -b '+ab_spktar+' -c '+outdir+' -d '+header+' -f '+header+'abdukka.tar.bz2'
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

sys.exit(0)
