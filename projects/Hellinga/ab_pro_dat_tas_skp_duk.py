#!/usr/bin/python

import os,sys
from utilities.small_utilities import Bye,chomp,junkName
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import pastry

inp=inpHand('Usage: ab_pro_dat_tas_spk_duk.py [options]\ncalculate a normal tasser run\n',
            ' -a _RA_fastaf sequence file in fasta format',
            ' -b _RA_inprosp tarred file with input for prospector',
            ' -c _RA_outdir output directory',
            ' -d _R_header single header (in place of list of headers)',
            ' -g _RA_expdatf exp.dat file',
            ' -i _RA_pairtar pairtar contains pair1.dat, pair3.dat',
            ' -j _RA_inddf in.dd file',
            ' -k _RA_seqdatf seq.dat file',
            ' -l _A_chaindatf chain.dat file (def:none)'
            )
inp.parse(locals(),sys.argv)

currd=os.getcwd()
root='/library/jose/code/python/combo_jobs/predictor'
rootH='/gpfs1/active/jose/code/projects/Hellinga'

#mkrmsinp and mkdat run
dat_abtar=outdir+'/'+header+'.dat_ab.tar'
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
ab_spktar=outdir+'/'+header+'.out.ab_spicker.tar'
if not os.path.exists(ab_spktar):
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+rootH+'/tas_spk_ab.py . ')
    cmd='./tas_spk_ab.py -a '+fastaf+' -b '+dat_abtar+' -c '+pairtar+' -d '+inprosp+' -e '+expdatf+' -f '+outdir+' -i '+inddf
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

#run dukka and pulchra
abduktar=outdir+'/'+header+'abdukka.tar'
if not os.path.exists(abduktar):
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+root+'/dukka_pulchra.py . ')
    cmd='./dukka_pulchra.py -b '+ab_spktar+' -c '+outdir+' -d '+header+' -f '+header+'abdukka.tar'
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

sys.exit(0)
