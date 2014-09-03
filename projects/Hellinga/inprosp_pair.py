#!/usr/bin/python

import os,sys
from utilities.small_utilities import Bye,chomp,junkName
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import pastry

inp=inpHand('Usage: inprosp_pairb.py [options]\ncalculate ab initio and regular structure prediction\n',
            ' -a _RA_fastaf sequence file in fasta format',
            ' -b _RA_outdir output directory',
            ' -c _R_header single header (in place of list of headers)',
            ' -d __co sequence cutoff for mkpair(def:1.0)',
            ' -e _AR_seqdatf pass seq.dat file so that we don\'t calculate it',
            ' -f _RA_blastout blastout tarred results from blast.py script',
            )
inp.parse(locals(),sys.argv)


if not co:co='1.0'

root='/library/jose/code/python/combo_jobs/predictor'
rootH='/gpfs1/active/jose/code/projects/Hellinga'

currd=os.getcwd()

#run inprosp.py only if inprosp does not exists
inprosp=outdir+'/'+header+'.in.prospector.tar'
if not os.path.exists(inprosp):
    wd=os.path.join(currd,junkName()) #create a subdirectory
    os.system('/bin/mkdir -p '+wd)
    os.chdir(wd)
    os.system('/bin/cp '+rootH+'/inprosp.py . ')
    cmd='./inprosp.py -a '+fastaf+' -b '+blastout+' -c '+outdir+' -d '+seqdatf
    pastry(cmd)
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd)

#run pair.py
wd=os.path.join(currd,junkName()) #create a subdirectory
os.system('/bin/mkdir -p '+wd)
os.chdir(wd)
os.system('/bin/cp '+root+'/pair.py . ')
cmd='./pair.py -a '+fastaf+' -b '+blastout+' -c '+inprosp+' -d '+co+' -e '+outdir
pastry(cmd)
os.chdir(currd)
os.system('/bin/rm -rf '+wd)

sys.exit(0)
