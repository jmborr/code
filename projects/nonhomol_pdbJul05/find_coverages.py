#!/usr/bin/python

import sys,os,re
from utilities.codedir import projectsdir
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand,deglobb
from tmScore.tm_score import gettmCA
from spicker.spickerResultsManager import spickOut

ih=inpHand('Usage: find_coverages.py [options]\nFile summaryTM.all file lists for each header the centroids and their rank. Thus go through each header, pick a particular centroid, and compare to native',
           '  -a _AR_listf file list of headers',
           '  -b _A_rdir root directory (default:projectsdir/createInputsTasser/x/xxxxx/abinitio)',
           '  -c _AR_outf output file containing first line of tm_score.py',
           '  -d __d0 d0 parameter (default d0=8.0)')
ih.parse(locals(),sys.argv)

if not d0: d0=8.0
else:  d0=float(d0)
if not rdir: rdir=projectsdir+'/createInputsTasser/x/xxxxx/abinitio'

#will match "REMARK lseq= 173 rmsd=10.08 tm=0.462 cov=0.66 rmsd_cov= 4.35\n", extract
#"lseq= 173 rmsd=10.08 tm=0.462 cov=0.66 rmsd_cov= 4.35"
headers=chomp(open(listf,'r').readlines())
f2=open(outf,'w')
for header in headers:    
        cdir=deglobb(rdir,header)
        p=spickOut(dir=cdir,nat='CA')
        f2.write(header+' '+p.outFirstLineFromAlignToNatByTM(d0=d0)[6:]+'\n')
sys.exit(0)
