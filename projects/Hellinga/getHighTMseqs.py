#!/usr/bin/python
"""Retrieve sequences having one TASSER-predicted structure with TM>0.4

Steps:
1- Read summary.all, for each entry retrieve ID, TM, closcID,
   and target name
2- Whie TM>0.4, retrieve closc file and do TM-score to target. Identify
   aligned residues and output fasta file.
"""

import os,sys
from utilities.small_utilities import Bye,chomp,junkName
from inputArgs.inputArgs import inpHand
from jobs.job import pastry
from seq.letters import valid_A
from seq.fastaManager import writeFastaToString
from utilities.codedir import codedir

inp=inpHand('Usage: getHighTMseqs.py [options]\n'+locals()['__doc__'],
            ' -a _RA_summf summary.all file',
            ' -b _RA_targetf  target PDB file',
            ' -c _RA_outf output file with compendium of fasta entries',
            ' -d __mintm minimum TM cut-off (def:0.4)',
            )
inp.parse(locals(),sys.argv)

#initialize defaults
if not mintm: mintm=0.4
rootd=os.path.dirname(summf)
buf=''
for line in open(summf,'r').readlines():
    if buf:buf+='\n\n' #separ
    if line[0]=='#': continue #this is a comment line
    items=line.split()
    id=items[0]; target=id[0:5]; closcID=items[17]; ctm=items[20]; tm=float(ctm)
    if tm<mintm: break
    finald=os.path.join(rootd,id)
    closcf='closc'+closcID+'.pdb'
    pastry('tar xf '+os.path.join(finald,id+'.out.ab_spicker.tar')+' '+closcf)
    cmd=codedir+'/f77/TM-score/TMscore '+closcf+' '+targetf+' | grep "::" -B 1'
    (gseq,alig)=chomp(os.popen(cmd,'r').readlines())
    pastry('/bin/rm '+closcf) #clean the extracted closc file
    seq=''
    for i in range(len(gseq)):
        X=gseq[i]
        if X in valid_A:
            if alig[i]==':': seq+=X  #aligned residues in upper case
            else: seq+=X.lower()     #non-aligned residues in lower case
    header=id+' TM='+ctm
    buf+=writeFastaToString(header,seq,n=50)

open(outf,'w').write(buf)
sys.exit(0)
