#!/usr/bin/python

''' Program to assign a secondary structure to one of the decoys. We
assign as the original secondary structure, except if the secondary
predicted from the amino acid coordinates says the residue should be
assigned as helical '''

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp,junkName,Bye
from secondary.pdb2ss import outputSec
from secondary.seqdat2ss import translateSec,nn2ssIII

ih=inpHand('Usage: genSeqForDecoys.py required [options] '+__doc__,
           ' -a _RA_origf original secondary sequence file seq.dat',
           ' -b _RA_pdbf pdb structure',
           )
ih.parse(locals(),sys.argv)

orig0=translateSec(origf,translator=nn2ssIII)
out=outputSec(pdbf)
orig=''
for i in range(len(orig0)):
    if out[i]=='H': orig+='H'
    else: orig+=orig0[i]
print orig



