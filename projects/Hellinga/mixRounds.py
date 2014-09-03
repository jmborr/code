#!/usr/bin/python
"""
Examples:

I a0216 a0216_0145 1 0.16 0.02702 14.80 0.26 1 0.16 0.02702 14.80 0.26 5 8 0.05 0.01015 7.66 4 2 0.08 0.01308 0.48

II a0216 a0216_0735 1 0.15 0.02084 10.36 0.41 1 0.15 0.02084 10.36 0.41 9 10 0.06 0.01042 6.42 9 10 0.06 0.01042 0.50

III a0512 0251_050 2 0.15 0.02131 14.77 0.28 2 0.15 0.02629 13.15 0.36 6 4 0.06 0.01084 5.30 6 4 0.06 0.01084 0.57
"""
import os,sys
from utilities.small_utilities import Bye,chomp,junkName
from utilities.codedir import scratchdir
from inputArgs.inputArgs import inpHand


inp=inpHand('Usage: mixRounds.py [options]\nextract sequence\n',
            ' -a _R_line line from one summary.all file in I_II_III directory',
            )
inp.parse(locals(),sys.argv)

root=scratchdir+'/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum'

round,parent,id=line.split()[0:3]

if round=='I': round=''
print root+'/hongpred_high_seqid_high_helical_content'+round+'/'+parent+'/'+id+'/'+id #fasta file

sys.exit(0)
