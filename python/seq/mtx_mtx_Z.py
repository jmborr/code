#!/usr/bin/python

"""
calculate profile similarity between two .mtx files by Z-score

"""

import sys,re,os
from inputArgs.inputArgs import inpHand
from seq.mtx import mtx1_mtx2_Z_rank
from seq.align import getYangManageResultsAlignmentFormat
#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: mtx_mtx_Z.py [options]',
               '  -a _RA_mtx1 first .mtx file',
               '  -b _RA_mtx2 second .mtx file',
               )
    ih.parse(locals(),sys.argv)

    #do global alignment for the two sequences, return pairs of
    #indexes [ (i1,j1), (i2,j2), ...,(in,jn)] where residue in i1 of
    #seq1 identical to residue j1 of seq2.    
    seq1=open(mtx1,'r').readlines()[1].strip()
    seq2=open(mtx2,'r').readlines()[1].strip()
    pairs=getYangManageResultsAlignmentFormat(seq1,seq2).listIndexPairsIdentical()

    print mtx1_mtx2_Z_rank(mtx1,mtx2,alignment=pairs)
    sys.exit(0)
               
