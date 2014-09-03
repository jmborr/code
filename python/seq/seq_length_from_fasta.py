#!/usr/bin/python

"""Output sequence length from a fasta file"""

import sys,os
from inputArgs.inputArgs import inpHand,deglobb
from seq.fastaManager import importFastaEntry
from utilities.small_utilities import chomp

def getLength(fastaf):
    """Output sequence length from a fasta file"""
    header,seq=importFastaEntry( open(fastaf,'r') )
    return len(seq)
    
#execute as standalone program
if __name__=='__main__':

    ih=inpHand('Usage: seq_length_from_fasta.py  [options]',
               '  -a _RA_fastaf (blogged) fasta file',
               '  -b _A_listf file with list of headers.',
               )
    ih.parse(locals(),sys.argv)

    #prepare all file names
    listh=[]
    if listf:
        for header in chomp(open(listf,'r').readlines()): listh.append( deglobb(fastaf,header) )
    else:
        listh=[fastaf,]
    #process each fasta file
    for fastaf in listh:
        print fastaf,getLength(fastaf)
        
        
