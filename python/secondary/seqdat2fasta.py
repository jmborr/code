#!/usr/bin/python
"""

Obtain fasta file from the typical seq.dat file

"""
import os,sys
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand
from seq.letters import three2one
from seq.fastaManager import writeFasta

def extractSeq(seqdat):
    seq=''
    inpp=open(seqdat,'r')
    line=inpp.readline()
    while line:
        seq+=three2one[ line.split()[1] ]
        line=inpp.readline()
    return seq
    
if __name__=='__main__':
    inpHand('Usage: seqdat2fasta.py [options]',
            ' -a _RA_seqdat seq.dat-like file',
            ' -b _R_header header of the fasta file',
            ' -c _A_outf output file (def:STDOUT)',
            ).parse(locals(),sys.argv)
    if not outf: outp=sys.stdout
    else: outp=open(outf,'w')
    
    writeFasta(header,extractSeq(seqdat),filename=outp)

