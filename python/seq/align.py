#!/usr/bin/python

import sys,os
from inputArgs.inputArgs import inpHand
from seq.get_seq import get_seq
from seq.alignYangManageResults import gappedSeq,alignment

def do_alignment(seq1,seq2,gapopen=-11,gapext=-1):
    from Bio.SubsMat.MatrixInfo import blosum62
    from Bio.pairwise2 import align
    return align.globalds(seq1,seq2,blosum62,gapopen,gapext,one_alignment_only=True)[0]

def getYangManageResultsAlignmentFormat(seq1,seq2,gapopen=-11,gapext=-1):
    alignments=do_alignment(seq1,seq2,gapopen,gapext)
    gs1=gappedSeq(alignments[0],gapsymbol='-')
    gs2=gappedSeq(alignments[1],gapsymbol='-')
    return alignment(gs1,gs2)
    
def seqId(alignments):
    id=0
    query=alignments[0]
    target=alignments[1]
    begin=alignments[3]
    end=alignments[4]
    for i in range(begin,end):
        if query[i]==target[i] and query[i]!='-':id+=1
    return float(id)/(end-begin)

    
#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: align.py [options]',
               '  -a _AR_input1 (either xxx.fasta or xxx.seq or xxx.pdb or a pdbheader(xxxxx) or just a protein sequence',
               '  -b _AR_input2 (either xxx.fasta or xxx.seq or xxx.pdb or a pdbheader(xxxxx)',
               '  -c  __gapopen (gap opening should be negative, def=-11)',
               '  -d  __gapext  (gap extension should be negative, def=-1)',
               '  -e __outname outfile name',
               '  -o _A_outd output directory(default:current dir)')
    ih.parse(locals(),sys.argv)
    
    if not input1 and not input2:
        sys.stderr.write('ERROR: -a and -b needed')
        sys.exit(1)
    if not gapopen: gapopen=-11
    if not gapext: gapext=-1
    seq1=get_seq(geninput=input1)
    seq2=get_seq(geninput=input2)

    alignments=do_alignment(seq1,seq2,gapopen=gapopen,gapext=gapext)
    #print alignments
    #print seqId(alignments)
    sys.exit(0)
    
