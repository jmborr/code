#!/usr/bin/python

import sys,os
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye
from letters import returnDict,symbolstr0,translate
from get_seq import get_seq
from Bio.pairwise2 import align
from seq.alignYangManageResults import gappedSeq,alignment

#default scoring matrix in string format
scmtx0='alpha:alpha:1 alpha:beta:-1 alpha:coil:-0.8 beta:beta:1.5 beta:coil:-0.7 coil:coil:0.6'

#   ****************************************

def resolveScoringMatrix(scmtx,symbolstr):
    """prepare scoring matrix of same format as Bio.SubsMat.MatrixInfo.blosum62

    format like this: {('H','H'):1,...}
    scmtr is a string like this: 'alpha:alpha:1 alpha:beta:-1 alpha:coil:-0.8...'
    """
    m={}
    for item in scmtx.split():
        x,y,z=item.split(':')
        d=returnDict(symbolstr)
        m[ (d[x],d[y]) ]=float(z)
    return m

#   ****************************************

#execute as standalone script
if __name__=='__main__':

    #default values for some input parameters
    gapopen0=-1.5
    gapext0=-0.10

    #read command line
    ih=inpHand('Usage: align.py [options]',
               ' -a _R_input1 (either header.fasta or seq.dat or just a sequence',
               ' -b _R_input2 (either header.fasta or seq.dat or just a sequence',
               ' -c __gapopen (gap opening should be negative, def='+`gapopen0`+')',
               ' -d __gapext  (gap extension should be negative, def='+`gapext0`+')',
               ' -e __symbolstr1 symbol string for input1 (def: "alpha:H beta:E coil:C gap:-)',
               ' -f __symbolstr2 symbol string for input2 (def: "alpha:H beta:E coil:C gap:-)',
               ' -g _A_outf outfile for alignment (def:STDOUT)',
               ' -i __symbolstr3 symbol string for alignment (def: same as that of input1)',
               ' -j __algmethod alignment method allowed in Bio.pairwise2 (def:globalds)',
               ' -l __scmtx scoring matrix (def:"'+scmtx0+'"  )',
               )
    ih.parse(locals(),sys.argv)

    from Bio.SubsMat.MatrixInfo import blosum62
    l=blosum62.values();l.sort()
    
    #initialize input parameters
    if not gapopen: gapopen=gapopen0
    else: gapopen=float(gapopen)
    if not gapext: gapext=gapext0
    else: gapext=float(gapext)
    if not symbolstr1: symbolstr1=symbolstr0
    if not symbolstr2: symbolstr2=symbolstr0
    if not symbolstr3: symbolstr3=symbolstr1
    if not algmethod: algmethod='globalds'
    if not scmtx: scmtx=resolveScoringMatrix(scmtx0,symbolstr1)
    else: scmtx=resolveScoringMatrix(scmtx,symbolstr1)
    if not outf: out=sys.stdout
    else: out=open(out,'w')

    #retrieve secondary structure sequences
    seq1=get_seq(input1,symbolstr1)
    seq2=get_seq(input2,symbolstr2)

    #do alignment and output
    if symbolstr2!=symbolstr1: seq2=translate(seq2,symbolstr2,symbolstr1) #reconcile symbol strings
    kargs=[seq1,seq2,1,0,gapopen,gapext]
    #getattr allows us to call algmethod as a method of class pairwise2.align
    results=getattr(align,algmethod)(seq1,seq2,scmtx,gapopen,gapext,one_alignment_only=True)[0]

    if symbolstr3!=symbolstr1: #force to ouptut symbol string
        translate(results[0],symbolstr1,symbolstr3)
        translate(results[1],symbolstr2,symbolstr3)

    out.write( alignment( gappedSeq(results[0]),gappedSeq(results[1]) ).info() )
    
    sys.exit(0)
