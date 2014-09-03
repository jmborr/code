#!/usr/bin/python
"""

Translate secondary structure (helix H, strand E, coil C) from seq.dat

"""
import os,sys
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand

sblseqdat='alpha:2 beta:4 coil:1 gap:-'

nn2ss={'1':'C', '2':'H', '4':'E'}
ss2nn={'C':'1', 'H':'2', 'E':'4'}

nn2ssII={'1':'-', '2':'H', '4':'E'}
ss2nnII={'-':'1', 'H':'2', 'E':'4'}

nn2ssIII={'1':'L', '2':'H', '4':'B'}
ss2nnIII={'L':'1', 'H':'2', 'B':'4'}

def make_transl(transl):
    '''
    do a translator dictionary from one string
    '''
    tmpd={}
    list=transl.split()
    n=len(list)
    i=0
    while i<n:
        tmpd[ list[i] ] = list[i+1]
        i+=2
    return tmpd

def translateSec(seqdat,translator=nn2ss):    
    #determine type of input for "seqdat" parameter
    pin=None
    if isinstance(seqdat,str):
        try:
            pin=open(seqdat,'r')
        except:
            sys.stderr.write('ERROR secondary.seqdatss.translateSec cannot open '+seqdat+'\n')
    else:
        pin=seqdat
    ss=''
    for line in pin.readlines():
        nn=line.split()[2]
        ss+=translator[nn]
    return ss

if __name__=='__main__':
    inpHand('Usage: pdb2ss.py [options]',
            ' -a _RA_seqdat seq.dat file',
            ' -b _A_outfile output file (def:STDOUT)',
            ' -c __transl translator (def:"1 C 2 H 4 E")',
            ).parse(locals(),sys.argv)
    
    #determine output
    if not outfile: pt=sys.stdout
    else: pt=open(outfile,'w')
    if not transl:transl=nn2ss
    else: transl=make_transl(transl)
    #write secondary structure
    pt.write(translateSec(seqdat,translator=transl)+'\n')
    pt.close()
    sys.exit(0)
