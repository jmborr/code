#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.small_utilities import three2one #three2one dict for aa codes
from utilities.small_utilities import chomp

ih=inpHand('Usage: mkfasta.py [options]',
           ' -a _R_header header for the first line of the fasta file',
           ' -b _A_seqdatf seq.dat-like file',
           ' -c _A_outf outpuf file (def: STDOUT)'
           )
ih.parse(locals(),sys.argv)

if not outf: out=sys.stdout
else: out=open(outf,'w')

try: seqdat=open(seqdatf,'r')
except IOError: ih.abort(' !!! seq.dat like file '+seqdatf+' non-existent')
outl='>'+header+'\n'
buf=''
for line in chomp(seqdat.readlines()):
    buf+=three2one[line.split()[1]]
while buf:
    outl=outl+buf[0:60]+'\n'
    buf=buf[60:]

out.write(outl) #handle to file        

sys.exit(0)
