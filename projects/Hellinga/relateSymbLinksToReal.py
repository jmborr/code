#!/usr/bin/python
'''

'''
import sys,os
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb,extractHeader

ih=inpHand('Usage: relateSymbLinksToReal.py [options]\nExample:./relateSymbLinksToReal.py -a /gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/thioredoxin/filtered/list -b /gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/thioredoxin/filtered',
           ' -a _AR_filtf list of symbolic links',
           ' -b _AR_filtd directory containing symbolic links (filtd)',
           ' -c _A_s2r symbolic to real file (def:filtd/symb2real.dat)',
           ' -d _A_r2s real to symbolic file (def:filtd/real2symb.dat)'
           )
ih.parse(locals(),sys.argv)

if not s2r: s2r=filtd+'/symb2real.dat'
if not r2s: r2s=filtd+'/real2symb2.dat'
bufs2r=''
bufr2s=''
for symb in chomp(open(filtf,'r').readlines()):
    str=chomp(os.popen('ls -l '+filtd+'/'+symb).readline()).split('/')[-1] #name of structure
    bufs2r+=symb+' '+str+'\n'
    bufr2s+=str+' '+symb+'\n'

open(s2r,'w').write(bufs2r)
open(r2s,'w').write(bufr2s)

sys.exit(0)
