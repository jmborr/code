#!/usr/bin/python

import os,sys,re
from utilities.small_utilities import chomp,junkName,Bye
from inputArgs.inputArgs import addAbsPath,inpHand
from tasser.repTraReader import repTraReader

ih=inpHand(
    'Usage: reptra2pdbs.py [options] ',
    '  -a _AR_rep rep*.tra file',
    '  -b _AR_templpdb template pdb file',
    '  -c _A_outd output directory (def:./tmp)',
    '  -d __prefix prefix for file names (def:a)'
    )
ih.parse(locals(),sys.argv)
if not outd: outd='./tmp'
if not prefix: prefix='a'

rTr=repTraReader(rep)
rTr. write2pdb(dir=outd,prefix=prefix,strf=templpdb)

sys.exit(0)
