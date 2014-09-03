#!/usr/bin/python

#import some modules to ease the pain
import sys,os,re
from inputArgs.inputArgs import inpHand
from pdbDomains import check_CAs

ih=inpHand(
    'Usage: check_CAs.py ',
    '  -a _A_pdbf pdb file'
    )
ih.parse(locals(),sys.argv)

check_CAs(pdbf);

