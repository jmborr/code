#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq

ih=inpHand('Usage: sleepy.py',
           ' -h __help this script runs (possibly) prospector, tasser and spicker, then use the move '
           )
ih.parse(locals(),sys.argv)

os.system('sleep 600') #sleep for 10 minutes
