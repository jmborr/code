#!/usr/bin/python
'''
Typical header line for each snapshot in rep1rmsdF.tra
     150    -1854.1     337     337 rmsd=  2.45
'''
import os,sys,re,copy
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName

ih=inpHand('Usage: genChainFromrep1rmsdF.tra.py',
           ' -a _RA_repf rep1rmsdF.tra file',
           ' -b _A_outf output file (def:chain.dat)',
           ).parse(locals(),sys.argv)

if not outf: outf='chain.dat'

repp=open(repf,'r')
buf=''
line=repp.readline()
n_thr=0
while line:
    n_thr+=1 ;print n_thr
    items=line.split() ; L=int(items[0])
    buf2=items[0]+'  '+items[-1]+'\n'
    for i in range(L): buf2+='%4d'%(i+1)+repp.readline()
    buf+=buf2
    line=repp.readline()
repp.close()

outp=open(outf,'w').write(`n_thr`+'\n'+buf)

sys.exit(0)
