#!/usr/bin/python
'''

'''
import os,sys,re
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName

ih=inpHand('Usage: filterByTMcutoff3.py',
           ' -a _R_x number of chunks',
           ' -b _RA_outf output file',
           ' -c _R_prefix prefix'
           ).parse(locals(),sys.argv)

x=int(x)
currd=os.getcwd()
tmpd=currd+'/filterByTMalignDir_'+prefix

#resolve all cross-overs.
lall=[]
for i in range(x-1):
    #merge all crossed-over filtered chunks for a given "i"
    rooti='%02d'%(i)
    l=open(tmpd+'/'+prefix+rooti+'.'+'%02d'%(i+1),'r').readlines()
    lfil=l
    for j in range(i+2,x):
        lfil=[]
        rootj='%02d'%(j)
        fij=tmpd+'/'+prefix+rooti+'.'+rootj
        for s in open(fij,'r').readlines():
            if s in l: lfil.append(s)
        l=lfil #keep only those common structures
    lall+=lfil
lall+=open(tmpd+'/'+prefix+'%02d'%(x-1)+'.1','r').readlines()
open(outf,'w').write(''.join(lall))

#finish with cleanup and exit with success
#os.system('/bin/rm -rf '+tmpd)
sys.exit(0)
