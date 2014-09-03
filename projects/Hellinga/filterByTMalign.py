#!/usr/bin/python
'''

'''
import os,sys,re
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp
from utilities.codedir import codedir
from jobs.job import job

ih=inpHand('Usage: filterByTMcutoff.py',
           ' -a __tmco TM cut off (def: 0.35)',
           ' -b _RA_inf input file with list of structures',
           ' -c _RA_outf output file',
           ' -d _A_strd directory containing the structures (def:/tmp/jose/structures)'
           ).parse(locals(),sys.argv)


if not tmco:tmco=0.35
else: tmco=float(tmco)

if not strd:strd='/tmp/jose/structures'

strs=chomp(open(inf,'r').readlines()) #read the list of structures

#initialize subset, tms, and list with 0.85 TM cutoff
tms={}

#within parenthesis formula of Betancourt & Skolnick, Biopolymers 59, 305-309 (2001)
os.system('/bin/cp '+codedir+'/f77/TM-align/TMalign .')
maxtm=tmco
subset=[]
nn=0
m=0
for pdb1 in strs:
    n=0
    m+=1
    #print pdb1,m
    fullpdb1=strd+'/'+pdb1 #full path
    #print 'fullpdb1=',fullpdb1
    homolog=0
    tmsTemp={pdb1:{}}
    for pdb2 in subset:
        fullpdb2=strd+'/'+pdb2 #full path
        tm=float(os.popen('TMalign '+fullpdb1+' '+fullpdb2+'|grep "TM-score="|cut -d\'=\' -f 4|cut -d\',\' -f 1').readline())
        n+=1
        #print pdb1,pdb2,tm
        if tm>maxtm: #we found a structural homologous in subset
            homolog=1
            break
        tmsTemp[pdb1][pdb2]=tm
        tmsTemp[pdb2]={}
        tmsTemp[pdb2][pdb1]=tm
    nn+=n
    #print n
    if not homolog:
        subset.append(pdb1)
        for key in tmsTemp.keys():
            if key not in tms.keys():
                tms[key]=tmsTemp[key]
            else:
                tms[key].update(tmsTemp[key]) #merge the two dictionaries

#print nn
open(outf,'w').write('\n'.join(subset)+'\n')

sys.exit(0)
