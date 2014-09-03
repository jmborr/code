#!/usr/bin/python
'''

'''
import os,sys,re
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp
from jobs.job import job

ih=inpHand('Usage: filterByTMcutoff.py',
           ' -a __tmco TM cut off (def: 0.35)',
           ' -b _RA_inf1 input file with list of structures',
           ' -c _RA_inf2 second input file with list of structures',
           ' -d _RA_outf output file',
           ' -e _A_strd directory containing structures (def:/tmp/jose/structures)'
           ).parse(locals(),sys.argv)


if not tmco:tmco=0.35
else: tmco=float(tmco)

if not strd: strd='/tmp/jose/structures'
strs1=chomp(open(inf1,'r').readlines()) #list of structures
strs2=chomp(open(inf2,'r').readlines()) #list of structures

maxtm=tmco
subset=[]
nn=0
m=0
for pdb1 in strs1:
    n=0
    m+=1
    #print pdb1,m
    fullpdb1=strd+'/'+pdb1 #full path
    #print 'fullpdb1=',fullpdb1
    homolog=0
    for pdb2 in strs2:
        fullpdb2=strd+'/'+pdb2 #full path
        tm=float(os.popen('TMalign '+fullpdb1+' '+fullpdb2+'|grep "TM-score="|cut -d\'=\' -f 4|cut -d\',\' -f 1').readline())
        n+=1
        #print pdb1,pdb2,tm
        if tm>maxtm: #we found a structural homologous in subset
            homolog=1
            break
    nn+=n
    #print n
    if not homolog:
        subset.append(pdb1)

#print only the filtered subset of inf1
open(outf,'w').write('\n'.join(subset)+'\n')
    
sys.exit(0)
