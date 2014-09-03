#!/usr/bin/python

'''******************** Pseudocode *****************
1 Input a trajectory file
2 Output a chain.dat file with as many "templates" as snapshots in the trajectory file
'''

import sys,os,re
from utilities.small_utilities import chomp
from random import randint
from inputArgs.inputArgs import inpHand

ih=inpHand(
    'Usage: tra2chain.py ',
    '  -a _AR_rep rep1.tra file',
    '  -b _AR_natf native file name',
    '  -c _A_chain chain.dat file (def: ./chain.dat)',
    )
ih.parse(locals(),sys.argv)

if not chain: chain='chain.dat'
pout=open(chain,'w')
N=1+int( chomp(os.popen('grep -e L= '+rep+'|wc -l').readline()) ) #number of decoys plus native
pout.write(`N`+' hard\n')

#store native file as first decoy
CAlines=os.popen('grep "ATOM" '+natf+'|grep " CA "').readlines()
L=len(CAlines)
pout.write(`L`+' 0.000'+'\n')
for i in range(L):
    prefix='%3d  '%(i+1)
    CAlines[i]=CAlines[i][32:54] #keep coordinates only
    pout.write(prefix+CAlines[i]+'\n')

#store rest of decoys, along with TM score
pattern=re.compile('TM=(0\.\d+)')
pin=open(rep,'r')
for i in range(N-1):
    mt=' %5.3lf'%(1-float(pattern.search(pin.readline()).group(1)))#1-TMscore
    pout.write(`L`+mt+'\n')
    for j in range(1,L+1):
        prefix='%3d '%(j)
        pout.write( prefix+pin.readline() )
pin.close()
pout.close()
