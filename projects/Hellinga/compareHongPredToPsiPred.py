#!/usr/bin/python
import sys,os,re
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import addAbsPath,inpHand
from secondary.seqdat2ss import translateSec,nn2ssIII


def compare(ss1,ss2):
    nn={'A':0,'H':0,'B':0,'L':0}
    id={'A':0.0,'H':0.0,'B':0.0,'L':0.0}

    N=len(ss2)
    for i in range(N):
        x=ss2[i]
        nn['A']+=1
        nn[x]+=1
        if ss1[i]==x:
            id[x]+=1
            id['A']+=1

    for key,value in nn.items():
        if value==0: nn[key]=1
        id[key]=id[key]/nn[key]
    return '%4.2lf %4.2lf %4.2lf %4.2lf'%(id['A'],id['H'],id['B'],id['L'])

ih=inpHand('Usage: compareHongPredToPsiPred.py [options]',
           ' -a _RA_seqdatf seq.dat file',
           ' -b _RA_hongf hongpred.dat file',
           ' -c _RA_origf original seq.dat file')
ih.parse(locals(),sys.argv)

ord=translateSec(origf,translator=nn2ssIII)
sqd=translateSec(seqdatf,translator=nn2ssIII)
hpd=open(hongf,'r').readline().strip()

print 'hongpred vs  psipred='+compare(hpd,sqd)
print ' psipred vs original='+compare(sqd,ord)
print 'hongpred vs original='+compare(hpd,ord)

sys.exit(0)
