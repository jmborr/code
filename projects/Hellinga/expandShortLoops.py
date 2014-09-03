#!/usr/bin/python
'''

'''
import os,sys,re,copy
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName
from secondary.seqdat2ss import translateSec
from random import randint

ih=inpHand('Usage: expandShortLoops.py [options]',
           ' -a _AR_seqdatf seq.dat file',
           ' -b __minLoopL minimum loop length (def: 3)',
           ' -c _A_outf output file (def:STDOUT)'
           )
ih.parse(locals(),sys.argv)

#default arguments
if not minloopL: minloopL=3
else: minloopL=int(minloopL)
if not outf: outp=sys.stdout
else: outp=open(outf,'w')

ss=translateSec(seqdatf) #;Bye('ss='+ss) #translate to a single sequence of 'C', 'H', 'E'
Ntarget=len(ss)
loops=re.compile('C+').findall(ss) #list of loops
nloops=len(loops)
#initialize pointers to loops according to length
smaller=[] #loops smaller than minloopL
equal=[]
bigger=[]
for i in range(nloops):
    l=len(loops[i])
    pair={'order':i,'length':l}
    if l<minloopL: smaller.append(pair)
    elif l==minloopL: equal.append(pair)
    else: bigger.append(pair)

#expand short loops at the expense of the longer loops
nsmall=len(smaller)
nbigger=len(bigger)
while nsmall and nbigger:
    n=randint(0,nsmall-1) #randomly pick one small loop
    smaller[n]['length']+=1 #increase its length by one
    if smaller[n]['length']==minloopL: #maybe the length of the loop became minloopL
        equal.append( smaller[n] )
        del smaller[n]
        nsmall-=1 #one less small loop
    m=randint(0,nbigger-1)
    bigger[m]['length']-=1
    if bigger[m]['length']==minloopL:
        equal.append( bigger[m] )
        del bigger[m]
        nbigger-=1

#join all loops
pairs=smaller+equal+bigger
#junk=0
#for pair in pairs: junk+=pair['length']
#Bye(`junk`)

#output new seq.dat file
buf='' #buffer
m=0  #running index on ss
mm=1 #running index on new seq.dat
iloop=0
ssnew=''
while m<Ntarget:
    if ss[m]=='C': #print iloop
        for pair in pairs:
            if pair['order']==iloop:
                for i in range(pair['length']): #length of the loop
                    ssnew+='C'
                    buf+='%5d   VAL%5d   10\n'%(mm,1)
                    mm+=1
                break #do not search remaining pairs
        while m<Ntarget and ss[m]=='C': m+=1 #do not pay attention to old loop
        iloop+=1
    elif ss[m]=='H':
        buf+='%5d   VAL%5d   10\n'%(mm,2)
        ssnew+='H'
        mm+=1
        m+=1
    else: #strand
        buf+='%5d   VAL%5d   10\n'%(mm,4)
        ssnew+='E'
        mm+=1
        m+=1

outp.write(buf)
#sys.stdout.write(ss+'\n'+ssnew+'\n')
sys.exit(0)
