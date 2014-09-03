#!/usr/bin/python


import sys,os,re
from utilities.small_utilities import chomp,junkName,Bye
from inputArgs.inputArgs import inpHand,deglobb

alph2=['(H+)','(E+)','^(-+)','(-+)$','H(-+)H','E(-+)E','H(-+)E','E(-+)H']

if __name__=='__main__':   
    inpHand('Usage: countElements.py\n',
            ' -a _AR_seqdat seq.dat file'
            ).parse(locals(),sys.argv)

xss={'2':'H', '4':'E', '1':'C'}
alph2=['(H+)','(E+)','^(-+)','(-+)$','H(-+)H','E(-+)E','H(-+)E','E(-+)H']
px={}
for x in alph2:
    px[x]=re.compile(x) #regular expression to detect contigous stretches

lines=open(seqdat,'r').readlines()
ss=''
for line in lines: ss+=xss[ line.split()[2] ]
print len( px['(H+)'].findall(ss) )
