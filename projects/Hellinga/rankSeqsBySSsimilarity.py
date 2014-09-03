#!/usr/bin/python
import os,sys,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName
from secondary.seqdat2ss import translateSec

halfmaxspan={'H':8,'E':6,'C':10} #half of maximum span of a secondary element

def read_weights(seqdat):
    cmd='cat '+seqdat+' | tr -s \' \' | cut -d\' \' -f 5' #;Bye(cmd)
    ws=chomp(os.popen(cmd).readlines())
    for i in range(len(ws)): ws[i]=float(ws[i])/10
    return ws
    
ih=inpHand('Usage: rankSeqsBySSsimilarity.py',
           ' -a _RA_seqdat0 seq.dat file of the original secondary assignments',
           ' -b _RA_seqdat  seq.dat file containing predicted secondary assignments for the sequence ',
           ' -c __header header'
           ).parse(locals(),sys.argv)

sid_a=0.0 #sequence identity over secondary elements
sid_b=0.0 #sequence identity over all sequence
ss0=translateSec(seqdat0)
L=len(ss0)
nss=0 #number of residues in secondary elements
ss=translateSec(seqdat)  #;Bye(ss0+'\n'+ss)
weights=read_weights(seqdat) #confidence in the secondary structure predicted assignments

score=0.0

pat=re.compile('[^C]+')
it=pat.finditer(ss0)
for match in it:
    (b,e)=match.span()
    nss+=e-b
    c=1+(b+e)/2
    type=ss0[b] #'H' or 'E'
    h=halfmaxspan[type]
    for i in range(b,c):
        if ss[i]!=type:
            score+=weights[i]*(1+float(i-b)/h)
        else: sid_a+=1
    for i in range(c,e):
        if ss[i]!=type:
            score+=weights[i]*(1+float(e-i)/h)
        else: sid_a+=1
sid_a=sid_a/nss

pat=re.compile('C+')
it=pat.finditer(ss0)
h=halfmaxspan['C']
for match in it:
    (b,e)=match.span()
    c=1+(b+e)/2
    for i in range(b,c):
        if ss[i]!=type:
            score+=weights[i]*(1+float(i-b)/h)
    for i in range(c,e):
        if ss[i]!=type:
            score+=weights[i]*(1+float(e-i)/h)

avconf=sum(weights)/L #average confidence in the secondary predictions
ssavconf=0
for i in range(L):
    if ss0[i]==ss[i]: sid_b+=1
    if ss[i]!='C': ssavconf+=weights[i]
ssavconf=ssavconf/nss
sid_b=sid_b/L

print '#header score, ss_seqid, seqid, confidence, ss_confidence'
print '%s %5.1lf %4.2lf %4.2lf %4.2lf %4.2lf'%(header,score,sid_a,sid_b,avconf,ssavconf)

    
sys.exit(0)

