#!/usr/bin/python
'''
outputs all pairwise TMaligns among best templates across all G8 sequences
'''
import sys,os,re
from utilities.small_utilities import chomp,junkName,Bye

from prospector.prospManager import prospOut
from TMalign.TMalignManageResults import TMalignOut

root='/gpfs1/scratch/jose/G8'

def save_best_template(header,junkf):
    templatef=root+'/out/'+header+'/'+header+'rap3orienrev4s.pdb'
    if not os.path.exists(templatef): return 0
    pin=open(templatef,'r')
    line=pin.readline()
    while(line[0:4]!='ATOM'): line=pin.readline()
    pout=open(junkf,'w')
    n=0
    while(line[0:4]=='ATOM'):
        pout.write(line)  ;  n+=1  ;  line=pin.readline()
        if not line:break #reached end of file
    pout.close()

headers=chomp(open(root+'/input/list','r').readlines())
junk1=junkName()
junk2=junkName()
junk3=junkName()
N=len(headers)
avtm=0.0  ;  n=0
for i in range(0,N-1):
    hi=headers[i] #; print hi
    l1=save_best_template(hi,junk1)
    if l1==0: continue
    for j in range(i+1,N):
        hj=headers[j] #; print '  '+hj
        l2=save_best_template(hj,junk2)
        if l2==0: continue
        if l1<l2: cmd='TMalign '+junk1+' '+junk2+' > '+junk3
        else:  cmd='TMalign '+junk2+' '+junk1+' > '+junk3
        os.system(cmd)
        tm=TMalignOut(junk3).tm ; print tm
        avtm+=tm  ;  n+=1
print avtm/n
os.system('/bin/rm '+junk1+' '+junk2+' '+junk3)
sys.exit(0)
