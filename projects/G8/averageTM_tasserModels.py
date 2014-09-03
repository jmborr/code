#!/usr/bin/python
'''
outputs all pairwise TMaligns among top models across all G8 sequences
'''
import sys,os,re
from utilities.small_utilities import chomp,junkName,Bye

from spicker.spickerYangResultsManager import spickOut,TMscore
from TMalign.TMalignManageResults import TMalignOut

root='/gpfs1/scratch/jose/G8'

headers=chomp(open(root+'/input/list','r').readlines()) ; N=len(headers)

#save top models in memory
junk1=junkName()
models={}
currd=os.getcwd() ; os.system('mkdir '+junk1) ; os.chdir(junk1)
for i in range(0,N):
    header=headers[i]
    os.system('/bin/rm combo*pdb closc*pdb rst.dat 2>/dev/null')
    os.system('tar jxf '+root+'/out/'+header+'/'+header+'_out.tasser.tbz2 -C .')
    p=spickOut(dir='.')
    if p.readError:
        models[header]=[]
    else:
        models[header]=open(p.combo[p.densest],'r').readlines()
os.chdir(currd) ; os.system('/bin/rm -rf '+junk1)

#calculate all pairs of TMalign
junk1=junkName()
junk2=junkName()
junk3=junkName()
avtm=0.0  ;  n=0
tms={} ; nn={}
for i in range(0,N):
    tms[headers[i]]=0.0  ;  nn[headers[i]]=0 #initialization
for i in range(0,N-1):
    hi=headers[i] #; print hi
    if models[hi]==[]: continue
    open(junk1,'w').writelines(models[hi])
    l1=len(models[header])
    for j in range(i+1,N):
        hj=headers[j] #; print '  '+hj
        if models[hj]==[]: continue
        open(junk2,'w').writelines(models[hj])
        l2=len(models[hj])
        if l1<l2: cmd='TMalign '+junk1+' '+junk2+' > '+junk3
        else:  cmd='TMalign '+junk2+' '+junk1+' > '+junk3
        os.system(cmd)
        tm=TMalignOut(junk3).tm ; print tm
        if tm<1.0:
            tms[hi]+=tm ; tms[hj]+=tm
            nn[hi]+=1   ; nn[hj]+=1
        avtm+=tm  ;  n+=1
print avtm/n

#normalize cummulative TMscores for each target
for i in range(0,N):
    header=headers[i]
    if nn[header]>0: tms[header]=tms[header]/nn[header]
#find target with highest cummulative TMscore
bestHeader='' ; besttm=0.0    
for i in range(0,N):
    header=headers[i]    
    if tms[header]>besttm:
        besttm=tms[header] ; bestHeader=header
print besttm,bestHeader

os.system('/bin/rm '+junk1+' '+junk2+' '+junk3)
sys.exit(0)
