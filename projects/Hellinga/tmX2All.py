#!/usr/bin/python
"""

"""
import os,sys,re
from jobs.job import pastry,genJobList
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName

tmco=0.4
scd='/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/tmAll2All'


str2link={}
for line in open(scd+'/structures/structure2link.txt','r').readlines():
    str,link=line.split()
    str2link[str]=link

buffer={}
lengths={}

#for link in str2link.values():
#    buffer[link]=[]
#    lengths[link]=0

#local directory to avoid gpfs traffic
localoutd='/tmp/jose/tmX2all_140'
os.system('/bin/mkdir -p '+localoutd)

pin=open(scd+'/tmAll2All.dat','r')
line=pin.readline()
n=0
while line:
    n+=1
    if len(line)<2:
        line=pin.readline()
        continue
    items=line.split()    
    if len(items)!=3:
        line=pin.readline()
        continue
    str1,str2,tm=items #;print tm
    key1=str1 ; key2=str2
#    key1=str2link[str1] ; key2=str2link[str2]
    for key in (key1,key2):
        if key not in lengths.keys(): lengths[key]=0
        if key not in buffer.keys(): buffer[key]=[]
    tmf=float(tm) #*(150.0/140)
    if tmf>tmco and str1[0]!='d' and str2[0]!='d':
#       buffer[key1].append(key2+' '+tm )
        lengths[key1]+=1
#       buffer[key2].append(key1+' '+tm )
        lengths[key2]+=1
#        for stri in (str1,str2):
#            link=str2link[stri]
#            listx=buffer[link]
#            if len(listx) > 200: #flush buffer component
#                open(localoutd+'/'+link+'.dat','a').write( '\n'.join(listx)+'\n')
#                buffer[link]=[]
    line=pin.readline()
    

#for link in buffer.keys():
for link in lengths.keys():
    n=lengths[link]
    open(localoutd+'/'+link+'.dat','w').write(`n`+'\n')
 #   buff=open(localoutd+'/'+link+'.dat','r').readlines()
 #   open(localoutd+'/'+link+'.dat','w').write(`n`+'\n'+buff)


#os.system('/bin/mv -f '+localoutd+' '+scd) #move local stuff

sys.exit(0)


