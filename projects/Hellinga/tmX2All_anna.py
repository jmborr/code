#!/usr/bin/python
"""

"""
import os,sys,re
from utilities.small_utilities import Bye,chomp,junkName
from utilities.codedir import scratchdir

tmco=0.4
scd=scratchdir+'/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/tmAll2All/anna_atasser_minimized'
listf=scd+'/structures/list'

#init buffer to hold hits
buffer={}
lengths={}
for header in open(listf,'r').readlines():
    header=header.strip()
    buffer[header]=[]
    lengths[header]=0

#local directory to avoid gpfs traffic
localoutd='/tmp/jose/tmX2all'
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
    header1,header2,tm=items #;print tm
    tmf=float(tm)
#    tmf=float(tm)*(150.0/140) #length trick (we align a 140-residue sequence to a 150-res one)
    if tmf>tmco:
        buffer[header1].append(header2+' '+tm)
        lengths[header1]+=1
        buffer[header2].append(header1+' '+tm)
        lengths[header2]+=1
        #flush buffer component
        for headeri in (header1,header2):
            listx=buffer[headeri]
            if len(listx) > 200:
                open(localoutd+'/'+headeri+'.dat','a').write( '\n'.join(listx)+'\n')
                buffer[headeri]=[]
    line=pin.readline()
    
for header in buffer.keys():
    f=localoutd+'/'+header+'.dat'
    #flush buffer
    listx=buffer[header]
    open(f,'a').write( '\n'.join(listx)+'\n')
    buffer[header]=[]
    #prepend number of hits
    n=lengths[header]
    buff=''.join( open(f,'r').readlines() )
    open(f,'w').write(`n`+'\n'+buff)

os.system('/bin/mv -f '+localoutd+' '+scd) #move local stuff

sys.exit(0)
