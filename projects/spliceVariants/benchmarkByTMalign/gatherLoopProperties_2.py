#!/usr/bin/python

"""same as gatherLoopProperties.py but also include info from disorder predictors"""

import os,sys,re
from utilities.small_utilities import chomp

############################################################
#    GLOBAL VARIABLES
############################################################
minLoopL=4
minProximity=5
pat={}

root='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'

############################################################
#    DEFINITION OF FUNCTIONS
############################################################

def writeComments():
    return '''#C1 header
#C2 start of loop
#C3 end of loop
#C4 loop length
#C5 predicted exposed surface from 0.00(totally exposed) to 0.85 (totally buried)
#C6 predicted proportion of helix content
#C7 predicted proportion of extended content
#C8 predicted proportion of coil content
#C9 average rmsd of the loop
#C10 extended loop length
#C11 average extended rmsd (rmsd of the extened loop)
#C12 cluster density for the loop
#C13 cluster density for the extended loop
#C14 rmsd to native of the loop
#C15 rmsd to native of the extended loop
#       C1          C2  C3  C4  C5    C6    C7     C8   C9   C10   C11   C12     C13      C14   C15
'''

pat['isloopinfo']=re.compile('\s*\d+\s+\d+\s+\d+')
def isloopinfo(line):
    if pat['isloopinfo'].match(line): return line
    return ''

pat['isloopinfoII']=re.compile('\s*\d+\.\d+\s+\d+\.\d+')
def isloopinfoII(line):
    if pat['isloopinfoII'].match(line): return line
    return ''

def addLocalLoopDensity(valid,ro):
    Rc_in=float(valid.split()[7]) #assume <rmsd>=Rc_in in position 7
    if Rc_in<0.5: Rc_in=0.5
    return '%7.5lf'%(ro/Rc_in)

def addGlobalLoopDensity(valid,ro):
    Rc_in=float(valid.split()[9]) #assume <rmsd>=Rc_in in position 7
    if Rc_in<0.5: Rc_in=0.5
    return '%7.5lf'%(ro/Rc_in)

def looplength(line):
    return int(line.split()[2]) #assume loop length in position

def termProx(L,line):
    b=int(line.split()[0])   #distance to N-terminal
    e=L-int(line.split()[1]) #distance to C-terminal
    return min(b,e)

############################################################
#    MAIN FUNCTION
############################################################

pat['B']=re.compile('B-+\n 1\s+(\d+)')
pat['N']=re.compile('Number of structure in use=\s+(\d+)')
loopinfos=[]
nheader=0
#read directory of each successful TASSER simulation
listf=root+'/success.list'
#listf=root+'/toy.list'
for line0 in chomp(open(listf,'r').readlines()):
    tmp=[]
    header=os.path.basename(line0)
    nheader+=1
    #print header,nheader

    rmsinpf=root+'/out/'+header+'/rmsinp'
    L=int(open(rmsinpf,'r').readline().split()[1]) #sequence length from rmsinp

    rstdat=root+'/out/'+line0+'/rst.dat'
    all=''.join(open(rstdat,'r').readlines())
    N=int(pat['B'].search(all).group(1)) #cluster size
    N0=int(pat['N'].search(all).group(1)) #number of structures done in spicker clustering
    ro=(1.0*N)/N0 #number density
    
    loopdatf=root+'/out/'+line0+'/loop.dat'
    #read loop info independent of native state
    for line in chomp(open(loopdatf,'r').readlines()):
        valid=isloopinfo(line)
        if valid:
            valid+=' '+addLocalLoopDensity(valid,ro)+' '+addGlobalLoopDensity(valid,ro)            
            tmp.append(valid)

    #add loop comparison to native state
    n=0
    for line in chomp(open(loopdatf,'r').readlines()):
        valid=isloopinfoII(line)
        if valid:
            tmp[n]=tmp[n]+'  '+valid
            n+=1
    #filter by loop length
    tmp2=[]
    for line in tmp:
        if looplength(line)>minLoopL:tmp2.append(line)
    tmp=tmp2
    #filter by loop proximity to either termini
    tmp2=[]
    for line in tmp:
        #print line0,line,termProx(L,line),L
        if termProx(L,line)>minProximity: tmp2.append(line)
    tmp=tmp2
    #add tmp to loopinfos
    #print '\n'.join(tmp)
    for i in range(len(tmp)): tmp[i]=line0+' '+tmp[i]
    loopinfos=loopinfos+tmp

print writeComments()+'\n'.join(loopinfos)

sys.exit(0)
