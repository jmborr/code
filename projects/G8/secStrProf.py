#!/usr/bin/python
'''
Construct secondary structure profile from the global sequence
alignment of the 70 homologs to the reference sequence, then
substitute the amino acid sequence by the secondary structure states

'''

import sys,os,re
from utilities.small_utilities import Bye
from seq.alignYangManageResults import gappedSeq,alignOut
from seq.letters import valid_A

root='/gpfs1/scratch/jose/G8'
numberToSymbol={'2':'H', '4':'E', '1':'C'}

#return sequence of secondary assignments by reading seq.dat file
def getSS(header):
    ss=''
    cmd='tar jxf '+root+'/out/'+header+'/'+header+ '_tasserinp.tbz2 seq.dat' #extract seq.dat
    if os.system(cmd): return ''
    lines=open('seq.dat','r').readlines()
    for line in lines: ss+=numberToSymbol[ line.split()[2] ]
    return ss
    
s1=root+'/input/fastas/g8000' #fasta file
ss1=getSS('g8000')            #sequence of secondary assignments
L=125 #sequence length of g8000

#initialize profile
profile={}
for i in range(0,L):
    profile[i]={'H':0.0,'E':0.0,'C':0.0,}
    profile[i][ss1[i]]+=1

N=70 #number of homologs to g8000
for i in range(1,N+1):
    header='g80'+'%02d'%(i)  #;  print header
    s2=root+'/input/fastas/'+header
    ss2=getSS(header)
    if not ss2:continue
    algL=os.popen('$HOME/bin/align '+s1+' '+s2+' 0').readlines() #store alignment in a list
    alg=alignOut(algL).alg #alignment object
    gs1=alg.gs[1]  ;  gs2=alg.gs[2]  #gapped sequences
    for i in range(0,L):
        j=gs1.si2gsi[i] #index from ungapped sequence to gapped sequence for gs1
        k=gs2.gsi2si[j] #index from gapped sequence to ungapped sequence for gs2
        #if position j in gs2 is occupied by an amino acid then k>=0
        if k>=0: profile[i][ ss2[k] ]+=1 #add corresponding secondary structure assignment

#normalize the profiles
for i in range(0,L):
    total=0.0
    for ss in profile[i].keys(): total+=profile[i][ss]
    for ss in profile[i].keys(): profile[i][ss]=profile[i][ss]/total
    print '%3d %5.2lf %5.2lf %5.2lf'%(i+1,profile[i]['C'],profile[i]['E'],profile[i]['H'])

os.system('/bin/rm seq.dat') #some clean-up
sys.exit(0)
    
