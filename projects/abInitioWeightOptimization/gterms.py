#!/usr/bin/python


'''******************** Pseudocode *****************
1 Input energy terms file
Assumed File format:
#amt       EHB1       EHB1a       EHB1b        EHB1c       EHB2        EHB3         EHB4       EHB5a        EHB5b    ESHORT2     ESHORT5      ESHORT5a    ESHORT5b    ESHORT5c   ESHORT6      ESHORT7    ESHORT8      eprofo      E_cord      E_cnum
# 1         2           3           4           5           6           7           8           9          10          11         12           13          14          15          16          17          18          19          20          21
0.000      9.8679      4.6943     36.0000      1.4000     66.5100    -99.6463    -20.6437      0.0000     -0.7891     -4.8842   -113.5359     15.0000    107.4666     -7.3506      9.9800     88.7000     47.5350      0.0000      5.4282     46.0000
2 calculate <amt>, <amt^2>, <amt*Ee>, <Ee>, <Ee*Ee">, Ene where <> denotes average over all decoys, and Ee the different energy terms. Ene are the different energy terms for the native state.
'''


import sys,os,re
from utilities.small_utilities import chomp
from random import randint
from inputArgs.inputArgs import inpHand

ih=inpHand(
    'Usage: tra2chain.py ',
    '  -a _AR_etermsf energy terms file',
    '  -b _A_outf output file (def: STDOUT)'
    )
ih.parse(locals(),sys.argv)
if not outf: pout=sys.stdout
else: pout=open(outf,'w')

M=20 #number of energy terms
M1=M-1

pin=open(etermsf,'r')
#read the first two lines, which are comments
line=pin.readline()
pout.write(line) #output the list of energy terms
line=pin.readline()
#The following line are the energy terms for the native state, and 1-TM==0 for native
enterms=pin.readline().strip().split()[1:]
for i in range(M): enterms[i]=float(enterms[i])
#load the rest of lines into memory
amt=[]
eterms=[] #eterms[i][j]: energy term "j" for decoy "i"
line=pin.readline()
while(line):
    list=line.split()
    for i in range(len(list)): list[i]=float(list[i]) #convert to float
    amt.append(list[0])     #1-TMscore
    eterms.append(list[1:]) #energy terms
    line=pin.readline()
#initialize averages
N=len(amt)     #number of decoys
avmt=0.0       #<mt>
avmt2=0.0      #<mt^2>
avmtE=[0.0]*M  #<mt*Ee>  e=1,..,M
avE=[0.0]*M    #<Ee>     e=1,..,M
corr=[]
for i in range(M): corr.append([0.0]*M) #<Ee*Ee"> e=1,..,M , e"=1,..,M

#calculate averages
for i in range(N): #cycle over all decoys
    avmt+=amt[i]
    avmt2+=amt[i]*amt[i]
    etermsi=eterms[i] #a handy reference
    for j in range(M): #cycle over all energy terms
        avmtE[j]+=amt[i]*etermsi[j]
        avE[j]+=etermsi[j]
        for k in range(M): corr[j][k]+=etermsi[j]*etermsi[k] #redundancy, since EeEe"==Ee"Ee

#normalize
avmt=avmt/N
avmt2=avmt2/N
for j in range(M):
    avmtE[j]=avmtE[j]/N
    avE[j]=avE[j]/N
    for k in range(M):
        corr[j][k]=corr[j][k]/N
        
#output
line=' %5.3lf %8.6lf'%(avmt,avmt2)
pout.write('#<mt> <mt^2>\n'+line+'\n')

line=''
for j in range(M): line+=' %14.7lf'%(avmtE[j])
pout.write('#<mt*Ee> (1..M)\n'+line+'\n')

line=''
for j in range(M): line+=' %11.4lf'%(avE[j])
pout.write('#<Ee> (1..M)\n'+line+'\n')

line=''
for j in range(M):
    for k in range(M):
        line+=' %17.6lf'%(corr[j][k])
pout.write('#<Ee*Ee"> (1,1)..(1,M),..(M,M)\n'+line+'\n')

line=''
for j in range(M): line+=' %11.4lf'%(enterms[j])
pout.write('#Ene (1..M)\n'+line+'\n')



pin.close()
sys.exit(0)
