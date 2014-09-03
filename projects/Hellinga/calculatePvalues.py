#!/usr/bin/python
'''
read all.instructures.dat file that is the merge of all a0112_xxxx/abinitio/comboXX.instructures.dat.

'''
import sys,os
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb

ih=inpHand('Usage: calculatePvalues.py [options]',
           ' -a _AR_allf all.instructures.dat file',
           ' -b __targetstr target structure (def:a0112)',
           ' -c __outf output file (def: STDOUT)'
           )
ih.parse(locals(),sys.argv)

if not targetstr:targetstr='a0112'
if not outf: outp=sys.stdout
else: outp=open(outf,'w')

#create the cummulative distribution of TM-scores for each structure
#against all combo models. We do NOT read allf into memory because it
#is a huge file. Rather, we read it line by line.
tmdist={}
structures=[]
allp=open(allf,'r')
line=allp.readline() 
while line:
    if line[0]!='#': #line is not a comment
        items=line.split()
        structure=items[2]
        tm=int(float(items[5])*100) #; print 'tm=',tm #two decimal places
        if structure not in structures:
            structures.append(structure)
            tmdist[structure]=[0]*100 #initialize the histogram
        tmdist[structure][tm]+=1            
    line=allp.readline()
allp.close()
for structure in structures:
    i=98
    while i>=0:
        tmdist[structure][i]+=tmdist[structure][i+1] #cummulative distribution
        i-=1
    area=tmdist[structure][0]
    for j in range(100): tmdist[structure][j]=float(tmdist[structure][j])/area
#for each tm-cutoff, calculate number of structures with a fraction of
#models having tm>tmco bigger than the fraction of models having
#tm>tmco to target structure
nstructures=len(structures)
buf=''
for tmco in range(100):
    over=-1 #begin with -1 because we are counting targetstr in the for loop
    f0=tmdist[targetstr][tmco]
    for structure in structures:
        if tmdist[structure][tmco]>=f0: over+=1
    pval=float(over)/nstructures
    buf+='%4.2lf %6.4lf %8.6lf\n'%(float(tmco)/100,pval,f0)

outp.write(buf)
