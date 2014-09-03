#!/usr/bin/python
"""
Shashi noticed that there are repeated sequences (identical sequences)
under different targets names. This code identifies if the sequences of two fasta files
are identical
"""
import sys,os
from utilities.small_utilities import chomp

fastas=chomp(open('/gpfs1/scratch/jose/newShortProteins/fastas.list','r').readlines())
myroot='/gpfs1/scratch/jose/newShortProteins/tasser'
out=open('/gpfs1/scratch/jose/newShortProteins/filteredByIdenticalFasta.list','w')

#import the sequences of all fasta files into a list
seqs=[]
for i in range(0,len(fastas)):
    f1=myroot+'/'+fastas[i][1]+'/'+fastas[i]+'/'+fastas[i]+'.fasta'
    inp=open(f1,'r')
    inp.readline() #forget about the header of the fasta file
    lines=chomp(inp.readlines())
    inp.close()
    seq=''.join(lines) #;print seq; sys.exit(1)
    seqs.append(seq)

#compare all the sequences
for i in range(0,len(seqs)-1):
    seq1=seqs[i]
    f1_repeated=0 #initialize flag   
    #print seq1
    for j in range(i+1,len(fastas)):
        if seq1==seqs[j]: #signal two identical sequences
            f1_repeated=1
            #print seq1+'\n'+seqs[j]
            #print fastas[i]+'\n'+fastas[j]
            #sys.exit(1)
            break
    if f1_repeated==0: out.write(fastas[i]+'\n') #not repeated sequence

sys.exit(1)
