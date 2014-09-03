#!/usr/bin/python

import os,sys,re
from utilities.small_utilities import chomp,junkName
from seq.alignYangManageResults import alignOut

#global variables
rootc='/gpfs1/active/jose/code/projects/spliceVariants/benchmark/preparing_database'
root='/gpfs1/scratch/jose/spliceVariants/benchmark/preparing_database'
ph=re.compile('>\s*(.*)\s*\n') #extract the header out of the fasta header line
paa=re.compile('[G,A,V,L,I,S,T,C,M,P,D,N,E,Q,K,R,H,F,Y,W,B,Z,X]') #check for amino acid letter

#return sequence and header of current fasta entry where pin is pointing to
def importFastaEntry(pin):
    seq='' ;  header=''

    #check end-of-file nad beginning of fasta entry
    while not header:
        header=pin.readline()
        if not header: return '','' #check end-of-file
        header=ph.match(header).group(1) #check beginning of fasta entry

    #read fasta sequence
    line=chomp( pin.readline() )
    while paa.match(line): #correspond to an amino acid letter
        seq+=line
        n=pin.tell()
        line=chomp( pin.readline() )

    if line: pin.seek(n) #go back one line, since we over-read in the previous loop
    return header,seq

#load fasta files into memory
pin=open(root+'/pdb.40.300','r')
while 1:
    header,seq=importFastaEntry(pin) #; print header,seq
    if not header: break #reached end-of-file
    print '%s %3d'%(header,len(seq))
pin.close()
