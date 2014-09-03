#!/usr/bin/python

from seq.fastaManager import importFastaEntry

def calc_wall_time(L=0,fastaf=''):

    """estimate running time for a TASSER simulation"""
    
    if fastaf: #we pass a fasta file
        header,seq=importFastaEntry( open(fastaf,'r') )
        L=len(seq)

    minwt=0.99
    if L<100: wt=minwt
    else: wt=minwt+(L-100.0)/200*5
    return '%5.2lf'%(wt,)


def number_HB(pdb):

    """  **********    U N F I N I S H E D !!!!  ***********
    output number of hydrogen bonds according to TASSER definition

    argument pdb can be a file, a string, a list, or a file pointer"""

    from  inputArgs.inputArgs import resolveArgument

    lines=[]
    resolveArgument(pdb,outfmt='_LIST_',out=lines)

    ca=[]
    for line in lines:
        if 'ATOM' in line and ' CA ' in line:
            ca.append( float(line[30:38]),float(line[38:46]),float(line[46:54]) )

    N=len(ca)
    for i in range(N-1):
        cai=ca[i]
        for j in range(N):
            caj=ca[j]


    
    


