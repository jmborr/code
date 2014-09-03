#!/usr/bin/python
"""

assign secondary structure (helix H, strand E, coil C, turn T) from CA atom
coordinates. The algorithm is the one employed in the assignment of
secondary structure in the TMalign program

"""
import os,sys
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand
from math import sqrt

#return square of the distance, using global variable xyz
def dis(xyz,i,j):
    d2=0.0
    for k in range(0,3):
        delta=xyz[k][i]-xyz[k][j]
        d2+=delta*delta
    return sqrt(d2)
    
def evalSec(dis13,dis14,dis15,dis24,dis25,dis35):
    ss='C'
    #check for helix
    delta=2.1
    if abs(dis15-6.37)<delta and abs(dis14-5.18)<delta and abs(dis25-5.18)<delta:
        if abs(dis13-5.45)<delta and abs(dis24-5.45)<delta and abs(dis35-5.45)<delta:
            return 'H' 
    #check for strand
    delta=1.42
    if abs(dis15-13) < delta and abs(dis14-10.4)<delta and abs(dis25-10.4)<delta:
        if abs(dis13-6.1)<delta and abs(dis24-6.1)<delta and abs(dis35-6.1)<delta:
            return 'E'
    #check for turn
    if dis15<8.0: return 'T'
    #leftover is coil
    return 'C'

#smooth the secondary structure assignment
def smooth(ss):
    L=len(ss)
    tt=ss[0:2]
    #--x-- => -----
    #if position "i" is helix or strand, then set as coil if
    #surrounding four positions are not of the same secondary
    #structure type
    for i in range(2,L-2):
        s=ss[i]
        if s=='H' or s=='E':
            if ss[i-2]!=s and ss[i-1]!=s and ss[i+1]!=s and ss[i+2]!=s: s='C'
        tt+=s
    tt+=ss[-2:]
    #--xx-- => ------
    #if position "i+2" and "i+4" is helix or strand, then set as coil
    #if surrounding four positions are not of the same secondary
    #structure type
    uu=tt[0:2]
    i=2
    while i<L-2:
        s=tt[i] ; s2=tt[i+1]
        if s==s2 and (s=='H' or s=='E'):
            if tt[i-2]!=s and tt[i-1]!=s and tt[i+2]!=s and tt[i+3]!=s:
                s='C' ; uu+=s ; i+=1
        uu+=s ; i+=1
    uu+=tt[-2:]
    #x-x => xxx
    tt=uu[0]
    for i in range(1,L-1):
        t=uu[i]
        for s in ('H','E'):
            if uu[i-1]==s and uu[i+1]==s:
                t=s
                break
        tt+=t
    tt+=uu[-1]
    return tt

def outputSec0(x,y,z):
    ss='CC'
    xyz=[x,y,z]
    L=len(x) #sequence length
    for i in range(2,L-2):
        dis13=dis(xyz,i-2,i  )
        dis14=dis(xyz,i-2,i+1)
        dis15=dis(xyz,i-2,i+2)
        dis24=dis(xyz,i-1,i+1)
        dis25=dis(xyz,i-1,i+2)
        dis35=dis(xyz,i  ,i+2)
        ss+=evalSec(dis13,dis14,dis15,dis24,dis25,dis35)
    return ss+'CC'
    
def outputSec(pdbfile,smooth=smooth):
    """method to evaluate secondary structure of all residues except for
    the first and last two residues

    pdbfile: PDB file or list of lines containing the pdb file
    """
    cas=[]
    if isinstance(pdbfile,list):
        for x in pdbfile:
            if ' CA ' in x: cas.append(x)
    else:
        cas=chomp(os.popen('grep -e "ATOM " '+pdbfile+'|grep -e " CA "').readlines())
    #retrieve coordinates
    x=[] ; y=[] ; z=[]
    for ca in cas:
        x.append( float(ca[30:38]) )
        y.append( float(ca[38:46]) )
        z.append( float(ca[46:54]) )
    return smooth( outputSec0(x,y,z) )

if __name__=='__main__':
    inpHand('Usage: pdb2ss.py [options]',
            ' -a _RA_pdbfile PDB file',
            ' -b _A_outfile output file (def:STDOUT)'         
            ).parse(locals(),sys.argv)

    #determine output
    if not outfile: pt=sys.stdout
    else: pt=open(outfile,'w')
    #write secondary structure
    pt.write(outputSec(pdbfile)+'\n')
    pt.close()
    sys.exit(0)
