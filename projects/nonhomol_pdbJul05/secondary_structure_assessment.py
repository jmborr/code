#!/usr/bin/python

import os,sys
from utilities.small_utilities import Bye,chomp

all='/gpfs1/scratch/jose/nonhomol_pdbJul05/analysis/summaryMoreBelow200.all'
root='/gpfs1/scratch/jose/createInputsTasser/'
root2='/gpfs1/archive/adrian/dat/pdb_jul05/seq'

print '#(1)header    (2)TM of best model    (3)fraction of residues with identical secondary structure assignment. The methods are seqdat2ss.py and pdb2ss.py    (4)(5)(6) fraction of C,H,E for seqdat2ss.py   (7)(8)( fraction of C,H,E for pdb2ss.py'
print '#  1    2    3   4     5    6    7    8    9'

for line in open(all,'r').readlines():
    if line[0]=='#': continue #comment line
    items=line.split()
    header=items[0]
    tm=items[19]
    pdbf=root+'/'+header[1]+'/'+header+'/abinitio/CA'
    ss_from_pdb=chomp(os.popen('pdb2ss.py -a '+pdbf).readline()) #sec.str. from PDB
    seqdatf=root2+'/'+header+'.SEQ'
    ss_from_seqdat=chomp(os.popen('seqdat2ss.py -a '+seqdatf).readline()) #sec.str. from seq.dat

    naa=len(ss_from_pdb)  ;  id=0.
    fract={'C':0. , 'H':0. , 'E':0. }
    fract2={'C':0. , 'H':0. , 'E':0. }
    for i in range(naa):
        x=ss_from_seqdat[i] ; fract[x]+=1
        y=ss_from_pdb[i]
        if y!='H' and y!='E': y='C'  #pdb2ss contains more symbols than C, H, E
        fract2[y]+=1
        if x==y: id+=1
    for key in fract.keys():
        fract[key]=fract[key]/naa
        fract2[key]=fract2[key]/naa
    print header+' %s %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf'%(tm,id/naa,fract['C'],fract['H'],fract['E'],fract2['C'],fract2['H'],fract2['E'] )
        
sys.exit(0)
