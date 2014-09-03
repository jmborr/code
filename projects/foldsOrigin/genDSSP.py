#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp
from Bio.PDB.DSSP import PDBParser,DSSP

def dsspFastaLike(header='',simpl=0):
    f='>'+header+'\n'
    g=''
    pdbf='/library/pdb/pdb'+header[0:4]+'.ent'
    cid=header[4].upper()
    if cid=='_': cid=' '
    p=PDBParser()
    s=p.get_structure(header[0:4],pdbf)  #structure
    m=s[0] #model
    d=DSSP(m,pdbf) #DSSP object
    for i in s[0][cid]:
        x=i.xtra
        if x:
            if x.has_key('SS_DSSP'):
                x=x['SS_DSSP']
                if simpl and x!='H' and x!='E' : x='-'
                g=g+x
    while g:
        f=f+g[0:60]+'\n'
        g=g[60:]
    return f
    
ih=inpHand('Usage: genDSSP.py [options]\nExample:',
           '  -a _AR_list file list of five letter codes (pdb id plus chain id)',
           '  -b __outf file (def:STDOUT)',
           '  -c __simpl simplified alphabet (H,E,-) (default:\'no\')'
	)
ih.parse(locals(),sys.argv)

if not outf: out=sys.stdout
else       : out=open(outf,'w')
if not simpl:  simpl=0
else:  simpl=1

for header in chomp(open(list,'r').readlines()):
    out.write( dsspFastaLike(header=header,simpl=simpl)+'\n' )
    
out.close()
