#! /usr/bin/env python

from pymol import cmd

def dat2bf(mol,file,col=1):
    """change bf of selected object with the contents of a data file
    mol: selected object
    file: file containing b-factors
    col: column in 'file' corresponding to the b-factors
    """
    pin=open(file)
    m=cmd.get_model("%s"%mol)
    for at in m.atom:
         l=pin.readline()
         while l and l[0] in ('#',): l=pin.readline() #skip comment lines
         b=float( l.split()[int(col)-1] )
         c=at.chain; r=at.resi; n=at.name
         cmd.alter( '/%s//%s/%s/%s'%(mol,c,r,n) , 'b=%f'%b )
    cmd.spectrum("b") #refresh the coloring
    
cmd.extend("dat2bf",dat2bf)

def dat2bfByres(mol,file,col=1,factor=1.0,comment_line_header='#'):
    """change bf of selected object with the contents of a data file
    mol: selected object
    file: file containing b-factors by residue
    col: column in 'file' corresponding to the b-factors
    """
    pin=open(file)
    m=cmd.get_model("%s"%mol)
    rp=None; b=None
    for at in m.atom:
         c=at.chain; r=at.resi; n=at.name
         if not rp or rp != r:
             rp=r
             l=pin.readline()  #read another residue
             while l and l[0]==comment_line_header: l=pin.readline() #skip comments
             b=float( l.split()[int(col)-1] ) * factor
         cmd.alter( '/%s//%s/%s/%s'%(mol,c,r,n) , 'b=%f'%b )
    cmd.spectrum("b") #refresh the coloring
    
cmd.extend("dat2bfByres",dat2bfByres)
