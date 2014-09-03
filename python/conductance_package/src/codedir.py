#!/usr/bin/python
import os,sys
from small_utilities import chomp

#variables derived from .bashrc
homedir=os.environ['HOME']
codedir=os.environ['CODED']
projectsdir=os.environ['PROJD']
scratchdir=os.environ['SCRTD']


pdbnonhom='/gpfs1/archive/adrian/dat/pdb_may06.apr06'
CAdir='/gpfs1/archive/adrian/dat/pdb_may06.apr06.aux/CA'
libseq=pdbnonhom+'/seq'
pdbdir='/library/pdb'
yanglib='/library/yzhang'
yangpdb='/library/yzhang/PDB'
nrfiltdb='/library/yzhang/nr/nr.filter'
pdbhomolrestricted='/gpfs1/UB/net/dell/01/users/skolnick/pdbhomoljul05'

scopd=scratchdir+'/db/scop/scop.1.71'
archdir='/gpfs1/archive/jose'
psipredd=codedir+'/python/secondary/psipred'
ncbidir=codedir+'/bin/blast'
nodeLib='/local/library/jose'
codeNodeLib=nodeLib+'/code'
nodecodedir=nodeLib+'/code'

prefixd={'homedir':homedir,'libseq':libseq,'pdbdir':pdbdir,
         'pdbnonhom':pdbnonhom,'pdbhomolrestricted':pdbhomolrestricted,
         'codedir':codedir,'projectsdir':projectsdir,'scratchdir':scratchdir}

if __name__=='__main__':
    print codedir,projectsdir
