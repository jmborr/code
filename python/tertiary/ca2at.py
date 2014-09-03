#!/usr/bin/python
"""given a CA pdb, build a minimized all-atom model
(1) pulchra only on CA atoms
(2) dukka's
(3) pulchra to reconstruct to all-atom
(4) amber protocol minimization"""

import os,sys
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName
from utilities.codedir import codedir
from jobs.job import pastry

#bring input
pulchra=os.path.join(codedir,'c/pulchra/p2')
dukkax='/gpfs1/u/jose/bin/dukka.x'
amberminimize=os.path.join(codedir,'bash/tertiary/amberminimize.sh')
pastry('/bin/cp '+pulchra+' '+dukkax+' '+amberminimize+' .')

def ca2at(file='',file_pointer='',list=''):
    #extract CA's only from input
    if file:
        list=chomp(open(file,'r').readlines())
    if file_pointer:
        list=chomp(file_pointer.readlines())
    caonly=[]
    for line in list:
        if " CA " in line: caonly.append(line)
    #create temporary working directory and work there
    currd=os.getcwd()
    wd=os.path.join(currd,junkName())
    os.system('/bin/mkdir '+wd)
    os.chdir(wd)
    #dump CA's onto file
    capdbf=junkName()
    open(capdbf,'w').write( '\n'.join(caonly) )
    #run pulchra only CA atoms
    if os.system('./p2 -b '+capdbf): return []
    #run dukka on refined CA's
    if os.system('./dukka.x '+capdbf+'.rebuilt'+' '+capdbf): return []
    #reconstruct all-atom
    initpdb=junkName()
    if os.system('./p2 best.pdb && /bin/mv best.pdb.rebuilt '+initpdb): return []
    #amber minimize
    if os.system('./amberminimize.sh '+initpdb+' .'): return []
    #read file and return
    minl=chomp(open(initpdb+'.min','r').readlines())
    os.chdir(currd)
    os.system('/bin/rm -rf '+wd) #clean-up
    return minl        
    
if __name__=='__main__':
    inpHand('Usage: ca2at.py [options]',
            ' -a _RA_pdbf PDB file with CA representation',
            ' -b _RA_outf output PDB file with all-atom representation',
            ).parse(locals(),sys.argv)

    minl=ca2at(file=pdbf)
    if not minl: sys.exit(1) #error
    open(outf,'w').write( '\n'.join(minl)+'\n' )
    sys.exit(0)
