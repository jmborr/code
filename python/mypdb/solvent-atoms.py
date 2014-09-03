#!/usr/bin/python
import sys,os
from utilities.small_utilities import junkName
from inputArgs.inputArgs import inpHand
"""Calculate the solvent atoms close to the protein surface
   ******************************************************
   * DEBUG NOTE: this script still needs to be finished *
   ******************************************************"""

massf=junkName()
cmd='pdb_to_mass '+pdbf+' '+massf
os.system(cmd)

resListf=junkName()
cmd="cat "+pdbf+" | awk '{print($5)}' > "+resListf
os.system(cmd)

closeAtomsf=junkName9)
cmd='cat '+resListf+' '+massf+' '+crdf+' | solvent-atoms.x -na '+str(NallAtoms)+\
     ' -np '+NprotAtoms+' > '+closeAtomsf
os.system(cmd)


if __name__=='__main__':
    ih=inpHand('Usage: solvent-atoms.py req-args [opt-args]',
               ' -a _AR_pdbf reference PDB file of the solvated system',
               ' -b _AR_crdf AMBER coordinate file to calculate its close solvent molcules',
               ' -c _AR_closeAtomsf output file to store close solvent molecules'
               )
    
    ih.parse(locals(),sys.argv)
    
    sys.exit(0)        
