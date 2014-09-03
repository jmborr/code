#!/usr/bin/python

""" calculate the percentage of residues assigned to be in coil state
over the constitutive region for a reference sequence and a variant.
Secondary structures of reference and variant derived from seq.dat predictions.
"""

import os,sys
from inputArgs.inputArgs import inpHand
from seq.alignYangManageResults import alignOut
from secondary.seqdat2ss import translateSec
from jobs.job import Bye

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: coilsPercentageConstitutive.py [options]',
               '  -a _R_var variant ID (example: "Q8R09202")',
               )
    ih.parse(locals(),sys.argv)

    #init reference DI
    root=var[0:6]
    ref=root+'00'

    #find sequence alignment between reference and variant
    scd='/gpfs1/scratch/jose/spliceVariants'
    cmd='grep '+var+' '+os.path.join(scd,'preparing_database/var0_var2.alg.idx | cut -d " " -f 1')
    pos=int( os.popen(cmd).readline().strip() ) #where to look in var0_var2.alg
    pin=open(os.path.join(scd,'preparing_database/var0_var2.alg'),'r')
    pin.seek(pos)
    #list of pair of ungapped-indexes for identical residues in the alignment
    pairs=alignOut(pin).alg.listIndexPairsIdentical()
    #find coil content over the identical pairs
    cmd='tar xf '+os.path.join(scd,'out',root,'00',ref+'.in.tasser.tar seq.dat')
    os.system(cmd)
    refss=translateSec('seq.dat') #one line of secondary structure symbols
    os.system('tar xf '+os.path.join(scd,'out',root,var[-2:],var+'.in.tasser.tar seq.dat'))
    varss=translateSec('seq.dat')
    id=0
    for pair in pairs:
        if refss[ pair[0] ]=='C': id+=1
        if varss[ pair[1] ]=='C': id+=1
    buf='%4.2f '%(0.5*float(id)/len(pairs),)
    #find coil content over the unaligned region (here assumed the variant region)
    idref=[]
    idvar=[]
    for pair in pairs:
        idref.append(pair[0])
        idvar.append(pair[1])          
    id=0 ; n=0
    for i in range(len(refss)):
        if i not in idref:
            n+=1
            if refss[i]=='C': id+=1
    for i in range(len(varss)):
        if i not in idvar:
            n+=1
            if refss[i]=='C': id+=1            
    buf+='%4.2f'%(float(id)/n)
    print buf
    os.system('/bin/rm seq.dat')
    sys.exit(0)
         
        
