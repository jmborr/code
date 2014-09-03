#!/usr/bin/python

import os,sys
from utilities.small_utilities import chomp,junkName,Bye
from TMalign.TMalignManageResults import TMalignOut
from seq.alignYangManageResults import gappedSeq

'''

For every (structural) alignment in alignmentsByTM.dat, find a PDB
structure in the pdb_aug06.apr06 database that will be template for
both aligned sequences.

For a structural alignment between aaaaa.pdb and bbbbb.pdb, instead of
looking all PDB structures, we gather all structures in aaaaa.hist and
bbbbb.hist in /gpfs1/active/skolnick/struhistaug06 and we will do
TMalign of these structures to aaaaa.pdb and bbbbb.pdb

For a structural alignment between aaaaa.pdb and bbbbb.pdb, let be L
the number of aligned residues. Then, for a particular template
ccccc.pdb in either aaaaa.hist or bbbbb.hist, we do
TMalign(aaaaa,ccccc) and TMalign(bbbbb,ccccc) with d0 parameter given
by L, ie., d0(L). Let Na_bc be the number of residues of aaaaa that
are aligned to both bbbbb and ccccc (similarly Nb_ac). We aim to find
the template ccccc such that Na_bc*Nb_ac is maximal. This will be the
common template that we assign to both aaaaa and bbbbb.

'''

SCroot='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
templd='/gpfs1/archive/adrian/dat/pdb_aug06.apr06' #template directory
pdbd=templd+'/CA'
seqd=templd+'/input'
histd='/gpfs1/active/skolnick/struhistaug06'
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'

junk=junkName()

print '#short  long  best_templ  max_both_shared'
align_in=open(SCroot+'/alignmentsByTM.dat','r')
for line in open(SCroot+'/alignmentsByTM.idx','r'):
    align_in.seek( int(line.split(':')[0]) ) #go to the position in the file
    alg=TMalignOut(align_in) #read one entry in the file containing the structural alignments
    if alg.readError: break  #read error or reached end of file
    nalg=alg.nAligned        #number of aligned residues
    short=alg.Chain[1][0:5]  #header of the short sequence
    long =alg.Chain[2][0:5]  #header of the long sequence
    alg.initIsAlignedArray() #initialize alg.isAlined dictionary
    arrayOfAligned={short:alg.isAligned[1],long:alg.isAligned[2]}
    #gather all templates in short.hist and long.hist into a single list
    templs=[] #list of templates
    for header in (short,long):
        histfile=histd+'/'+header+'.hist'
        hist_in=open(histfile,'r')
        for line in open(histfile,'r').readlines()[2:]: #first two lines are useless
            templ=line.split()[0]
            if templ not in templs and templ not in (short,long):
                templs.append(templ) #add new template
    #for every template in templs, calculate TMalign to long and
    #short, then calculate the residues of short and long that aligned
    #to the the template    
    max_both_shared=0.0
    best_templ=''
    for templ in templs:
        pdbtempl=templd+'/CA/'+templ+'.pdb'
        both_shared=1
        for header in (short,long):
            cmd='TMalign '+templd+'/CA/'+header+'.pdb '+pdbtempl+'  -L '+`nalg`+' > '+junk
            os.system(cmd)
            alg2=TMalignOut(junk)
            alg2.initIsAlignedArray()
            shared=0.0 #number of residues of header aligned in both alg and alg2
            for i in range( len(arrayOfAligned[header]) ):
                shared+=arrayOfAligned[header][i]*alg2.isAligned[1][i]
            both_shared=both_shared*shared
        both_shared=both_shared/(nalg*nalg)
        if both_shared > max_both_shared:
            max_both_shared = both_shared
            best_templ = templ
    print '%s %s %s %5.3lf'%(short,long,best_templ,max_both_shared)
    
align_in.close()
os.system('/bin/rm '+junk)
sys.exit(0)
