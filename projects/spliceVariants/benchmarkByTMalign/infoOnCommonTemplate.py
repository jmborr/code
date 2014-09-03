#!/usr/bin/python

import os,sys
from utilities.small_utilities import chomp,junkName,Bye
from TMalign.TMalignManageResults import TMalignOut
from seq.alignYangManageResults import gappedSeq

SCroot='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
templd='/gpfs1/archive/adrian/dat/pdb_aug06.apr06' #template directory
pdbd=templd+'/CA'
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'

junk=junkName()

print '#          A               B           C      D'
print '#short xxxxxxxxxxxxx----------------  xxxx--------'
print '#                E             F             G                         H'
print '#long  xxxxxxxxxxxxxxxxxxxxxx-------  xxxxxxxxxxxxxxxxxx----------------------------'
print '#'
print '#A: constitutive   & short/templ (residues aligned in both in short/long and short/templ'
print '#C: short_variant  & short/templ'
print '#E: constitutive   & long/templ'
print '#G: long_variant   & long/templ'
print '#'
print '#C1=constitutive  C2=short  C3=nshort  C4=C1/C3  C5=shared[short]  C6=C5/C1  C7=var[short]  C8=C7/(C3-C1)  C9=long  C10=nlong  C11=C1/C10  C12=shared[long]  C13=C12/C1  C14=var[long]  C15=C14/(C10-C1)'
print '#'
print '#C1   C2   C3   C4   C5   C6   C7   C8   C9   C10  C11  C12  C13  C14  C15'

for line in open(scd+'/findCommonTemplate.dat','r').readlines():
    if line[0]=='#': continue #line is a comment
    (short,long,templ)=line.split()
    spdb=pdbd+'/'+short+'.pdb' ; lpdb=pdbd+'/'+long+'.pdb' ; tpdb=pdbd+'/'+templ+'.pdb'
    pdbs={short:spdb,long:lpdb}
    #gather info on the alignment between short and long 
    cmd='TMalign '+spdb+' '+lpdb+' > '+junk  ;  os.system(cmd)
    alg=TMalignOut(junk)
    constitutive=alg.nAligned               #number of aligned residues
    alg.initIsAlignedArray() #initialize alg.isAlined dictionary
    arrayOfAligned={short:alg.isAligned[1],long:alg.isAligned[2]}
    nshort=len( arrayOfAligned[short] ) #length of short sequence
    nlong=len( arrayOfAligned[long] )   #lenthg of long sequence
    lengths={short:nshort,long:nlong}
    #compare to template
    shared={}
    var={}
    for header in (short,long):
        pdb=pdbs[header]
        cmd='TMalign '+pdb+' '+tpdb+'  -L '+`constitutive`+' > '+junk  ;  os.system(cmd)
        alg2=TMalignOut(junk)
        alg2.initIsAlignedArray()
        shared[header]=0 #num. residues of header aligned to both the template and the other header
        var[header]=0    #num. residues of header aligned to template but not to other header
        for i in range( len(arrayOfAligned[header]) ):
            shared[header]+=arrayOfAligned[header][i]     * alg2.isAligned[1][i]
            var[header]   +=(1-arrayOfAligned[header][i]) * alg2.isAligned[1][i]
    
    if nshort==constitutive: #treat this special case
        a=0.0
    else:
        a=(1.0*var[short])/(nshort-constitutive)
        
    if nlong==constitutive: #treat this special case
        b=0.0
    else:
        b=(1.0*var[long])/(nlong-constitutive)
       
    print '%3d %s %3d %5.3lf %3d %5.3lf %3d %5.3lf %s %3d %5.3lf %3d %5.3lf %3d %5.3lf'%(constitutive,short,nshort,(1.0*constitutive)/nshort,shared[short],(1.0*shared[short])/constitutive,var[short],a,long,nlong,(1.0*constitutive)/nlong,shared[long],(1.0*shared[long])/constitutive,var[long],b)

sys.exit(0)
