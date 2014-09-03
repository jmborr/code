#!/usr/bin/python

"""
Given a short, long, and template structures aligned like this:

          (A)       |  (B)  |   (C)  |       (D)         | (C) | (D) |
long  --------------|-------|----------------------------|     |     |
templ --------------|       |--------|                   |-----|     |
short --------------|-------|        |                   |-----|-----|

We have four categories:

(A) "constitutive" regions covered by the template (Here
"constitutive" mean short and long are aligned by TM-align)

(B) "constitutive" regions not covered by the template. By
construction we tried to minimize these

(C) "variant" regions covered by the template (Here "constitutive"
mean short and long are not aligned by TM-align)

(D) "variant" regions not covered by the template

By scanning though all alignments, we calculate calculate histogram of
sizes of chunks for the four different types of combination of
constitutive/variant, covered/uncovered options

We also output histograms of chunks for constitutive or variant
fragments.

"""

import os,sys,re
from utilities.small_utilities import chomp,junkName,Bye
from TMalign.TMalignManageResults import TMalignOut
from seq.alignYangManageResults import gappedSeq

#Global variables
SCroot='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
templd='/gpfs1/archive/adrian/dat/pdb_mar07' #template directory
pdbd=templd+'/CA'
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
junk=junkName()
alignments={}
histograms={'constitutive':[0]*300,'variant':[0]*300,'cons+cov':[0]*300,'cons-cov':[0]*300,
       'var+cov':[0]*300,'var-cov':[0]*300,}#"cons-cov" means const. region not-covered by template
rex=re.compile('1+')
outf='./region_sizes.dat'

trios=open(scd+'/findCommonTemplate.dat','r').readlines()
remnant=len(trios)
for trio in trios:
    missingpdbflag=False
    #if not remnant%100: print '%4d'%(remnant),'left'
    print remnant
    if trio[0]=='#': continue #trio is a comment
    (short,long,templ)=trio.split()
    spdb=pdbd+'/'+short+'.pdb' ; lpdb=pdbd+'/'+long+'.pdb' ; tpdb=pdbd+'/'+templ+'.pdb'
    for pdb in (spdb,lpdb,tpdb):
        if not os.path.exists(pdb):
            sys.stderr.write(os.path.basename(pdb)+' missing from pdb_mar07\n')
            missingpdbflag=True
    if missingpdbflag: continue
    #gather info on the alignment between short and long
    cmd='TMalign '+spdb+' '+lpdb  #;os.system(cmd)#;Bye(cmd)
    alignments['short,long']=TMalignOut( os.popen(cmd) ).alg
    #print alignments['short,long'].presentAlignedFragments() ;Bye('')
    for chunk in alignments['short,long'].presentAlignedFragments():
        histograms['constitutive'][len(chunk)]+=1
    #print alignments['short,long'].extractGaps() ;Bye('')
    for chunk in alignments['short,long'].extractGaps():
        histograms['variant'][len(chunk)]+=1
    #align template to short and to long
    nconst=alignments['short,long'].nAligned()
    cmd='TMalign '+spdb+' '+tpdb+'  -L '+`nconst` #;os.system(cmd)#;Bye(cmd)
    alignments['short,templ']=TMalignOut( os.popen(cmd) ).alg
    cmd='TMalign '+lpdb+' '+tpdb+'  -L '+`nconst` #;os.system(cmd)#;Bye(cmd)
    alignments['long,templ']=TMalignOut( os.popen(cmd) ).alg
    #lines containing 1 or 0. "1" if the condition is true
    #find line containing 1 for different conditions
    isAligned={ 'const,covered':'', 'const,uncovered':'',
                'var,covered':{'short':'','long':''},
                'var,uncovered':{'short':'','long':''},
                }
    x=alignments['short,long'].isAlignedLines(algtype='aligned')
    isAligned['short,const']=x[1]
    isAligned['long,const']=x[2]
    #Bye(isAligned['short,const']+'\n'+isAligned['long,const'])
    isAligned['short,covered']=alignments['short,templ'].isAlignedLines(algtype='aligned')[1]
    isAligned['long,covered']=alignments['long,templ'].isAlignedLines(algtype='aligned')[1]
    #Bye(isAligned['short,covered']+'\n'+isAligned['long,covered'])
    for i in range(len(isAligned['short,const'])):
        x=int(isAligned['short,const'][i])
        y=int(isAligned['short,covered'][i])
        isAligned['const,covered']+=`x*y` #1 if residue is constitutive and covered by template
        isAligned['const,uncovered']+=`x*(1-y)`
        isAligned['var,covered']['short']+=`(1-x)*y`
        isAligned['var,uncovered']['short']+=`(1-x)*(1-y)`
    #Bye(isAligned['const,covered']+'\n'+isAligned['const,uncovered'])
    #Bye(isAligned['var,covered']['short']+'\n'+isAligned['var,uncovered']['short'])
    for i in range(len(isAligned['long,const'])):
        x=int(isAligned['long,const'][i])
        y=int(isAligned['long,covered'][i])        
        isAligned['var,covered']['long']+=`(1-x)*y`
        isAligned['var,uncovered']['long']+=`(1-x)*(1-y)`
    #Bye(isAligned['var,covered']['long']+'\n'+isAligned['var,uncovered']['long'])
    #update histograms
    for chunk in rex.findall( isAligned['const,covered'] ):        
        histograms['cons+cov'][len(chunk)]+=1
    for chunk in rex.findall( isAligned['const,uncovered'] ):
        histograms['cons-cov'][len(chunk)]+=1
    for type in ('short','long'):
        for chunk in rex.findall( isAligned['var,covered'][type] ):
            histograms['var+cov'][len(chunk)]+=1
        for chunk in rex.findall( isAligned['var,uncovered'][type] ):
            histograms['var-cov'][len(chunk)]+=1
    remnant=remnant-1
#print histograms
outp=open(outf,'w')
for key in histograms.keys():
    hist=histograms[key]
    outp.write('&'+key+'\n')
    for i in range(300):
        outp.write('%3d %4d\n'%(i,hist[i]))
    


    
     
   
    
