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

We find the different regions in the alignment corresponding to these categories, and output the length of the biggest region for each categorie.
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
rex=re.compile('1+')
outp=open('./region_big_sizes.dat','w')

trios=open(scd+'/findCommonTemplate.dat','r').readlines()
remnant=len(trios)
outp.write('#C-C: constitutive and covered (by the template) region\n')
outp.write('#V+C: variant and covered (by the template) region\n')
outp.write('#short long templ C+C C-C V+C V-C\n')
for trio in trios:
    missingpdbflag=False
    #if not remnant%100: print '%4d'%(remnant),'left'
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
    #output biggest chunk
    buf=trio.strip()+' '
    for category in ('const,covered','const,uncovered'):
        maxl=0
        for chunk in rex.findall( isAligned[category] ):
            l=len(chunk)
            if l>maxl: maxl=l
        buf+='%3d '%(maxl,)        
    for category in ('var,covered','var,uncovered'):
        maxl=0
        for type in ('short','long'): #dump both short and long on same category
            for chunk in rex.findall( isAligned[category][type] ):
                l=len(chunk)
                if l>maxl: maxl=l
        buf+='%3d '%(maxl,)
    outp.write(buf+'\n')
    remnant=remnant-1

sys.exit(0)
