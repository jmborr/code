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
from utilities.codedir import scratchdir

#Global variables
SCroot='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
templd='/gpfs1/archive/adrian/dat/pdb_mar07' #template directory
pdbd=templd+'/CA'
pdbd_2=templd+'/aux'
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
junk=junkName()
alignments={}
rex=re.compile('1+')
trioIDs=('short','long','templ')
file_header='#C-C: constituve region not covered by template\n#V+C: variant region covered\n'

def getIsCovered(trio,whichone,L):
    """return a string of "1" and "0". "1" is that position is covered by related template"""
    short,long,templ=trio.split()
    if whichone=='short': which=short
    else: which=long
    templpdb=os.path.join(scratchdir,'spliceVariants/benchmarkByTMalign/out',short+'_'+long,which,
                          which+'rap3orienrev5s.pdb')
    if not os.path.exists(templpdb): return ''
    isCovered=['0']*L
    for l in open(templpdb,'r').readlines():
        if ' CA ' in l:
            isCovered[ int(l[22:26])-1 ]='1'
    return ''.join(isCovered)
                          
    
trios=open(scd+'/findCommonTemplate.dat','r').readlines()
remnant=len(trios)
for trio in trios:
    missingpdbflag=False
    if not remnant%100: print '%4d'%(remnant),'left'
    if trio[0]=='#': continue #trio is a comment
    headers={} #short,long, and templ headers
    (headers['short'], headers['long'], headers['templ'])=trio.split()
    pdbs={} #find associated pdb files to previously identified headers
    for trioID in ('short','long'): #we don't need pdb file of template
        pdbf=os.path.join(pdbd,headers[trioID]+'.pdb')
        if not os.path.exists(pdbf):
            pdbf=os.path.join(pdbd_2,headers[trioID]+'.pdb')
            if not os.path.exists(pdbf):
                sys.stderr.write(os.path.basename(pdb)+' missing from pdb_mar07\n')
                missingpdbflag=True
        pdbs[trioID]=pdbf
    if missingpdbflag: continue #skip to next trio
    #gather info on the alignment between short and long
    cmd='TMalign '+pdbs['short']+' '+pdbs['long']  #;os.system(cmd)#;Bye(cmd)
    alignments['short,long']=TMalignOut( os.popen(cmd) ).alg
    
    isAligned={ 'const,covered':{'short':'','long':''},
                'const,uncovered':{'short':'','long':''},
                'var,covered':{'short':'','long':''},
                'var,uncovered':{'short':'','long':''},
                }
    x=alignments['short,long'].isAlignedLines(algtype='aligned')
    isAligned['short,const']=x[1]
    isAligned['long,const']=x[2]
    #Bye(isAligned['short,const']+'\n'+isAligned['long,const'])
    L={'short':len(isAligned['short,const']), 'long':len(isAligned['long,const']) }
    isAligned['short,covered']=getIsCovered(trio,'short',L['short'])
    isAligned['long,covered']=getIsCovered(trio,'long',L['long'])
    #Bye(isAligned['short,covered']+'\n'+isAligned['long,covered'])
    for i in range(len(isAligned['short,const'])):
        x=int(isAligned['short,const'][i])
        y=int(isAligned['short,covered'][i])
        isAligned['const,covered']['short']+=`x*y` #1 for res constitutive and covered by template
        isAligned['const,uncovered']['short']+=`x*(1-y)`
        isAligned['var,covered']['short']+=`(1-x)*y`
        isAligned['var,uncovered']['short']+=`(1-x)*(1-y)`
    #Bye(isAligned['const,covered']+'\n'+isAligned['const,uncovered'])
    #Bye(isAligned['var,covered']['short']+'\n'+isAligned['var,uncovered']['short'])
    for i in range(len(isAligned['long,const'])):
        x=int(isAligned['long,const'][i])
        y=int(isAligned['long,covered'][i])
        isAligned['const,covered']['long']+=`x*y` #1 for res constitutive and covered by template
        isAligned['const,uncovered']['long']+=`x*(1-y)`        
        isAligned['var,covered']['long']+=`(1-x)*y`
        isAligned['var,uncovered']['long']+=`(1-x)*(1-y)`
    #Bye(isAligned['var,covered']['long']+'\n'+isAligned['var,uncovered']['long'])
    #output files
    #create single line with all four categories
    allcateg={'const,covered':'1','const,uncovered':'2','var,covered':'3','var,uncovered':'4'}
    inverallcat={'1':'C+C', '2':'C-C', '3':'V+C', '4':'V-C'}
    for key in ('short','long'):
        header=headers[key]
        mix=''
        for i in range(L[key]):
            if isAligned['const,covered'][key][i]=='1': mix+=allcateg['const,covered']
            elif isAligned['const,uncovered'][key][i]=='1': mix+=allcateg['const,uncovered']
            elif isAligned['var,covered'][key]=='1': mix+=allcateg['var,covered']
            else: mix+=allcateg['var,uncovered']
        #output chunk by chunk
        outf=os.path.join(scratchdir,'spliceVariants/benchmarkByTMalign/out',
                          headers['short']+'_'+headers['long'],
                          header,'region.dat')
        outp=open(outf,'w')
        outp.write(file_header)
        n=0
        while(n<L[key]):
            cat=mix[n]
            m=n+1
            while(m<L[key]):
                if mix[m]!=cat: break
                m+=1
            outp.write('%3d %3d %s\n'%(n+1,m,inverallcat[cat]))#shift +1 because of C-array issue
            n=m
    remnant=remnant-1

sys.exit(0)
