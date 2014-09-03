#!/usr/bin/python
'''
Filtering a set of structures by TM score. It contains several steps

(1) Read a list of structures, then shuffle the list

(2) Divide the list into N chunks of size about 200 structures. Send
    every chunk as a job to be filtered. The filtered chunk will
    contain a set of structures such that no two structures have a TM
    score bigger than a given cutoff.

(3) Filter one chunk agains other chunk. We have N*(N-1)/2
    possibilities, which amounts to same number of jobs. For example,
    when we filter chunk 3 against chunk 8, then we obtain a filtered
    list from chunk 3 such that no structure in this list has a structural
    homologous to structures in chunk 8.

(4) For chunk "i", gather all filtered lists from all (i,j) j=i+1,..,N
    possible pairings. Then compare all this i+1,..,N lists and keep
    only those structures showing in all the lists, forming a final
    filter list for chunk "i". This ensures that these structures have
    no structural homologos in chunk i+1 through chunk N.  We do this
    for chunk 1 to chunk N-1, merging all the final filtered
    lists. Finally, we also add the filtered list from chunk N
    obtained in step 3

'''
import os,sys,re
from jobs.job import pastry,genJobList
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName
from jobs.job import job
from TMalign.TMalignManageResults import TMalign
from random import shuffle

ih=inpHand('Usage: filterByTMcutoff.py',
           ' -a __tmco TM cut off (def: 0.45)',
           ' -b _RA_inf input file with list of structures',
           ' -c _RA_outf output file',
           ' -d _RA_shared "/tmp/jose/structures:structures.tbz2" line for library of structures to be imported to each node upon job submission',
           ' -e _R_prefix prefix (not more that three characters)',
           ' -f __randomize run routine to randomize the list of structures and split the jobs (def: yes)'
           ).parse(locals(),sys.argv)

joblist=genJobList()

if not tmco:tmco='0.45'
randomize=0
if not randomize or randomize[0] in ('y','Y'): randomize=1
prefix=prefix[0:3] #take first three characters

strd=re.compile('([/\w]+)').match(shared).group(1) #directory containing the structures

#temporary directory
currd=os.getcwd()
tmpd=currd+'/filterByTMalignDir_'+prefix
pastry('/bin/mkdir -p '+tmpd)


#divide list of structures into several chunks
if randomize:
    l=open(inf,'r').readlines()
    shuffle(l) #mix all structures by random shuffling
    chunksize=200
    x=len(l)/chunksize
    chunksize+=(len(l)-chunksize*x)/x #refine chunksize
    b=0
    e=0
    for i in range(x-1):
        e=b+chunksize
        f=tmpd+'/'+prefix+'%02d'%(i)+'.0'
        open(f,'w').write(''.join(l[b:e]))
        b=e
    f=tmpd+'/'+prefix+'%02d'%(x-1)+'.0'    
    open(f,'w').write(''.join(l[b:]))

#send chunks to be filtered each one with filterByTMalign.py
jobID={} #store job ID's
for i in range(x):
    jobname=prefix+'%02d'%(i)
    fin =tmpd+'/'+jobname+'.0'  #one of the chunks
    fout=tmpd+'/'+jobname+'.1' #filtered chunk
    jobID[jobname]=job(exe='filterByTMalign.py', exed=currd, shared=shared, args='-a '+tmco+' -b '+fin+' -c '+fout+' -d '+strd).qsub(jobname,tmpd,wallt='0.249',ddisk='1500',mem_limit='200',joblist=joblist)
    if jobID[jobname]=='':Bye(jobname+' could not be submitted')
    print jobID[jobname]


#cross-over each filtered chunk
jobID2={}
for i in range(x-1):
    jobnamei=prefix+'%02d'%(i)
    fini=tmpd+'/'+jobnamei+'.1'
    for j in range(i+1,x):
        jobnamej=prefix+'%02d'%(j)
        finj=tmpd+'/'+jobnamej+'.1'
        jobnameij=prefix+'%02d'%(i)+'.'+'%02d'%(j)
        fout=tmpd+'/'+jobnameij
        #job will not run until chunks have been filtered
        dependencies=[]
        if jobID[jobnamei]!='COMPLETED': dependencies.append(jobID[jobnamei])
        if jobID[jobnamej]!='COMPLETED': dependencies.append(jobID[jobnamej])
        jobID2[jobnameij]=job(exe='filterByTMalign2.py', exed=currd, shared=shared, args='-a '+tmco+' -b '+fini+' -c '+finj+' -d '+fout+' -e '+strd).qsub(jobnameij,tmpd,wallt='0.249',ddisk='1500',mem_limit='200',joblist=joblist,depend=dependencies)
        if jobID2[jobnameij]=='':Bye(jobnameij+' could not be submitted')

#resolve all cross-overs Run filterByTMalign3.py, which contains the code below
#lall=[]
#for i in range(x-1):
#    #merge all crossed-over filtered chunks for a given "i"
#    rooti='%02d'%(i)
#    l=open(tmpd+'/'+prefix+rooti+','+'%02d'%(i+1),'r').readlines()
#    l=lfil
#    for j in range(i+2,x):
#        lfil=[]
#        rootj='%02d'%(j)
#        fij=tmpd+'/'+prefix+rooti+'.'+rootj
#        for s in open(fij,'r').readlines():
#            if s in l: lfil.append(s)
#        l=lfil #keep only those common structures
#    lall+=lfil
#lall+=open(tmpd+'/'+prefix+'%02d'%(x-1),'r').readlines()
#open(outf,'w').write(''.join(lall))

#finish with cleanup and exit with success
#os.system('/bin/rm -rf '+tmpd)
sys.exit(0)
