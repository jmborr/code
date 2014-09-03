#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq
from utilities.small_utilities import chomp
from jobs.job import job,genJobList

ih=inpHand('Usage: preditc_struct_list.py',
           ' -a _A_subj list of successfully submitted jobs (log file)',
           ' -b __memlimit minimum memory reservation (def=500 (Mb))',
           ' -c __nidp number of independent runs for each target (def=1) to be run on different nodes',
           ' -d __hd header (in place of flag -l)',
           ' -l _A_listf (list of five-letter codes, in place of flag -d)',
           ' -m _R_move list of output files to be moved (in double quotes, globbing allowed)',
           ' -o _AR_outd output directory (substitute /x/xxxxx by appropriate)',
           ' -p _A_prospector prospector job',
           ' -q __queue submission queue(default=\'batch\')',
           ' -t _A_tasser tasser job         (leave bare header-related flag at the end!)',
           ' -u __submode submission mode of each job (sub,qsub,rsub,info) def=qsub',
           ' -v __xflags submission may need these extra flags. Enclose them within double quotes',
           ' -x _R_identifier one word for job-name purposes and for logfile name. Will append xxxxx',
           ' -z _A_spicker spicker job       (leave bare header-related flag at the end!)'
           )
ih.parse(locals(),sys.argv)

if not submode: submode='qsub'
if not queue: queue='batch'
if not memlimit: memlimit='500'
if subj and (not os.path.exists(subj)): os.system('touch '+subj)
if not nidp: nidp=1 #single run for each header
else: nidp=int(nidp)
if not listf and not hd: ih.abort()
#prepare common flags for all jobs
args0=''
#retain only the common flags that are actually passed. These flags
#lack the argument of the header flag
common={'p':prospector,'t':tasser,'z':spicker}
for flag in common.keys():
    if not common[flag]: del common[flag]
#instantiate a template job
j=job(exe='preditc_strutc.py', exed='/gpfs1/active/jose/code/python/combo_jobs', args='' )
if listf:
    if not os.path.exists(listf):
        sys.stderr.write('ERROR: no '+os.basename(listf)+' found in '+os.dirname(listf))
        stdout.write('program exiting')
        sys.exit(0)

#store the names of the jobs running in the cluster
jl=genJobList()

#submit one job at a time, change j.args for each job
if hd: listl=[hd,] #create a list out of the single header
else: listl=chomp(open(listf,'r').readlines())
                  
for header in listl: #chomp removes '\n'
    subd=header[1]
    branch='/'+subd+'/'+header
    for iidp in range(nidp):
        if nidp > 1: runID='%02d' % (iidp)
        else: runID=''
        j.args=''
        
        for flag in common.keys():
            j.args=j.args+' -'+flag+' "'+common[flag]+'"'
            
        j.args=j.args+' -x '+identifier+runID+header #set identifier for jobname
        j.args=j.args+' -o "'+outd+'"'

        if nidp>1: j.args=j.args+' -a '+`iidp`+' ' #pass iidp, if non-zero
        j.args=j.args+' -m "'+move+'"'             #pass the list of files to move        
        j.args=re.compile('/x/xxxxx').sub(branch,j.args) #substitute  /x/xxxxx by branch        
        j.args=re.compile('xxxxx').sub(header,j.args)    #substitute  xxxxx by header   

        #print 'j.args=',j.args

        #getattr returns reference to object method (it is also an attribute of
        #jobs j) associated to string submode. Then, () obligues execution
        if getattr(j,submode)(queue=queue,mem_limit=memlimit,extraflags=xflags,joblist=jl):
            if subj: os.system('echo '+header+' >>'+subj)
