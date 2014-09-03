#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb,isglobbed
from utilities.codedir import projectsdir,codedir,nodecodedir,scratchdir
from utilities.small_utilities import chomp,todayDate,Bye
from jobs.job import job,pastry,genJobList

#check if there's a logfile with succesful outcome
def jobcompleted(logfilename):
    if os.path.exists(logfilename):
        if not os.system('tail -1 '+logfilename+'|grep "exit mode = 0" > /dev/null'): return 1
    return 0

#generate log file name
def logname(identifier,header,outd):
    branch='/'+header[1]+'/'+header
    outd2=deglobb(outd,header)
    logfile=outd2+'/'+identifier+header+'.log'
    return logfile

#check input files exists. It is assume an input line of the form
#[/bin/mv] [-f] input1 input2 inputN targetDir
def check_input_files(inpcmd,logfile):
    inputfine=1
    list=inpcmd.split()  ;del list[0]  ;del list[len(list)-1] #list of input files
    #cycle through all input files
    for input in list:
        if input[0]=='-': continue  #it's a flag, not an input file
        #if the input file does not exists, exit and report in the logfile
        #we don't use os.path.exists because this function does not support unix bloging
        if os.system('ls '+input+' &>/dev/null')!=0:
            inputfine=0
            if os.path.exists(logfile): #pre-existing logfile with same name                
                #if indicates success of previous job, then it doesn't matter that input
                #file does not exists. We will submit the job and then job.qsub will be
                #smart enough to NOT queue the job
                if os.system('grep "exit mode = 0" '+logfile): return 1
                else:
                    pastry('/bin/rm '+logfile)
            pastry('touch '+logfile)
            pastry('echo "input file '+input+' missing. Can not submit job" >> '+logfile)
    return inputfine

ih=inpHand('Usage: preditc_struct_list.py',
           ' -a _A_list (list of headers)',
           ' -b __header (single header, in place of -a option)',
           ' -c _AR_outd output directory where logfile will go (will substitute /x/xxxxx by appropriate)',
           ' -d __wallt walltime, in days (default: 7.0)',
           ' -e __xflags qsub may need these extra flags. Enclose them within double quotes',
           ' -f __inpcmd input command line (in double quotes, x/xxxxx globbing)',
           ' -g __execcmd executable command line, run from the temp dir (in double quotes)',
           ' -h __help outputs a help message (no arguments needed)',
           ' -i __outcmd output command line  (in double quotes)',
           ' -j __submode submission mode of each job (sub,bsub,rsub) def=qsub',
           ' -k __memlimit minimum memory reservation (def=500 (Mb))',
           ' -l __filesize minimum disk space reservation (def=1000 (Mb))',
           ' -m __shared0 libline "dir1:tbz21,dir2:tbz22", globbing allowed (def: None)',
           ' -n _R_identifier one word for jobname and logfile name. Will append xxxxx',
           ' -o __fake do we create qsub script but do not send the job? (def:no)',
           ' -p __ngroup send N jobs to a single node to be run concurrently (def:1)',
           ' -q __joblist generate a current list of jobs in queing server (def:yes)'
           )
ih.parse(locals(),sys.argv)
    
if not submode: submode='qsub'
if not wallt: wallt='7.0'
if not memlimit: memlimit='500'
if not filesize: filesize='100'
if not header: header=None
if not shared0: shared0=None
if not ngroup:ngroup=1
else: ngroup=int(ngroup)

#joblist avoids calling qstat for every single job that we qsub
if not joblist or joblist[0] in ('y','Y'):
    joblist=genJobList()
else: joblist=[]

#print 'list=',list ; print 'outd=',outd ; print 'inpcmd=',inpcmd ; print 'outcmd=',outcmd

j=job(exe='generic_job.py', exed=nodecodedir+'/python/combo_jobs', args='' )
jj=job(exe='generic_job.py', exed=nodecodedir+'/python/combo_jobs', args='' )

if not list and not header: ih.abort() #we need at least on the two inputs
if header: listl=[header,] #create a single item list
else:      listl=chomp(open(list,'r').readlines())

#We need scaping of ",`,$ when inside double quotes, because unfortunately, these
#characters will be interpreted.
#example: junk.py -a "echo "Hello" "  #confusion with the "
#example: junk.py -a "echo $wd"       #$wd will be substituted with whatever value
#example; junk.py -a "wd=`cat junk`"  #wd will be initialized, ie, the argument within quotes
#                                     #will not be passed literally, but interpreted

execcmd=execcmd.replace('\\','\\\\')  #this substitution MUST come in the first place!
execcmd=execcmd.replace('"','\\"')
execcmd=execcmd.replace('`','\`')
execcmd=execcmd.replace('$','\$')

#cicle through all headers in list, prepare a job for each header
ithgroup=0
ninthebunch=0   #current number of jobs assigned to the current group of jobs
groupexeccmd='' #line containing ngroup scripts, one for each of the ngroup jobs
bunchlist=''
for header in listl:
    branch='/'+header[1]+'/'+header
    jobname=identifier+header
    logname=jobname+'.log'
    unixcommands=''
    if inpcmd : unixcommands=' -a "'+inpcmd +'" '
    if execcmd: unixcommands=unixcommands+' -b "'+execcmd+'" '
    if outcmd : unixcommands=unixcommands+' -c "'+outcmd +'" '
    j.args=unixcommands
    j.args=deglobb(j.args,header) #deglobb switches xxxxx to header, and yyyyy to xxxxx
    if shared0 and isglobbed(shared0): j.shared=deglobb(shared,header)

    outdH=deglobb(outd,header) #directory where logname will go
    logf=outdH+'/'+logname     #full file name of log file
    if jobcompleted(logf): #there's a log file with a successful outcome
        sys.stdout.write('COMPLETED: '+jobname+'\n')
        continue
    if not os.path.exists(outdH): pastry('/bin/mkdir -p '+outdH)
    allclear=True #flag input files are not missing
    if inpcmd: #check that input files do exist
        inpcmdH=deglobb(inpcmd,header)
        if not check_input_files(inpcmdH,logf): allclear=False
    if allclear:
        if ngroup>1: #we're bundling jobs in bundles of ngroup jobs
            #create a script for this single job and add to the list
            groupexeccmd+=j.scriptGen(jobname,outdH,mem_limit=memlimit,libRev='',submode='sub')+' ; '
            ninthebunch+=1
            bunchlist+=jobname+'  '
            if ninthebunch==ngroup or header==listl[-1]: #bundle is full or last header
                ithgroup+=1
                jj.args=' -b "'+groupexeccmd+'" '
                jobname=identifier+'G'+`ithgroup`
                outdir=scratchdir+'/qsub/'+todayDate()
                if shared0 and not isglobbed(shared0): jj.shared=shared0 #same lib for all jobs
                sys.stdout.write('\nlist of jobnames included in '+jobname+':\n'+bunchlist+'\n')
                getattr(jj,submode)(jobname,outdir,wallt=wallt,mem_limit=memlimit,
                                   file_size=filesize,extraflags=xflags,fake=fake)
                ninthebunch=0  #begin another bundle
                groupexeccmd=''
                bunchlist=''
        else:
            getattr(j,submode)(jobname,outdH,wallt=wallt,mem_limit=memlimit,
                               file_size=filesize,extraflags=xflags,fake=fake,
                               joblist=joblist)

sys.exit(0)
