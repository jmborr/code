#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.codedir import projectsdir,codedir
from utilities.small_utilities import chomp
from jobs.job import job,pastry

ih=inpHand('Usage: preditc_struct_list.py',
           ' -a _A_list (list of headers)',
           ' -b __header (single header, in place of -a option)',
           ' -c _AR_outd output directory where logfile will go (will substitute /x/xxxxx by appropriate)',
           ' -d __queue submission queue(default=batch)',
           ' -e __xflags qsub may need these extra flags. Enclose them within double quotes',
           ' -f __inpcmd input command line (in double quotes, x/xxxxx globbing)',
           ' -g __execcmd executable command line, run from the temp dir (in double quotes)',
           ' -h __help outputs a help message (no arguments needed)',
           ' -i __outcmd output command line  (in double quotes)',
           ' -j __submode submission mode of each job (sub,bsub,rsub) def=qsub',
           ' -k __memlimit minimum memory reservation (def=500 (Mb))',
           ' -l __filesize minimum disk space reservation (def=1000 (Mb))',
           ' -m __shared libline "dir1:tbz21,dir2:tbz22", globbing allowed (def: None)',
           ' -n _R_identifier one word for jobname and logfile name. Will append xxxxx',
           )
ih.parse(locals(),sys.argv)
    
if not submode: submode='qsub2'
if not queue: queue='batch'
if not memlimit: memlimit='500'
if not filesize: filesize='1000'
if not header: header=None
if not shared: shared=None

#print 'list=',list ; print 'outd=',outd ; print 'inpcmd=',inpcmd ; print 'outcmd=',outcmd

#generate log file name
def logname(identifier,header,outd):
    branch='/'+header[1]+'/'+header
    outd2=deglobb(outd,header)
    logfile=outd2+'/'+identifier+header+'.log'
    return logfile

#check input files exists. It is assume an input line of the form
#[/bin/mv] [-f] input1 input2 inputN targetDir
def check_input_files(header,identifier,inpcmd,outd):
    inputfine=1
    list=inpcmd.split()  ;del list[0]  ;del list[len(list)-1]
    #cycle through all input files
    for input in list:
        if input[0]=='-': continue  #it's a flag, not an input file
        #if the input file does not exists, exit and report in the logfile
        print 'input=',input
        if os.system('ls '+input+' 2>/dev/null >/dev/null'):
            inputfine=0
            logfile=logname(identifier,header,outd)
            pastry('mkdir -p '+outd)
            if os.path.exists(logfile): #pre-existing logfile with same name                
                #if indicates success of previous job, then it doesn't matter that input
                #file does not exists. We will submit the job and then job.qsub will be
                #smart enough to NOT queue the job
                if os.system('grep "exit mode = 0" '+logfile): return 1
                else:
                    pastry('/bin/rm '+logfile)
            pastry('touch '+logfile)
            pastry('echo "input file '+input+' missing. Can not submit job" >> '+logfile)
            inputfine=0
    return inputfine


j=job(exe='generic_job2.py', exed='/gpfs1/active/jose/code/python/combo_jobs', args='' )

if not list and not header: ih.abort() #we need at least on the two inputs
if header: listl=[header,] #create a single item list
else:      listl=chomp(open(list,'r').readlines())

#cicle through all headers in list, prepare a job for each header
for header in listl:
    branch='/'+header[1]+'/'+header
    j.args=' -x '+identifier+header

    #We need scaping of ",`,$ when inside double quotes, because unfortunately, these
    #characters will be interpreted.
    #example: junk.py -a "echo "Hello" "  #confusion with the "
    #example: junk.py -a "echo $wd"       #$wd will be substituted with whatever value
    #example; junk.py -a "wd=`cat junk`"  #wd will be initialized, ie, the argument within quotes
    #                                     #will not be passed literally, but interpreted

    if inpcmd :
        inpcmd=inpcmd.replace('"','\\"')
        unixcommands=' -a "'+inpcmd +'" '

    if execcmd:
        execcmd=execcmd.replace('"','\\"')
        unixcommands=unixcommands+' -b "'+execcmd+'" '

    if outcmd :
        outcmd=outcmd.replace('"','\\"')
        unixcommands=unixcommands+' -c "'+outcmd +'" '

    #escaping of several special characters
    unixcommands=unixcommands.replace('`','\`')
    unixcommands=unixcommands.replace('$','\$')

    #print unixcommands ; sys.exit(1)
    j.args=j.args+unixcommands
    j.args=j.args+' -d '+outd
    j.args=deglobb(j.args,header)

    print j.args;sys.exit(1)
    
    if shared: j.shared=deglobb(shared,header)

    outdH=deglobb(outd,header)
    inpcmdH=deglobb(inpcmd,header)
    if check_input_files(header,identifier,inpcmdH,outdH) >0: #submit job if input files do exist
        getattr(j,submode)(queue=queue,mem_limit=memlimit,file_size=filesize,extraflags=xflags)
    else:
        pastry('echo "exit mode = 1" >> '+logname(identifier,header,outdH))

sys.exit(0)
