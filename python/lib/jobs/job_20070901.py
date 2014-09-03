#!/usr/bin/python

"""implementation of class 'job'

"""

import os,re,sys
from utilities.codedir import codedir,scratchdir,nodecodedir #we have write permissions to codedir
from utilities.small_utilities import chomp,Bye,todayDate #remove '\n' of strings and lists
from time import sleep
from updateLib.updateLib import resolveLastRevision

sharedexe=nodecodedir+'/python/lib/jobs/shared.py' #script to install/remove temporary libraries
updnodelibexe=nodecodedir+'/python/updateLib/updateLib.py' #script to update nodecodedir
updnodelibexeII='/gpfs1/active/jose/nodelibmirror/code/python/updateLib/updateLib.py'

#####################################################
##########  O B S O L E T E    C O D E   ############
#####################################################
#add this directory in related search_in_dirs.list file
#def add_dir_for_backup(dir,depth=1):
#    projectD,subD=re.compile('jose/([\w,\.]+)/*([\w,\./]*)').search(dir).groups()
#    if not subD: subD='.'
#    searchf=scratchdir+'/'+projectD+'/search_in_dirs.list'
#    #check if the directory is already marked for backup
#    lines=chomp(os.popen('grep '+subD+' '+searchf).readlines())
#    if lines: #find the maximal depth scheduled for search of recent files
#        maxdepth=1
#        for line in lines:
#            depth2=int(line.split()[1])
#            if depth2>maxdepth:maxdepth=depth2
#        if maxdepth>depth: return None #nothing to do
#    pout=open(searchf,'a+').write(subD+' '+`depth`+'\n')    
#
#return the name of a log file, and handle special cases
#def gen_lf(logrequest,outd,identifier='job'):    
#    if logrequest=='_ASK_' or logrequest=='_LOG_' or not logrequest:
#        lf=identifier+'.log'
#        if logrequest=='_ASK_':
#            os.system('echo '+outd+'/'+lf)
#            sys.exit(0) #signal job should stop
#        else:
#            lf=os.getcwd()+'/'+lf #local to current working directory
#    else: stderr.write('Error: logrequest either of \'_ASK_\' or \'_LOG_\'')
#    return lf

def PBSdirectives(args):
    """include necessary PBS directives in script"""
    options=['-d /tmp',
             '-r n -c n -m n -j oe',
             '-N '+args['jobname'],
             '-o '+args['logf'],             
             #'- '+args[''],
             ]
    resource_list=['mem='+args['mem_limit']+'mb',
                   'file='+args['file_size'],
                   'ddisk='+args['ddisk'],
                   'walltime='+args['wallt'],
                   ]
    options.append( '-l '+','.join(resource_list) )
    return '#PBS '+'\n#PBS '.join(options)

def checkInput(inpl):
    """check that a list of input files exists.

    Returns a list of non-existant files. Otherwise returns empty string."""
    nonexistent=[]
    for inpf in inpl:
        if not os.path.exists(inpf): nonexistent.append(inpf)
    if nonexistent:
        return 'Input files:\n'+'\n  '.join(nonexistent)+'\ndo not exist'
    return ''
    
def add_job_dependencies(depend):
    """add job dependencies. Assume \"afterok\" type only"""
    if isinstance(depend,list):
        d='-W depend='
        for jobid in depend: d+='afterok:'+jobid+','
        return d[:-1] #remove trailing comma
    return ''

def formatWallTime(wallt):
    """assume wallt is a string with the number of days the job will run"""
    d=float(wallt)
    h=24*d ; m=(h-int(h))*60 ; s=(m-int(m))*60
    h=int(h) ; m=int(m) ; s=int(s)
    return '%d:%02d:%02d'%(h,m,s)
            
def genJobList():
    """generate a list of my jobnames in the cluster. Iterate if qstat call fails"""
    nit=10
    reply=[]
    while len(reply)<10: #assume qstat fails if less that 10 jobs on cluster
        reply=chomp(os.popen('qstat|expand|tr -s \' \'|cut -d\' \' -f 1,2,5').readlines())
        nit+=1
        if nit>10: break
    return reply
        
def abort(message,exit_code=1):
    sys.stderr.write(message+'\n')
    sys.exit(exit_code)
    
def pastry(command):
    if os.system(command):#unsuccessful command
        sys.stderr.write('ERROR while executing: '+command+'\n')
        sys.exit(1) #exit the calling (child) process
    return '0'

class job:
    """class to send jobs to PBS

    Attributes:
    name: name of the job
    exe: commands to execute
    inpl: list with input files
    outl: list with output files
    outd: ouptut directory for output files
    shared: libraries to be imported to /tmp
    time: some time stamp for the job
    file: filename if we want to save the job as bytecode
    story: an explanation of the job
    """
#    ####################################################
#    #########  O B S O L E T E    C O D E   ############
#    ####################################################
#    #absolute path to the logfile associated with this job. From this, we
#    #can produce a job name for queueing purposes (subname)    
#    def logfile(self):
#        return chomp(os.popen(self.exed+'/'+self.exe+' '+self.args+' -y _ASK_').readline())
#    
#    def subname(self,logfilename=''):
#        if not logfilename: logfilename=self.logfile()
#        return os.path.basename(logfilename).replace('.log','')
#    
#    def bsub(self,queue='standard',tmpd='/tmp/jose',zzz=2,extraflags='',
#             mem_limit=300):
#        pass #obsolete method
#    
#    def rsub(self,node='',tmpd='/tmp/jose',zzz=5,extraflags=''):
#        "optional args:tmpd='/tmp/jose',zzz=5,node=get_free_node()"
#        pass #obsolete method

    fields=('name','exe','inpl','outl','outd','shared','time','file','story')
    def __init__(self,name,exe,inpl=[],outl=[],outd='',shared='',time='',
                 file='',story=''):
        "optional args:inpl=[],outl=[],outd='',shared='',story=''"
        self.error=checkInput(inpl)
        self.name=name
        self.exe=exe
        self.inpl=inpl
        self.outl=outl
        self.outd=outd
        self.shared=shared
        self.story=story
        self.scriptf=''
        self.time=time
        self.file=file
        if self.error: sys.stderr.write(self.error+'\n')

    def __repr__(self):
        return "%(name)s instantiated on %(time)s" % vars(self)
    
    def info(self):
        line=''
        for field in self.fields:
            try: #it may be that instantiated object lacks a particular field
                attr=getattr(self,field)
            except AttributeError:
                continue #jump to next for iteration
            if isinstance(attr,list): attr=' '.join(attr)
            line=line+field.upper()+' = '+attr+'\n'
        return line

    def report_status(self,jobn,joblist=[]):
        if not joblist: joblist=genJobList()
        for joblisted in joblist:
            if jobn in joblisted:
                jobid,name,status=joblisted.split()
                if status!='C': return jobid
        return ''
        
    def completed(self,logfilename):
        if os.path.exists(logfilename):
            if not os.system('tail -1 '+logfilename+
                             '|grep "exit mode = 0" > /dev/null'):
                #print 'job is already completed !'
                return 'yes'
        return ''

    def updateNodeLib(self,libRev):
        """Run updateLib.py script to update node library. We pass a revision file"""
        foo='updateNodeLib () {\n'
        foo+='#update /library/jose\n'
        foo+=' e=1\n' #initialize exit status as error
        #transient_cmd is a command we want executed for some time
        transient_cmd= '   /bin/cp '+updnodelibexeII+' '+updnodelibexe+'\n'
        transient_cmd+='   find /library/jose/code/python -name "__init__.pyc" -exec /bin/rm {} \;'
        if libRev:
            foo+=' if [ -e '+updnodelibexe+' ]; then\n'
            foo+=transient_cmd+'\n'
            foo+='   '+updnodelibexe+' -a '+libRev +' && e=0\n'
            foo+=' else\n'
            foo+='   /bin/cp -r /gpfs1/active/jose/nodelibmirror/* /library/jose/ \n'
            foo+='   e=0\n'
            foo+=' fi\n'
        else:
            foo+='echo "No node libraries to update" && e=0\n'
        foo+='return $e\n'
        foo+='}'
        return foo
        
    def shared_temporal_libraries(self):
        """Run the executable that imports/deletes libraries to /tmp"""
        foo='shared_temporal_libraries () {\n'
        foo+='#bring or delete necessary temporary libraries\n'
        foo+=' e=1\n' #initialize exit status as error
        if self.shared:
            cmd=sharedexe+' -a "'+self.shared+'" -b $$ -c $1' #$1== _PREPARE_ or _CLEANUP_
            foo+=' if '+cmd+'; then e=0;fi\n'
        else:
            foo+=' echo "no shared temporal libraries to $1" && e=0\n'
        foo+=' return $e\n'
        foo+='}'
        return foo
        
    def cleanup_exit(self,submode='qsub'):
        """remove working directory, and maybe imported libraries"""
        foo='cleanup_exit () {\n'
        foo+=' shared_temporal_libraries _CLEANUP_\n'
        foo+=' hostname && echo $wd\n'
        foo+=' basename $LOCAL_LOG\n'
        foo+=' echo "exit mode = $1"\n'
        foo+=' /bin/mkdir -p `dirname $REMOTE_LOG`\n'
        if submode=='sub':
            foo+=' /bin/mv $LOCAL_LOG $REMOTE_LOG\n'
        foo+=' /bin/rm -rf $wd\n'
        foo+=' exit $1\n'
        foo+='}'            
        return foo
    
    def scriptGen(self,tmpd='/tmp/jose',libRev='last',submode='qsub',
                  redirect=1,PBSoptions=''):
        """Generate a script to be used with qsub or sub 
        wallt=\7.0\' walltime, in days (estimated time to complete the
        job) tmpd=\'/tmp/jose\' directory where generic subdir will be
        created zzz=1 time in seconds between consecutive jobs
        submissions"""
        jobname=self.name
        outdir=self.outd
        qsubdir=scratchdir+'/qsub/'+todayDate() #subdirectory to deposit the script
        if not os.path.exists(qsubdir):  pastry('/bin/mkdir -p '+qsubdir)
        script=qsubdir+'/'+jobname+'.sh' #full script file name

        if len(jobname) > 15:
            sys.stderr.write('Error: job name '+jobname+' cannot exceed 15 characters')
            return ''
        if not os.path.exists(outdir): os.system('/bin/mkdir -p '+outdir)
        buf=''
        ulimit=int(float(mem_limit)*1024) #maximum resident memory size (Kb) to prevent swapping
        wd=tmpd+'/${PBS_JOBID}'
        #wd=tmpd+'/'+ re.compile('\W').sub('',self.name) +'_$$' #working directory
        logname=jobname+'.log'
        local_log=wd+'/'+logname
        remote_log=outdir+'/'+logname
        buf=  '#!/bin/bash\n\n'
        buf+= PBSoptions+'\n\n'
        buf+= '#bash function to update library\n'
        buf+= self.updateNodeLib(libRev)+'\n\n'
        buf+= '#bash function to import temporary libs\n'
        buf+= self.shared_temporal_libraries()+'\n\n'
        buf+= '#bash function to clean exit\n'
        buf+= self.cleanup_exit(submode=submode)+'\n\n'
        buf+= 'echo "'+script+'"\n' #write script name withing script body
        buf+= 'hostname\n'          #node where job will be run
        buf+= 'echo $PBS_JOBID\n'
        buf+= 'ulimit -m '+`ulimit`+' #maximum memory\n'
        buf+= 'source ~/.bash_profile >/dev/null #environment variables\n'
        buf+= 'wd='+wd+' #working directory\n'
        buf+= '/bin/mkdir -p $wd\n'
        buf+= 'export LOCAL_LOG="'+local_log+'"\n'
        buf+= '/bin/touch $LOCAL_LOG\n'
        if submode=='sub' and redirect:
            buf+='exec &> $LOCAL_LOG #redirect STODOUT, STDERR to LOCAL_LOG\n' 
        buf+= 'export REMOTE_LOG="'+remote_log+'"\n'

        but+= '#clean up old log file\n'
        buf+= 'if [ -f $REMOTE_LOG ]; then\n' 
        buf+= ' /bin/rm -f $REMOTE_LOG\n'
        buf+= 'fi\n\n'

        buf+= 'trap "cleanup_exit 1" TERM #in case of killing job\n\n'

        buf+= '#update node code library && import libraries\n'
        buf+= 'if !('
        buf+= 'updateNodeLib && ' 
        buf+= 'shared_temporal_libraries _PREPARE_'
        buf+= ');then\n'
        buf+= '  cleanup_exit 1\n'
        buf+= 'fi\n\n'
        
        buf+= '/bin/cp '+' '.join(self.inpl)+' $wd #bring input files\n' 
        buf+= 'cd $wd\n\n'
        buf+= '#Test command success\n'
        buf+= 'exs=0 #variable holding script exit status\n'
        buf+= 'if !('
        buf+= self.exe
        buf+= ');then\n'
        buf+= '  exs=1\n'
        buf+= 'fi\n\n'
        buf+= '#move even partial results (exs=1)\n'
        buf+= '/bin/mv '+' '.join(self.outl)+' '+outdir+'\n'
        buf+= 'cleanup_exit $exs'

        open(script,'w').write(buf)
        pastry('chmod u+x '+script)

        return script

    def sub(self,jobname,outdir,tmpd='/tmp/jose',mem_limit='500',zzz=1.0,libRev='last',
            fake='',redirect=1,joblist=[]):
        logname=jobname+'.log'
        logf=outdir+'/'+logname  #;Bye(logf)
        if not self.completed(logf):
            libRev=resolveLastRevision(libRev)
            self.scriptf=self.scriptGen(jobname,outdir,tmpd=tmpd,mem_limit=mem_limit,libRev=libRev,
                                        submode='sub',redirect=redirect)
            cmd=self.scriptf+' &'
            if not fake:
                os.system(cmd)
                sleep(zzz)
                print 'SUB:'+cmd+'\n'
            else:
                print 'FAKE:'+cmd+'\n'
        else:
            print 'COMPLETED '+jobname

            
    def qsub(self,wallt='7.0',tmpd='/tmp/jose',mem_limit='300',
             file_size='',ddisk='100',zzz=0.2,extraflags='',joblist=[],depend=[],libRev='last',
             fake='',redirect=1):
        '''wallt=\7.0\'  walltime, in days (estimated time to complete the job)
           tmpd=\'/tmp/jose\' directory where generic subdir will be created
           zzz=1 time in seconds between consecutive jobs submissions
           ddisk: reserved space under /tmp/jose
           file_size: maximum file size than can be created
           mem_limit: maximum resident memory that can be used
        '''
        jobname=self.name
        if not file_size: file_size=ddisk
        if depend!=[]: extraflags+=' '+add_job_dependencies(depend)
        wallt=formatWallTime(wallt)
        logname=jobname+'.log'
        logf=self.outd+'/'+logname  #;Bye(logf)
        jobid=self.report_status(jobname,joblist=joblist) #;Bye('status='+status)
        PBSoptions=PBSdirectives(locals())
        if not jobid: #job name was not found on the queing system, or it is "C" signaled
            if not self.completed(logf): #there's no logfile with success exit
                libRev=resolveLastRevision(libRev) #obtain last revision file for code library
                self.scriptf=self.scriptGen(tmpd=tmpd,libRev=libRev,submode='qsub',
                                            redirect=redirect,PBSoptions=PBSoptions)
                cmd='qsub '+extraflags+' '+self.scriptf
                if not fake:
                    jobid=os.popen(cmd).readline().rstrip();
                    sleep(zzz)
                    print 'QSUB:'+jobid+':'+cmd+'\n'
                else:
                    print 'FAKE:'+jobid+':'+cmd+'\n'
                return jobid
            else:
                print 'COMPLETED '+jobname  #;Bye(cmd)
                return 'COMPLETED' #return COMPLETED keyword
        else:
            print 'SUBMITTED '+jobname
            return jobid
        
    def store(self,database=codedir+'/python/db/job'):
        import tempfile
        from time import asctime,time,localtime
        from pickle import dump
        self.time=asctime(localtime(time()))
        tempfile.tempdir=database ;#print database
        self.file=tempfile.mktemp() ; #print self.file #;sys.exit()
        dump(self,open(self.file,'w'))
        return database+'/'+self.file

    def help(self):
        return os.popen(self.exed+'/'+self.exe+' -h').readlines()
        
if __name__=='__main__':
    pass

