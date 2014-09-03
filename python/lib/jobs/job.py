#!/usr/bin/python
import os,re,sys
from utilities.codedir import codedir,scratchdir,nodecodedir #we have write permissions to codedir
from utilities.small_utilities import chomp,Bye,todayDate #remove '\n' of strings and lists
from utilities.returnVal import returnVal
from time import sleep
from updateLib.updateLib import resolveLastRevision

sharedexe=nodecodedir+'/python/lib/jobs/shared.py' #script to install/remove temporary libraries
updnodelibexe=nodecodedir+'/python/updateLib/updateLib.py' #script to update nodecodedir
updnodelibexeII='/gpfs1/active/jose/nodelibmirror/code/python/updateLib/updateLib.py'

#add job dependencies. Assume "afterok" type only
def add_job_dependencies(depend):
    if isinstance(depend,list):
        d='-W depend='
        for jobid in depend: d+='afterok:'+jobid+','
        return d[:-1] #remove trailing comma
    return ''

#assume wallt is a string with the number of days the job will run
def formatWallTime(wallt):
    d=float(wallt)
    h=24*d ; m=(h-int(h))*60 ; s=(m-int(m))*60
    h=int(h) ; m=int(m) ; s=int(s)
    return '%d:%02d:%02d'%(h,m,s)
    
    
#add this directory in related search_in_dirs.list file
def add_dir_for_backup(dir,depth=1):
    projectD,subD=re.compile('jose/([\w,\.]+)/*([\w,\./]*)').search(dir).groups()
    if not subD: subD='.'
    searchf=scratchdir+'/'+projectD+'/search_in_dirs.list'
    #check if the directory is already marked for backup
    lines=chomp(os.popen('grep '+subD+' '+searchf).readlines())
    if lines: #find the maximal depth scheduled for search of recent files
        maxdepth=1
        for line in lines:
            depth2=int(line.split()[1])
            if depth2>maxdepth:maxdepth=depth2
        if maxdepth>depth: return None #nothing to do
    pout=open(searchf,'a+').write(subD+' '+`depth`+'\n')    
        
#generate a list of my jobnames in the cluster. Iterate if qstat call fails
def genJobList():
    nit=9 #extra number of qstat invocations
    pt=os.popen('qstat|grep jose|expand|tr -s \' \'|cut -d\' \' -f 1,2,5')
    reply=chomp(pt.readlines())
    while nit and pt.close() : #invoque qstat more if unsuccessful
        sys.stderr.write('qstat ERROR, '+`nit`+' times left to invoque again\n')
        pt=os.popen('qstat|grep jose|expand|tr -s \' \'|cut -d\' \' -f 1,2,5')
        reply=chomp(pt.readlines())
        nit-=1
    if not nit: Bye('ERROR could not qstat successfully for '+`nit+1`+' times')
    return reply

    
#return the name of a log file, and handle special cases
def gen_lf(logrequest,outd,identifier='job'):    
    if logrequest=='_ASK_' or logrequest=='_LOG_' or not logrequest:
        lf=identifier+'.log'
        if logrequest=='_ASK_':
            os.system('echo '+outd+'/'+lf)
            sys.exit(0) #signal job should stop
        else:
            lf=os.getcwd()+'/'+lf #local to current working directory
    else: stderr.write('Error: logrequest either of \'_ASK_\' or \'_LOG_\'')
    return lf
    
def abort(message,exit_code=1):
    sys.stderr.write(message+'\n')
    sys.exit(exit_code)
    
def pastry(command):
    if os.system(command):#unsuccessful command
        sys.stderr.write('ERROR while executing: '+command+'\n')
        sys.exit(1) #exit the calling (child) process
    return '0'

class job:
    """
    class to send jobs to queing system
        
    name: name of the job
    exe: name of the executable file
    args: all arguments (with flags) passed to the executable
    exed: directory where exe is located
    time: some time stamp for the job
    file: filename if we want to save the job as bytecode
    story: a one-liner explanation of the job
    outcmd: /bin/mv line to move output to output directory
    inputs: list of absolute paths to input files/executables
    """
    fields=('name','exe','args','shared','exed','time','file','story','outcmd')
    def __init__(self,name='',exe='',args='',shared='',exed='',time='',
                 file='',story='',outcmd='',inputs=[],incmd=''):
        "optional args self,name='',exe='',exed='',args='',shared='',story=''"
        self.name=name
        self.exe=exe
        self.exed=exed
        self.args=args
        self.shared=shared
        self.story=story
        self.outcmd=outcmd
        self.error=''
        self.scriptf=''
        self.qsubrv=None
        self.incmd=incmd
        self.inputs=inputs  #list of input files to be moved to the working directory
        notfound=self.isinputs()
        if notfound: self.error+='Error: '+','.join(notfound)+' not found\n'
        if not self.name : #executable's name, without extensions
            self.name=re.compile('(\w+)[\.\w*]*$').match(self.exe).group(1)
        if not os.path.exists(self.exed):
            self.error+='Error: '+self.exed+'/'+self.exe+' does not exists!\n'
        if not os.path.exists(self.exed+'/'+self.exe):
            self.error+='Error: '+self.exed+'/'+self.exe+' does not exists!';
            

    def isinputs(self):

        """check input files indeed do exist"""
        
        notfound=[]
        isthere=os.path.exists
        if self.inputs:
            for input in self.inputs:
                if not isthere(input): notfound.append(input)
        return notfound

    
    def __repr__(self):
        return "%(name)s instantiated on %(time)s" % vars(self)
    
    def info(self):
        line=''
        for field in self.fields:
            try: #it may be that instantiated object lacks a particular field
                attr=getattr(self,field)
            except AttributeError:
                continue #jump to next for iteration
            line=line+field.upper()+' = '+attr+'\n'
        return line

    #running the executable with -y _ASK_ should output to STDOUT the
    #absolute path to the logfile associated with this job. From this, we
    #can produce a job name for queueing purposes (subname)    
    def logfile(self):
        return chomp(os.popen(self.exed+'/'+self.exe+' '+self.args+' -y _ASK_').readline())

    def subname(self,logfilename=''):
        if not logfilename: logfilename=self.logfile()
        return os.path.basename(logfilename).replace('.log','')

    def report_status(self,jobn,joblist=[]):
        if not joblist: joblist=genJobList()
        for joblisted in joblist:
            if jobn==joblisted.split()[1]:
                jobid,name,status=joblisted.split()
                if status!='C': return jobid
        return ''
        
    def completed(self,logfilename):
        if os.path.exists(logfilename):
            if not os.system('tail -1 '+logfilename+'|grep "exit mode = 0" > /dev/null'):
                return 'yes'
        return ''

    def temporaryCommands(self):
        """execute some commands for a while"""
        cmd='#temporary commands\n'
        cmd+='/bin/mv '+nodecodedir+'/python/tasser/utilities.py '+nodecodedir+'/python/tasser/tasser_utilities.py\n'#200700824
        cmd+='/bin/rm '+nodecodedir+'/python/tasser/utilities.pyc' #200700824
        return cmd
        
    def updateNodeLib(self,libRev):
        """Run updateLib.py script to update node library. We pass a revision file"""
        foo='updateNodeLib () {\n'
        foo+='#update /library/jose\n'
        foo+=' e=1\n' #initialize exit status as error
        #transient_cmd is a command we want executed for some time
        if libRev:
            foo+=' if [ -e '+updnodelibexe+' ]; then\n'
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
        
    #bring executable sharedexe and run the executable that checks the existence
    #of shared libraries on the node
    def shared_temporal_libraries(self):        
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
        

    #remove working directory, and maybe imported libraries
    def cleanup_exit(self,submode='qsub'):
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

    def PBSdirectives(self,directives):
        """include neccesary directives in script"""
        buf='#directives for PBS\n'
        buf='#PBS -d /tmp\n'
        #PBS -l other=library,mem=839mb,ddisk=500mb,walltime=12:31:00
        #PBS -o /dev/null
        return buf
    
    #generate a script to be used with qsub or sub 
    def scriptGen(self,jobname,outdir,tmpd='/tmp/jose',mem_limit='300',libRev='last',
                  submode='qsub',redirect=1,PBSoptions=''):
        '''

        wallt=\7.0\' walltime, in days (estimated time to complete the
        job) tmpd=\'/tmp/jose\' directory where generic subdir will be
        created zzz=1 time in seconds between consecutive jobs
        submissions

        '''
#        qsubdir='/tmp/jose/qsub/'+todayDate() #subdirectory to deposit the script
        qsubdir=scratchdir+'/qsub/'+todayDate() #subdirectory to deposit the script
        if not os.path.exists(qsubdir):  pastry('/bin/mkdir -p '+qsubdir)
        script=qsubdir+'/qsub.sh_'+jobname #full script file name

        if len(jobname) > 15:
            sys.stderr.write('Error: job name '+jobname+' cannot exceed 15 characters')
            return ''
        if not os.path.exists(outdir): os.system('/bin/mkdir -p '+outdir)
        buf=''
        ulimit=int(float(mem_limit)*1024) #maximum resident memory size (Kb) to prevent swapping
        wd=tmpd+'/'+ re.compile('\W').sub('',self.name) +'_$$' #working directory
        logname=jobname+'.log'
        local_log=wd+'/'+logname
        remote_log=outdir+'/'+logname
        command='./'+self.exe+' '+self.args
        buf=  '#!/bin/bash\n\n'        
        buf+= self.PBSdirectives(PBSoptions)+'\n\n'
        buf+= self.updateNodeLib(libRev)+'\n\n' #bash function to update library
        buf+= self.shared_temporal_libraries()+'\n\n' #bash function to bring temporary libs
        buf+= self.cleanup_exit(submode=submode)+'\n\n' #bash function to clean exit
        buf+= 'echo "'+script+'"\n' #write script name withing script body
        buf+= 'hostname\n'          #node where job will be run
        buf+= 'echo $PBS_JOBID\n'
        buf+= 'ulimit -m '+`ulimit`+' #maximum resident memory\n'
        buf+= 'source ~/.bash_profile >/dev/null #environment variables\n'
        buf+= 'wd='+wd+' #working directory\n'
        buf+= '/bin/mkdir -p $wd\n'
        buf+= 'export LOCAL_LOG="'+local_log+'"\n'
        buf+= '/bin/touch $LOCAL_LOG\n'
        if submode=='sub' and redirect:
            buf+='exec &> $LOCAL_LOG\n' #redirect STODOUT, STDERR to LOCAL_LOG
        buf+= 'export REMOTE_LOG="'+remote_log+'"\n'

        buf+= 'if [ -f $REMOTE_LOG ]; then\n'
        buf+= ' /bin/rm -f $REMOTE_LOG\n'
        buf+= 'fi\n\n'

        buf+= 'trap "cleanup_exit 1" TERM\n\n'        #in case of killing job

        buf+= self.temporaryCommands()+'\n\n' #some temporary commands

        buf+= '#update node code library && copy necessary libraries to /tmp/jose\n'
        buf+= 'if !('
        buf+= 'updateNodeLib && ' #update node code library
        buf+= 'shared_temporal_libraries _PREPARE_' #copy necessary libraries to tmpd
        buf+= ');then\n'
        buf+= '  cleanup_exit 1\n'
        buf+= 'fi\n\n'

        buf+= 'cd $wd\n\n'
        buf+= '#bring all necessary inputs and perform input commands AT THE END\n'
        cmd='/bin/cp '+self.exed+'/'+self.exe+' '+' '.join(self.inputs)+' $wd'
        if self.incmd: cmd=cmd+' && '+self.incmd
        buf+= 'if !('+cmd+');then\n'
        buf+= '  cleanup_exit 1\n'
        buf+= 'fi\n'
        
        buf+= 'if !('  #we need to test command aside previous tests
        buf+= command
        buf+= ');then\n'
        buf+= '  cleanup_exit 1\n'
        buf+= 'fi\n\n'

        buf+= self.outcmd+'\n' #move output files
        
        buf+= 'cleanup_exit 0'

        open(script,'w').write(buf)
        pastry('chmod u+x '+script)

        return script

    ##################################################################################
    
    def dormant(self,ibatch):

        """check for qsub error and take a nap if we submitted many consecutive jobs"""
        from art.ascii import nap
        nbatch=100      #number of jobs to be sumitted without extended sleeping
        hibernation=240 #4 minutes nap after submitting nbatch jobs
        
        if self.qsubrv.st:
            sys.stderr.write(self.error+'\n')
            return ibatch
        if self.qsubrv.msg[0:4]=='QSUB': ibatch+=1
        if ibatch==nbatch:
            sys.stdout.write(nap+' taking a nap for '+`int(hibernation/60)`+' minutes\n')
            sleep(hibernation)
            ibatch=0
        return ibatch

    ##################################################################################
    
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
            
    ##################################################################################

    def qsub(self,jobname,outdir,wallt='7.0',tmpd='/tmp/jose',mem_limit='300',
             file_size='',ddisk='200',zzz=0.6,extraflags='',joblist=[],depend=[],libRev='last',
             fake='',redirect=1):
        '''wallt=\7.0\'  walltime, in days (estimated time to complete the job)
           tmpd=\'/tmp/jose\' directory where generic subdir will be created
           zzz=3 time in seconds between consecutive jobs submissions
           ddisk: reserved space under /tmp/jose
           file_size: maximum file size than can be created
           mem_limit: maximum resident memory that can be used
        '''
        #print 'shared=',self.shared  ;  sys.exit(1)
        #print self.args ; sys.exit(1)

        if self.error: #some initialization error
            self.qsubrv=returnVal(st=1,errMsg=self.error)
            return self.qsubrv
        if not file_size: file_size=ddisk
        if depend!=[]: extraflags+=' '+add_job_dependencies(depend)
        wallt=formatWallTime(wallt)
        logname=jobname+'.log'
        logf=outdir+'/'+logname  #;Bye(logf)
        jobid=self.report_status(jobname,joblist=joblist) #;Bye('status='+status)
        PBSoptions='' #for the future we will put PBS options in the body of the script
        retVal=returnVal()
        if not jobid: #job name was not found on the queing system, or it is "C" signaled
            if not self.completed(logf):
                libRev=resolveLastRevision(libRev) #obtain last revision file for code library
                self.scriptf=self.scriptGen(jobname,outdir,tmpd=tmpd,mem_limit=mem_limit,
                                           libRev=libRev, submode='qsub',redirect=redirect,
                                            PBSoptions=PBSoptions)
                cmd='qsub -r n -c n -m n -l mem='+mem_limit+'mb,file='+file_size+'mb,ddisk='+ddisk+'mb,walltime='+wallt+' -j oe -o '+logf+' -N '+jobname+' '+extraflags+' '+self.scriptf
                if not fake:
                    pt=os.popen(cmd)
                    jobid=pt.readline().rstrip()
                    #os.system('/bin/rm '+self.scriptf) #script already been loaded to memory
                    sleep(zzz)
                    st=pt.close()
                    emsg=''
                    if st: emsg='QSUB ERROR '+`st`
                    self.qsubrv=returnVal(st=st,errMsg=emsg,msg='QSUB:'+jobid+':'+cmd,val=jobid)
                else:
                    self.qsubrv=returnVal(st=1,errMsg='FAKE',msg='FAKE:'+jobid+':'+cmd,val=jobid)
            else:
                self.qsubrv=returnVal(st=0,msg='COMPLETED '+jobname,val=jobid)
        else:
             self.qsubrv=returnVal(st=0,msg='SUBMITTED '+jobname,val=jobid)
        self.error+=retVal.errMsg
        print self.qsubrv.msg
        return self.qsubrv
    
    ##################################################################################

    def bsub(self,queue='standard',tmpd='/tmp/jose',zzz=2,extraflags='',
             mem_limit=300):
        pass #obsolete method
    
    def rsub(self,node='',tmpd='/tmp/jose',zzz=5,extraflags=''):
        "optional args:tmpd='/tmp/jose',zzz=5,node=get_free_node()"
        pass #obsolete method
    
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
    j=job(exe='mkseq.py', exed='/gpfs1/active/jose/code/python/tasser/tasser1.0/mkinput', args='-i 101m_ -o junkdir' )
    j.store()
