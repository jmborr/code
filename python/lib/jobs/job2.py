#!/usr/bin/python
import os,re,sys
from utilities.codedir import codedir #we have write permissions to codedir
from utilities.small_utilities import chomp #remove '\n' of strings and lists
from time import sleep

#generate a list of my jobnames in the cluster
def genJobList():
    return chomp(os.popen('qstat|grep jose|tr -s \' \'|cut -d\' \' -f2').readlines())
    
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
    
def pastry(command):
    if os.system(command):#unsuccessful command
        sys.stderr.write('ERROR while executing: '+command+'\n')
        os.system('sync')
        sys.exit(1) #exit the calling (child) process
    return '0'

class job:
    "class to send jobs to LFS system"
    fields=('name','exe','args','exed','time','file','story')
    #name: name of the job
    #exe: name of the executable file
    #args: all arguments (with flags) passed to the executable
    #exed: directory where exe is located
    #time: some time stamp for the job
    #file: filename if we want to save the job as bytecode
    #story: a one-liner explanation of the job
    def __init__(self,name='',exe='',args='',exed='',time='',
                 file='',story=''):
        "optional args self,name='',exe='',exed='',args='',story='None'"
        self.name=name
        self.exe=exe
        self.exed=exed
        self.args=args
        self.story=story
        if not self.name : #executable's name, without extensions
            self.name=re.compile('(\w+)[\.\w*]*$').match(self.exe).group(1)

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

    def report_status(self,submode='all',jobn='',joblist=''):
        if not jobn: jobn=self.subname()
        if joblist: #we passed the list of jobnames present in the cluster 
            if jobn in joblist: return 'from qstat: '+jobn+' present'
            else: return ''
        if submode=='qsub' or submode=='all':
            #expand translates TAB to single space
            #tr -s ' ' will squeeze many spaces into a single one
            reply=chomp(os.popen('qstat|grep '+jobn+'|expand|tr -s \' \'|cut -d\' \' -f 2,5').readline())
            #print 'from qstat: reply=',reply
            if reply:
                return 'from qstat: '+reply
        if submode=='bsub' or submode=='all':
            output=os.popen('bjobs -J '+jobn+' 2>/dev/null').readlines();
            #print 'from bjobs: output=',output
            if output:
                return 'from bsub:',output[1].split()[2] #either RUN, SUSP
        return ''
        
    def completed(self,logfilename=''):
        if not logfilename: logfilename=self.logfile()
        if os.path.exists(logfilename):
            #print 'logfilename=',logfilename
            if not os.system('tail -1 '+logfilename+
                             '|grep "exit mode = 0" > /dev/null'):
                #print 'job is already completed !'
                return 'yes'
        return ''

    def sub(self,tmpd='/tmp/jose',zzz=9,extraflags=''):        
        "optional args: tmpd='/tmp/jose,zzz=9'"
        if not os.path.exists(self.exed):
            print 'Error: '+self.exed+' does not exists!'; return None
        if not os.path.exists(self.exed+'/'+self.exe):
            print 'Error: '+self.exed+'/'+self.exe+' does not exists!';
            return None
        #print "HELLO\n"
        command='./'+self.exe+' '+self.args
        #we escape every " because we submit with bash -c "..", hence the
        #trouble with bash -c "  exec \"args\" "
        command=re.compile('"').sub('\\"',command)
        #cmd will create temporary working directory, copy the executable to
        #it, run it, then remove the temporary directory
        subdir=re.compile('\W').sub('',self.name)
        cmd='source ~/.bash_profile  >/dev/null; '
        cmd=cmd + 'wd='+tmpd+'/'+subdir+'_$$;/bin/mkdir -p $wd;'
        cmd=cmd + '/bin/cp '+self.exed+'/'+self.exe+' $wd;'
        cmd=cmd + 'bash -c "cd $wd;'
        cmd=cmd +command
        cmd=cmd + ';wait; /bin/rm -rf $wd; '
        cmd=cmd +' " &'
        os.system(cmd)
        #print 'SUBMITTED: ',cmd
        sleep(zzz)
        return cmd
    
    def qsub(self,queue='batch',tmpd='/tmp/jose',mem_limit='300',file_size='500',zzz=9,
             jobname='',logname='',outdir='',outputmove='',extraflags='',joblist=''):
        "optional args: queue='batch', tmpd='/tmp/jose',zzz=3"
        if not os.path.exists(self.exed):
            print 'Error: '+self.exed+' does not exists!'; return None
        if not os.path.exists(self.exed+'/'+self.exe):
            print 'Error: '+self.exed+'/'+self.exe+' does not exists!';
            return None
        #logf contains absolute path to the remote log file
        if outdir and logname:
            logf=outdir+'/'+logname
        else:
            logf=self.logfile()
            outdir=os.path.dirname(logf)  #; print 'outdir=',outdir,'logf=',logf
        if not jobname:
            jobname=self.subname(logfilename=logf)  #; print jobname
        if len(jobname) > 15:
            sys.stderr.write('Error: job name '+jobname+' cannot exceed 15 characters')
            return None
        status=self.report_status('all',jobn=jobname,joblist=joblist) ; #print 'status=',status
        if not status: #job name was not found on the queing system
            if not self.completed(logfilename=logf):
                subdir=re.compile('\W').sub('',self.name)
                wd=tmpd+'/'+subdir+'_$$'
                local_log=wd+'/'+os.path.basename(logf)
                command='./'+self.exe+' '+self.args
                cmd=      'if '
                cmd=cmd + 'source ~/.bash_profile >/dev/null && '#set-up envir.
                cmd=cmd + 'wd='+wd+' && '
                cmd=cmd + 'export LOCAL_LOG="'+local_log+'" && '
                cmd=cmd + 'export REMOTE_LOG="'+logf+'" && '
                cmd=cmd + '/bin/rm -f $REMOTE_LOG && '
                cmd=cmd + '/bin/mkdir -p $wd && '
                cmd=cmd + '/bin/cp '+self.exed+'/'+self.exe+' $wd && '
                cmd=cmd + 'cd $wd && '
                cmd=cmd + command
                cmd=cmd + '\\nthen\\nmode=0\\nelse\\nmode=1\\nfi\\n'
                cmd=cmd + 'basename $LOCAL_LOG && '
                cmd=cmd + 'echo "exit mode = $mode" && '              
                cmd=cmd + 'mkdir -p `dirname $REMOTE_LOG`\\n'
                if outputmove:
                    cmd=cmd+'if [ "$mode" -eq 0 ]\\nthen\\n'
                    cmd=cmd+'/bin/mv -f '+outputmove+' '+outdir+'\\n'
                    cmd=cmd+'fi\\n'
                cmd=cmd + '/bin/rm -rf $wd && ' #remove temporary work dir
                cmd=cmd + 'exit $mode'
                cmd=' echo -e \''+cmd+'\' | qsub -r n -c n -m n -l mem='+`mem_limit`+'mb,file='+file_size+'mb -j oe -o '+logf+' -N '+jobname+' -q '+queue+' '+extraflags
                                
                jobid=os.popen(cmd).readline().rstrip(); #fist line is jobid. rstrip removes '\n'
                print 'QSUB:',cmd
                sleep(zzz)
                return jobid
            else:
                cmd = 'COMPLETED '+self.subname()                
                print cmd
                return None
        else:
            #print status
            cmd = status
            print cmd
            return None
    
    def bsub(self,queue='standard',tmpd='/tmp/jose',zzz=2,extraflags='',
             mem_limit=300):
        "optional args: queue='standard', tmpd='/tmp/jose',zzz=2"
        if re.compile('str').search(`type(mem_limit)`):
            mem_limit=int(mem_limit) #cast to int if we passed a string
        if not os.path.exists(self.exed):
            print 'Error: '+self.exed+' does not exists!'; return None
        if not os.path.exists(self.exed+'/'+self.exe):
            print 'Error: '+self.exed+'/'+self.exe+' does not exists!';
            return None
        logf=self.logfile()  #absolute path to the remote log file
        jobname=self.subname(logfilename=logf) #;print jobname
        status=self.report_status('all',jobn=jobname) #; print 'status=',status
        if not status: #job not in queueing server
            if not self.completed(logfilename=logf):
                subdir=re.compile('\W').sub('',self.name)
                wd=tmpd+'/'+subdir+'_$$'
                local_log=wd+'/'+os.path.basename(logf)
                command='./'+self.exe+' '+self.args
                cmd='source ~/.bash_profile && '#set-up environment variables
                cmd=cmd + 'wd='+wd+' && '
                cmd=cmd + 'export LOCAL_LOG="'+local_log+'" && '
                cmd=cmd + 'export REMOTE_LOG="'+logf+'" && '                
                cmd=cmd + '/bin/mkdir -p $wd && '
                cmd=cmd + '/bin/cp '+self.exed+'/'+self.exe+' $wd && '
                cmd=cmd + 'cd $wd && '
                cmd=cmd + command + ' && '
                cmd=cmd + 'echo "exit mode = 0" >> $LOCAL_LOG && '
                cmd=cmd + 'mkdir -p `dirname $REMOTE_LOG` && '
                cmd=cmd + '/bin/mv -f $LOCAL_LOG $REMOTE_LOG && '
                cmd=cmd + '/bin/rm -rf $wd && '
                cmd=cmd + 'exit 0 ; '
                cmd=cmd + '/bin/rm -rf $wd && '
                cmd=cmd + 'exit 1  '

                cmd='bsub -R "rusage[mem='+`mem_limit`+']" -o '+local_log+\
                     ' -e '+local_log+' -J '+jobname+' -q '+queue+' '+\
                     extraflags+' \''+cmd+'\''
                os.system(cmd);
                print 'BSUB:',cmd
                sleep(zzz)
            else:
                print jobname,'completed'
        else:
            print status,jobname

    def rsub(self,node='',tmpd='/tmp/jose',zzz=5,extraflags=''):
        "optional args:tmpd='/tmp/jose',zzz=5,node=get_free_node()"
        #for some unknown reason, writing the following line as a default
        #argument causes get_a_free_node()
        #node=get_a_free_node()
        #if not node:
        #    stderr.write("No free nodes in the cluster. I didn't send the job!\n")
        #    return None
        if not os.path.exists(self.exed):
            print 'Error: '+self.exed+' does not exists!'; return None
        if not os.path.exists(self.exed+'/'+self.exe):
            print 'Error: '+self.exed+'/'+self.exe+' does not exists!';
            return None
        subdir=re.compile('\W').sub('',self.name)
        command='./'+self.exe+' '+self.args
        cmd='source ~/.bash_profile ; '
        cmd=cmd + 'wd='+tmpd+'/'+subdir+'_$RANDOM;/bin/mkdir -p $wd;'
        cmd=cmd + '/bin/cp '+self.exed+'/'+self.exe+' $wd;cd $wd;'
        cmd=cmd + command + ';wait; '
        cmd=cmd + '/bin/rm -rf $wd'
        #send to background, -n flag is neccessary for job to run
        cmd='rsh -n '+node+' \'' + cmd + '\' &'
        #print cmd
        os.system(cmd)

        sleep(zzz); print self.name+' sent to '+node; sleep(0.5)
        return node
        
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
