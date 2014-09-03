#!/usr/bin/python
import os,re,sys
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.small_utilities import chomp
from jobs.job import pastry
from fcntl import flock,LOCK_EX

"""Install or remove temporal libraries under /tmp/jose
   Example: /tmp/jose/pdb:pdb.tbz2,/tmp/jose/pdbx:pdbx.tbz2
   If we expand pdb.tbz2 on /tmp/jose, then we will generate
   /tmp/jose/pdb database
   /tmp/jose is the root directory
"""

#global variables
lockfiled='/tmp/jose/lockfiles' #store info on which processes are using which library
pattern={'pid':re.compile('pid=(\d+)')}
#********************************************************
def print_header():
    """wellcoming message"""
    buf= '\n'
    buf+='          **************************\n'
    buf+='          *** entering shared.py ***\n'
    sys.stdout.write(buf)
    sys.stdout.flush()
#********************************************************
def print_bye():
    """exiting program message"""
    buf= '          *** exiting (shared.py) ***\n'
    buf+='          ***************************\n\n'
    sys.stdout.write(buf)
    sys.stdout.flush()
#********************************************************
def does_pid_exists(pid):
    """check if a PID does exists"""
    return len( os.popen('ps -fu jose | grep -P "jose\s+'+pid+'\s+"').readlines() )
#********************************************************
def get_dict(listlib):
    dict={}
    for pair in listlib.split(','):
        dir,tbz2=pair.split(':')
        tbz2=addAbsPath(tbz2) #expands codedir,scratchdir, and the like
        dict[dir]=tbz2
    return dict
#********************************************************
#********************************************************
def prepare_libs(listlib,jobid):
    """bring libraries under /tmp/jose
    The libraries may already be in place and other processes are
    using them"""    
    print_header()
    dict=get_dict(listlib) #get the dir:tbz2 pairs in a dictionary
    currdir=os.getcwd()
    for dir in dict.keys():
        tbz2=dict[dir]
        root=chomp(os.popen('/usr/bin/dirname '+dir).readline()) #up once in the directory tree
        if not os.path.exists(root): pastry('mkdir -p '+root)    #root dir wasn't not present
        suffix=tbz2.replace('/','.')

        busyfile=os.path.join(lockfiled,'busy.'+suffix) #list of processes using the library
        busyp=open(busyfile,'a+') #file descriptor for update. Will create if not present
        flock(busyp.fileno(),LOCK_EX) #ask for exclusive lock, otherwise wait until lock released

        #update busyfile contents and unpack library if necessary
        busylines=busyp.readlines()
        nbusylines=len(busylines)
        if busylines: #check integrity of existing library and running processes
            if not os.path.exists(dir): #library does not exists but there are listed processes
                sys.stderr.write('ERROR in shared.prepare_libs: library "'+dir+'" does not exists but busyfile "'+busyfile+'" is not empty. We are cleaning the busyfile and unpacking the library\n')
                busylines=[] #mark as no processes needing the library
            else: #library exists, but do listed proccesses exists?
                for busyline in busylines:
                    pid=pattern['pid'].search(busyline).group(1) #retrieve PID
                    if not does_pid_exists(pid): #pid was not properly removed
                        busylines.remove(busyline)
                if len(busylines) != nbusylines: #some of the listed processes were garbage
                    sys.stderr.write('ERROR in shared.prepare_libs: some of the listed PID in the busyfile did not exists anymore\n')
                    busyp.truncate(0)
                    busyp.write(''.join(busylines)) #update busyfile with only running processes
                else:
                    sys.stdout.write('Existing library '+tbz2+' in '+dir+'\n')
        if not busylines: #this is signal to unpack library
            sys.stdout.write('shared.prepare_libs unpacking '+tbz2+' to '+dir+'\n')
            busyp.truncate(0) #remove garbage lines, if any
            pastry('tar jxf '+tbz2+' -C '+root)                    
        busyp.write('pid='+jobid+'\n') #add proccess ID to the list
        busyp.close() #release the lock
    print_bye()
    return 0 #marks success in UNIX style
#********************************************************   
#********************************************************   
def cleanup_libs(listlib,jobid):
    """Attempt to remove passes libraries in /tmp/jose
    Libraries will not be removed if other running processes are doing
    use of them"""    
    print_header()
    print 'cleanup_libs entered: '+jobid+' cleaning libraries '+listlib
    dict=get_dict(listlib)
    for dir in dict.keys():
        tbz2=dict[dir]
        root=chomp(os.popen('/usr/bin/dirname '+dir).readline())
        suffix=tbz2.replace('/','.')

        busyfile=os.path.join(lockfiled,'busy.'+suffix) #list of processes using the library
        busyp=open(busyfile,'a+') #file descriptor for update. Will create if not present
        flock(busyp.fileno(),LOCK_EX) #ask for exclusive lock, otherwise wait until lock released

        #remove jobid from list of processes
        busylines=busyp.readlines()
        busyp.truncate(0) #remove contents
        for busyline in busylines: #search and remove jobid from the list
            if jobid in busyline: busylines.remove(busyline)
        if busylines: #there are remaining processes
            busyp.writelines(busylines)
        else: #no remaining processes, thus remove library
            pastry('/bin/rm -rf '+dir)
        busyp.close() #release lock
    print_bye()
    return 0
#********************************************************
#********************************************************
if __name__=='__main__':
    inpHand('Usage:shared.py [options]',
            ' -a _R_listlib libline (example:dir1:tbz21,dir2:tbz22)',
            ' -b _R_jobid jobid',
            ' -c _R_task (either _PREPARE_ or _CLEANUP_)'
            ).parse(locals(),sys.argv)

    exit_status=0
    if   task=='_PREPARE_' : exit_status=prepare_libs(listlib,jobid)
    elif task=='_CLEANUP_' : exit_status=cleanup_libs(listlib,jobid)
    else: exit_status=1
    sys.exit( exit_status )
