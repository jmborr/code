#!/usr/bin/python

import os,sys,pickle,copy,re
from time import sleep
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,Bye
from fcntl import flock,LOCK_EX

nodeLib='/library/jose'
noderevdir=nodeLib+'/code/python/updateLib/revisions'
mirror='/gpfs1/active/jose/nodelibmirror'
revdir=mirror+'/code/python/updateLib/revisions'
lockfiled='/tmp/jose/lockfiles'
busyfile=os.path.join(lockfiled,'busy.updateLib')
updateExe='/gpfs1/active/jose/code/python/updateLib/updateLib.py'

rsyncExitKeys=[0,1,2,3,4,5,10,11,12,13,14,20,21,22,23,24,30]
rsyncExitValues={0:'Success',1:'Syntax or usage error',2:'Protocol incompatibility',3:'Errors selecting input/output files, dirs',4:'Requested action not supported: an attempt was made to  manipulate 64-bit files on a platform that cannot support them; or an option was specified that is supported by the client and not by the server',5:'Error starting client-server protocol',10:'Error in socket I/O',11:'Error in file I/O',12:'Error in rsync protocol data stream',13:'Errors with program diagnostics',14:'Error in IPC code',20:'Received SIGUSR1 or SIGINT',21:'Some error returned by waitpid()',22:'Error allocating core memory buffers',23:'Partial transfer due to error',24:'Partial transfer due to vanished source files',30:'Timeout in data send/receive',999:'Unidentified rsync exit value'}
##########################################################################
def print_header():
    """wellcoming message"""
    buf= '\n'
    buf+='          ****************************\n'
    buf+='          *** entring updateLib.py ***\n'
    buf+='          ****************************\n'
    sys.stdout.write(buf)
    sys.stdout.flush()
##########################################################################
def print_bye():
    """exiting program message"""
    buf= '          *** exiting updateLib.py ***\n'
    buf+='          ****************************\n'
    sys.stdout.write(buf)
    sys.stdout.flush()
##########################################################################
def resolveLastRevision(libRev):
    '''
    find last revision, pointed to by "last" symbolic link within revdir directory
    '''
    if libRev=='last': #this job requires resolving "last"
        return os.popen('ls -l '+revdir+'/last').readline().strip().split('/')[-1]
    return libRev

#remove duplicates from a list
def unique(lst):
    d = {}
    for item in lst: d[item] = None
    return d.keys()

def remove_comment_lines(files):
    culled=[]
    for file in files:
        if file[0]!='#': culled.append(file)
    return unique(culled)

#will gather all files for revisions newer than olldate revision
def parse_files(files,olddate='last'):
    if olddate=='last': return remove_comment_lines(files)
    olddate=int(olddate)
    culled=[]
    if not files:
        sys.stderr.write('ERROR in updateLib.parse_files: no files to parse\n')
        return culled
    for file in files:
        if file[0:2]=='#20': #it's a revision version date
            if int(file[1:9])<olddate: break #we reached a too old revision
        if file[0]!='#': culled.append(file)
        elif file[0:2]=='#/': culled.append(file[1:])
    return unique(culled) #remove duplicate file names
    
#return a revision file name with todays date plus a "library version" of two digits
def genRevName(db):
    libversion='00' #initialize
    date=os.popen('date +%F').readline().strip().replace('-','') #YYYYMMDD
    #check if there are files with same date. list[-1] is latest libversion
    cmd='ls -1 '+db+'/'+date+'* 2>/dev/null'
    list=os.popen(cmd).readlines()    
    if list: libversion='%02d'%(1+int(list[-1].strip()[-2:])) #last two digits of list[-1]
    return date+libversion
    
class revfile:
    def __init__(self,files=[],other={}):
        #pointer to previous (p) and next (n) revisions
        self.point={'p':'','c':'','n':''}
        self.files=files
        #extensions if we need to add other properties later
        self.ext={}
        if other:
            for property in other.keys(): self.ext[property]=other[property]

#    def info(self):
#        all=[]
        
#        return all

#    def __repr__(self):
#        return self.info()
             
    #add revision object to the end of the linked list of revisions
    def insertAtBack(self,files,db=revdir,other={}):
        self.files=files
        self.point['c']=genRevName(db) #give a file name to the revision object
        #load last revision. Create a new linked list if there is no revision        
        link=os.path.join(db,'last')
        if os.path.exists(link):
            last=pickle.load(open(link,'r'))
            #update pointers, since self will become new "last"
            self.point['p']=last.point['c'] #previous of self is last
            last.point['n']=self.point['c'] #next of last is self
            #overwrite last with new info
            fprev=os.path.join(db,last.point['c'])  ;   pickle.dump(last,open(fprev,'w'))
            #store new revision
            flast=os.path.join(db,self.point['c'])  ;   pickle.dump(self,open(flast,'w'))
            #update "current" and "last" symbolic links
            os.system('/bin/rm '+link+' && /bin/ln -s '+flast+' '+link)
            link=os.path.join(db,'current')
            os.system('/bin/rm '+link+' && /bin/ln -s '+flast+' '+link)            
        else: #create a new linked lists of revisions
            flast=os.path.join(db,self.point['c'])
            pickle.dump(self,open(flast,'w'))
            os.system('/bin/ln -s '+flast+' '+link) #symbolic link "last" point to end of list
            link=os.path.join(db,'first')
            os.system('/bin/ln -s '+flast+' '+link) #symbolic link "first" point to start of list
            link=os.path.join(db,'current')
            os.system('/bin/ln -s '+flast+' '+link) #symbolic link "current" point to end of list
        print 'Revision '+self.point['c']+' created'
        return self.point['c']

#copy files from prefix to libdir. Change the root directory from prefix to libdir
def transferFiles(files,oldprefix='/gpfs1/active/jose',newprefix=mirror):
    localfiles=[]
    if not files:
        sys.stderr.write('ERROR in updateLib.transferFiles: no files to transfer\n')
        return localfiles
    pattern=re.compile(oldprefix+'/([\S]+)') #retrieve filename withouth the "prefix" part
    for file in files:
        if not os.path.exists(file):
            sys.stderr.write(file+' does not exists!\n')
            continue
        match=pattern.match(file)
        if match:
            suffix=match.group(1) #filename except for prefix
            fname=os.path.basename(suffix) #name of file
            dir=newprefix+'/'+os.path.dirname(suffix)
            if not os.path.exists(dir): os.system('/bin/mkdir -p '+dir)
            localfile=dir+'/'+fname
            os.system('/bin/cp -r -p '+file+' '+localfile)
            localfiles.append(localfile)
        else: return [] #failure
    return localfiles #success Now the list has libdir as root directory

def transferToNode(filesToUpdate,oldprefix=mirror,newprefix=nodeLib):
    pattern=re.compile(oldprefix+'/([\S]+)') #retrieve filename withouth the "prefix" part
    updatedFiles={ 'success':[], 'failure':[] }
    for file in filesToUpdate:
        match=pattern.match(file)
        if match:
            suffix=match.group(1) #filename except for prefix
            fname=os.path.basename(suffix) #name of file
            dir=newprefix+'/'+os.path.dirname(suffix)
            if not os.path.exists(dir): os.system('/bin/mkdir -p '+dir)
            localfile=dir+'/'+fname
            val=os.system('/bin/cp -r -p '+file+' '+localfile)
            if val: updatedFiles['failure'].append(file) #failed copy
            else:   updatedFiles['success'].append(file) #successful copy
        else: updatedFiles['failure'].append(file) #failed match
    for key in updatedFiles.keys():
        if updatedFiles[key]==[]: updatedFiles[key]=['none',] 
    return updatedFiles

##########################################################################
            
def updateSomeNodes():
    '''
    Run updateLib.py on fec01 to fec04 and b1 to b14
    '''
    #update front nodes and some computing nodes
    updatecmd=' "source ~/.bash_profile && '+updateExe+' -a last" &>/dev/null'
    #update front nodes
    sys.stdout.write('Updating ')
    for i in range(1,5):
        sys.stdout.write(' fec0'+`i`) ;  sys.stdout.flush()
        os.system('/usr/bin/ssh fec0'+`i`+updatecmd)
    sys.stdout.write('\n') ;  sys.stdout.flush()
    #update interactive computing nodes
    sys.stdout.write('Updating ')
    for i in range(1,15):
        sys.stdout.write(' b'+`i`) ;  sys.stdout.flush()
        os.system('/usr/bin/rsh b'+`i`+updatecmd)
    sys.stdout.write('\n') ;  sys.stdout.flush()
##########################################################################
def commentlist(version,listf):
    '''
    We comment the files listed in listf that were updated. We also
    include a header with the new version
    '''
    lines=open(listf,'r').readlines()
    pt=open(listf,'w') #update now listf by commmenting updated files
    pt.write('#'+version+'\n')
    for line in lines:
        if line[0]!='#': line='#'+line
        pt.write(line)
    pt.close()
##########################################################################    
def add_new_revision(files,listf,db=revdir):
    """create a new revision file with a list of files to update.
    Bring the files to /library/jose in fec01, then rsync to fec04."""
    files=transferFiles(files)
    exit_mode=1 #failure (we follow UNIX convention)
    if files: #if files were transferred, then add a new revision file
        version=revfile([]).insertAtBack(files,db=db)
        updateSomeNodes() #update fec01 to fec04 and b1 to b14
        commentlist(version,listf)
        exit_mode=0 #success (we follow UNIX convention)
    else:
        sys.stderr.write('ERRORin updateLib.add_new_revision (no files were transferred)')
    return exit_mode
    
#update the revision file in the library, and bring whatever files are required
def updateNodeLib(rev,noderevdir=noderevdir,revdir=revdir):
    print_header()       
    if not os.path.exists(lockfiled): os.system('/bin/mkdir -p '+lockfiled)

    #open busyfile and put an exclusive lock. If other process is
    #doing the update then it will have a lock in place. In this case,
    #the program will wait until the other update is finished and
    #therefore the lock removed.    
    busyp=open(busyfile,'a+') #file descriptor for update
    flock(busyp.fileno(),LOCK_EX) #ask for exclusive lock, otherwise wait until lock released

    #resolve last revision on remote mirror
    if rev=='last': rev=resolveLastRevision(rev)

    #update busyfile contents
    busyp.truncate(0) 
    busyp.write('pid='+`os.getpid()`+' rev='+rev+'\n') #mark lock with pid and revision

    current=os.path.join(noderevdir,'current')   #simbolic link to local revision file
    if not os.path.exists(current):              #we have not set up the revisions directory
        if not os.path.exists(noderevdir): os.system('/bin/mkdir -p '+noderevdir)
        localDate=0                          #this will ensure we load all remote revisions
    else:
        currNodeRev=pickle.load( open(current,'r') ) #load local revision into object
        localDate=int(currNodeRev.point['c'])        #name of the local revision file
    if localDate>=int(rev): #no need to update
        sys.stdout.write('local library "/tmp/jose/code" is up to date\n')
        busyp.close()     #liberate lock
        print_bye()
        return 0  #nothing to do and exit update as success
    #bring all revision files to local repository. We don't just bring
    #"rev" because we may need revisions that are older than "rev",
    #but newer than "current".
    revs=[] #list of revision files that are used for the update
    cmd='/usr/bin/rsync -az --delete '+revdir+'/ '+noderevdir+'/ 2>/dev/null >/dev/null'
    exitVal=os.system(cmd)
    if exitVal not in rsyncExitKeys: exitVal=999 
    if not exitVal: #rsync transfer was a success
        if rev=='last':
            rev=os.popen('ls -l '+noderevdir+'/last').readline().split('/')[-1].strip()
        targetRev=pickle.load( open(os.path.join(noderevdir,rev),'r') )
        revs.append(targetRev.point['c'])             #initalize list of revisions to update
        filesToUpdate=targetRev.files[:]                        #initialize
        prevRevFile=targetRev.point['p']                        #name of previous revision
        while prevRevFile and localDate < int(prevRevFile):     #earlier revisions are needed
            revs.append(prevRevFile)
            prevRev=pickle.load( open(os.path.join(noderevdir,prevRevFile),'r') )        
            filesToUpdate+=prevRev.files
            prevRevFile=prevRev.point['p']                      #name of previous revision
        filesToUpdate=unique(filesToUpdate)                     #remove list duplicates
        updatedFiles=transferToNode(filesToUpdate)              #update files in the local node
        if os.path.exists(current): os.system('/bin/rm '+current)     
        os.system('/bin/ln -s '+os.path.join(noderevdir,rev)+' '+current) #update link
        buf= 'updating node code library using revisions: '+' '.join(revs)+'\n'
        buf+='Updated files:\n'+'\n  '.join(updatedFiles['success'])+'\n'
        buf+='Failed updated files:\n'+'\n  '.join(updatedFiles['failure'])+'\n'
        sys.stdout.write(buf)
        busyp.close()     #liberate lock
        print_bye()
        return 0 #success (we follow UNIX convention)
    else:
        sys.stderr.write('ERROR from updateLib::updateNodeLib: could not update node code library. rsync reports '+rsyncExitValues[exitVal]+'\n')
    busyp.close()     #liberate lock
    print_bye()
    return 1 #failure (we follow UNIX convention)
        
if __name__=='__main__':
    inpHand('Usage: updateLib.py [options]',
            ' -a __rev revision file name (def: None (will do nothing))',              
            ' -b _A_listf filename with list of files to update. Will create new revision',
            ' -c __olddate Will create a new revision but will add also all files for revisions newer than the YYYYMMDD date, if these revisions are also included in the list of files to update (def:None)',
            ).parse(locals(),sys.argv)

    if not olddate: olddate='last'
    
    if listf: #We want to add a new revision file
        files=chomp(open(listf,'r').readlines())
        files=parse_files(files,olddate=olddate)
        sys.exit( add_new_revision(files,listf) )        
    else:     #We want to update the node library
        if not rev: sys.exit(0) #nothing to do
        sys.exit( updateNodeLib(rev) )
