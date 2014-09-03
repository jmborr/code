#!/usr/bin/python
import os,re,sys
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.small_utilities import chomp #remove '\n' of strings and lists
from jobs.job import pastry,abort
from time import sleep
from random import randint

"""Example: /tmp/jose/pdb:pdb.tbz2,/tmp/jose/pdbx:pdbx.tbz2
   If we expand pdb.tbz2 on /tmp/jose, then we will generate
   /tmp/jose/pdb database
   /tmp/jose is the root directory
"""

def get_dict(listlib):
    dict={}
    for pair in listlib.split(','):
        dir,tbz2=pair.split(':')
        tbz2=addAbsPath(tbz2) #expands codedir,scratchdir, and the like
        dict[dir]=tbz2
    return dict

def prepare_libs(listlib,jobid):
    dict=get_dict(listlib) #get the dir:tbz2 pairs in a dictionary
    currdir=os.getcwd()

    for dir in dict.keys():

        sl=`randint(0,59)` ; print jobid,'sleeping',sl
        pastry('sleep '+sl) #avoid several jobs starting at the same time
        
        tbz2=dict[dir]
        root=chomp(os.popen('/usr/bin/dirname '+dir).readline()) #up once in the directory tree
        if not os.path.exists(root): pastry('mkdir -p '+root)    #root dir wasn't not present
        suffix=tbz2.replace('/','.')
        
        #Wait if library is being created by some other job running on the same node.
        #the name of the expandname is "expand" plus the name of the compressed library, but
        #we change '/' by '.'
        expandname='expand'+suffix       ; expandfile=root+'/'+expandname #lib. being created
        n=0
        while os.path.exists(expandfile):
            print jobid,'waiting for library to be expanded'
            sleep(60) #check every minute
            n+=1
            if n>60:
                abort('Error: I\'ve wait for one hour for '+expandfile+' to be removed')

        busyname='busy'+suffix+'.'+jobid ; busyfile  =root+'/'+busyname   #lib. being used
        pastry('touch '+busyfile)   ; print 'touch '+busyfile #signal this job will use the library
        
        #Wait if library is being removed by other job running in same node
        cleanname='clean'+suffix         ; cleanfile =root+'/'+cleanname  #lib. being removed
        n=0
        while os.path.exists(cleanfile):
            sleep(60)  #check every minute
            n+=1
            if n>60: #probably the cleanfile is a remnant of a finished job
                pastry('/bin/rm '+cleanfile)
                
        #Create the database if it doesn't exists
        if not os.path.exists(dir):     #library not present
            os.chdir(root)                ; print jobid,'expanding library from root'
            pastry('touch '+expandfile)
            #print expandfile+' marks that we\'re unzipping the database'
            pastry('tar jxf '+tbz2)       ; print 'unzip the database '+tbz2
            pastry('/bin/rm '+expandfile)
            print 'we finished unzipping the database, remove '+expandfile
            os.chdir(currdir)             #; print 'go back to previous working directory'

def cleanup_libs(listlib,jobid):
    dict=get_dict(listlib)
    for dir in dict.keys():
        tbz2=dict[dir]
        root=chomp(os.popen('/usr/bin/dirname '+dir).readline())
        suffix=tbz2.replace('/','.')
        busyname='busy'+suffix+'.'+jobid ; busyfile  =root+'/'+busyname   #lib. being used
        cleanname='clean'+suffix         ; cleanfile =root+'/'+cleanname  #lib. being removed

        pastry('touch '  +cleanfile)  #; print cleanfile+' marks lib may be in erase proccessing'
        pastry('/bin/rm '+busyfile)  ;  print 'removing busyfile '+busyfile

        #removing the  library
        busyfilelist=chomp(os.popen('ls '+root+'/busy'+suffix+'* 2>/dev/null').readlines())
        if busyfilelist: #files marking the library is being used by other jobs
            #for each busyfile, check if related job is actually running
            for busyfile in busyfilelist:
                jobid=split('.')[1] #PID process of the job running
                ps=os.system('ps '+jobid).readlines()
                if len(ps)==1: pastry('/bin/rm '+busyfile) #there is no job realted to the busyfile
        else:
            pastry('/bin/rm -rf '+dir) ; print 'removing library'+dir
        pastry('/bin/rm '+cleanfile)   ; print 'removing clean file '+cleanfile
        
        
ih=inpHand('Usage:shared.py [options]',
           ' -a _R_listlib libline (example:dir1:tbz21,dir2:tbz22)',
           ' -b _R_jobid jobid',
           ' -c _R_task (either _PREPARE_ or _CLEANUP_)'
           )
ih .parse(locals(),sys.argv)

if   task=='_PREPARE_' : prepare_libs(listlib,jobid)
elif task=='_CLEANUP_' : cleanup_libs(listlib,jobid)
else: sys.exit(1)
sys.exit(0)
