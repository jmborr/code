#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,junkName
from amber.amber9 import ptraj

codedir=os.environ['CODED']
isthere=os.path.exists
joink=os.path.join
basename=os.path.basename
template="""#!/bin/bash
#PBS -A BIP003
#PBS -N _JOBNAME_
#PBS -j oe
#PBS -o _LOGDIR_/${PBS_JOBID}._JOBNAME_.log
#PBS -l walltime=_WALLTIME_,size=_SIZE_

SANDER='aprun -n_SIZE_ /ccs/home/jmborr/code/software/amber9/exe/pmemd'

TYP=_TYP_
PROD=_PROD_

#submit simulation
echo \"current script _SCRIPTNAME_\"
cd $PROD/curr
echo \"$SANDER -O -i in -p top -c rst -ref ref -o out -x crd -v vel -r rst.new &\"
$SANDER -O -i in -p top -c rst -ref ref -o out -x crd -v vel -r rst.new &
pid=$!
echo \"pid=$pid TYP=$TYP $PROD\"
echo

#check if job finished.
wait $pid
echo \"arranging $PROD\"
cd $PROD/curr
/bin/mv out     ../prev/$TYP._ENDFR_.out
/bin/mv crd     ../prev/$TYP._STARTFR___ENDFR_.crd
/bin/mv vel     ../prev/$TYP._STARTFR___ENDFR_.vel
/bin/mv rst.new ../prev/$TYP._ENDFR_.rst
/bin/rm rst
/bin/ln -s      ../prev/$TYP._ENDFR_.rst rst
echo $TYP._ENDFR_.rst >> ../rst.list
qstat -f ${PBS_JOBID}  #output job info
"""
def getwalltime(ncores):
    """return walltime according to jaguar scheduling policies"""
    if ncores<=256:    walltime='1:30:00'
    elif ncores<=512:  walltime='2:30:00'
    elif ncores<=1024: walltime='4:00:00'
    elif ncores<=6000: walltime='12:00:00'
    else: walltime='24:00:00'
    return walltime

def getLastframe(infile,maxtime):
    """find total number of coordinate frames for a given time(ns)"""
    buf=''.join(open(infile,'r').readlines())
    dt=float(re.compile(' dt=(\d+\.\d+)').search(buf).group(1)) #ps/step
    ntwx=int(re.compile('ntwx=(\d+)').search(buf).group(1)) #steps/frame
    lastframe=int(float(maxtime)*1E3/(dt*ntwx))
    print 'lastframe=',lastframe
    return lastframe

def getDeltaFrame(infile):
    """find number of coordinate frames in a single run"""
    buf=''.join(open(infile,'r').readlines())
    nstlim=int(re.compile('nstlim=(\d+)').search(buf).group(1))
    ntwx=int(re.compile('ntwx=(\d+)').search(buf).group(1)) #steps/frame
    deltaframe=int(nstlim/ntwx)
    print 'deltaframe=',deltaframe
    return deltaframe

def setUp(prodd,TYP,jobname,maxtime,ncores):
    """set up the production directory"""
    subds=['curr','prev','sub'] #create subdirectories
    for subd in subds: os.system('/bin/mkdir '+joink(prodd,subd))
    inputs=['in','top','ref','rst'] #create symbolic links
    for input in inputs:         
        file='../init/'+TYP+'.'+input  #file receptor of the link
        currd=joink(prodd,'curr')
        os.system('cd '+currd+' && ln -s '+file+' '+input)
    #obtain lastframe and deltaframe
    infile=joink(prodd,'init',TYP+'.in')
    lastframe=getLastframe(infile,maxtime)
    dfr=getDeltaFrame(infile)
    #generate scripts for qsub
    scriptd=joink(prodd,'sub')
    Nruns=lastframe/dfr #number of scripts
    currfr=0
    for i in range(Nruns):
        startframe='%05d'%(currfr+1)
        endframe='%05d'%(currfr+dfr)
        script=joink(scriptd,startframe+'_'+endframe+'.sub')
        buf=template
        repl={ '_JOBNAME_':jobname,
               '_LOGDIR_':joink(prodd,'curr'),
               '_WALLTIME_':getwalltime(int(ncores)),
               '_SIZE_':ncores,
               '_TYP_' :TYP,
               '_PROD_':prodd,
               '_SCRIPTNAME_':basename(script),
               '_STARTFR_':startframe,
               '_ENDFR_':endframe,
               }
        for (key,val) in repl.items(): buf=buf.replace(key,val)
        currfr=currfr+dfr #advance deltaframe
        #capping
        if currfr<lastframe:
            startframe='%05d'%(currfr+1)
            endframe='%05d'%(currfr+dfr)
            nextscript=joink(scriptd,startframe+'_'+endframe+'.sub')
            buf+='qsub '+nextscript+'\n'
        buf+='exit 0\n'
        open(script,'w').write(buf)

    firstscript=joink(scriptd,'00001_'+'%05d'%(dfr)+'.sub')
    print 'type: qsub '+firstscript

    
if __name__=='__main__':
    ih=inpHand('Usage: qsubpmemd2.py',
               ' -a _RA_prodd absolute path to production directory',
               ' -b _R_TYP protein ID',
               ' -c _R_jobname jobname',
               ' -d _R_maxtime maxtime in nanoseconds',
               ' -e _R_ncores number of cores',
               )
    ih.parse(locals(),sys.argv)

    setUp(prodd,TYP,jobname,maxtime,ncores)
