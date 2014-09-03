#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,junkName

codedir=os.environ['CODED']
isthere=os.path.exists
joink=os.path.join
basename=os.path.basename
template="""#!/bin/bash
#PBS -A BIP003
#PBS -N _JOBNAME_
#PBS -j oe
#PBS -o _LOGDIR_/${PBS_JOBID}._JOBNAME_.log
#PBS -l walltime=_WALLTIME_,size=_SIZETOTAL_

SANDER='aprun -n_SIZE_ /ccs/home/jmborr/code/software/amber9/exe/pmemd'

TYPs=(_TYPs_)
PRODs=(_PRODs_)
M=_M_

echo \"current script _SCRIPTNAME_\"
#submit simulations
for i in `seq 0 $M`;do
  TYP=${TYPs[$i]}
  PROD=${PRODs[$i]}
  cd $PROD/curr
  echo \"$SANDER -O -i in -p top -c rst -ref ref -o out -x crd -v vel -r rst.new &\"
  $SANDER -O -i in -p top -c rst -ref ref -o out -x crd -v vel -r rst.new &
  pids[$i]=$!
  echo \"pids[$i]=${pids[$i]} TYP=$TYP $PROD\"
  echo
done

#check if jobs finished.
for i in `seq 0 $M`;do
  TYP=${TYPs[$i]}
  PROD=${PRODs[$i]}
  wait ${pids[$i]}
  echo \"arranging $PROD\"
  cd $PROD/curr
  /bin/mv out     ../prev/$TYP._ENDFR_.out
  /bin/mv crd     ../prev/$TYP._STARTFR___ENDFR_.crd
  /bin/mv vel     ../prev/$TYP._STARTFR___ENDFR_.vel
  /bin/mv rst.new ../prev/$TYP._ENDFR_.rst
  /bin/rm rst
  /bin/ln -s      ../prev/$TYP._ENDFR_.rst rst
  echo $TYP._ENDFR_.rst >> ../rst.list
done
qstat -f ${PBS_JOBID}  #output job info
"""
def getwalltime(sizetotal):
    """return walltime according to jaguar scheduling policies"""
    if sizetotal<=256:    walltime='1:00:00'
    elif sizetotal<=512:  walltime='2:30:00'
    elif sizetotal<=1024: walltime='4:00:00'
    elif sizetotal<=6000: walltime='12:00:00'
    else: walltime='24:00:00'
    return walltime

ih=inpHand('Usage: qsubM2pmemd.py',
           ' -a _RA_jobf job file containing jobs to be performed',
           ' -b _R_service type of service (setProd|setScripts)',
           ' -c _R_ncores number of CPU per job',
           ' -d __jobname jobname (def: None)',
           ' -e __maxtime maxtime in nanoseconds',
           ' -f __lastframe maximum number of frames (lastframe)',
           ' -g __dfr number of frames per simulation (deltaframe)',
           )
ih.parse(locals(),sys.argv)


logdir=os.path.dirname(jobf)

#read jobfile
props={'TYP':[], 'PROD':[]}
N=0
for line in open(jobf,'r').readlines():
    if line[0]=='#': continue
    for item in line.split():
        prop,val=item.split('=')
        props[prop].append(val)
    N+=1
sizetotal=int(ncores)*N

if service=='setProd':
    from amber.qsubpmemd2 import setUp
    if not maxtime: Bye('provide maxtime')
    for i in range(N):
        PROD=props['PROD'][i]
        TYP=props['TYP'][i]
        retv=setUp(PROD,TYP,jobname,maxtime,ncores)
elif service=='setScripts':
    req={jobname:'need jobname',
         lastframe: 'need lastframe',
         dfr: 'need delta-frame',
         }
    for var,msg in req.items():
        if not var: Bye(msg)
    lastframe=int(lastframe) ; dfr=int(dfr)
    Nruns=lastframe/dfr #number of scripts
    currfr=0
    for i in range(Nruns):
        startframe='%05d'%(currfr+1)
        endframe='%05d'%(currfr+dfr)
        script=joink(logdir,startframe+'_'+endframe+'.sub')
        buf=template
        repl={ '_JOBNAME_':jobname,
               '_LOGDIR_':logdir,
               '_WALLTIME_':getwalltime(sizetotal),
               '_SIZETOTAL_':str(sizetotal),
               '_SIZE_':ncores,
               '_TYPs_' :' '.join(props['TYP']),
               '_PRODs_':' '.join(props['PROD']),
               '_M_':str(N-1),
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
            nextscript=joink(logdir,startframe+'_'+endframe+'.sub')
            buf+='qsub '+nextscript+'\n'
        buf+='exit 0\n'
        open(script,'w').write(buf)
    firstscript='00000_'+'%05d'%(dfr)+'.sub'
    print 'type:\n cd '+logdir+'\n qsub '+firstscript
