#!/usr/bin/python
import re,os,sys,pdb

####################################
#                                  #
# NOTE: modify 'parms' AND 'opts'  #
# before running the script        #
#                                  #
####################################

"""
This script generates:
  - input files under each LAMBDA/curr directory
  - job submission scripts *.sub
  - Prod.sh script to send the previous job submission scripts to the queue

 Required directory hierarchy and files:
 Prod/
  |__sub/ (generate.py)
  |__files/ (ptraj.crd.in,ptraj.vel.in)
  |   |__templates/ (dh0.evb,groupfile,in,p.evb,r.evb,th0.evb)
  |__init/ (dh0.top,dhfr.pdb,dhfr.stripped.r.pdb,dhfr.stripped.r.top,in,p.top,README.txt,r.top,th0.top)
  |__L-0.00/
  |   |__curr/
  |   |__prev/
  |
  |__L-1.00/
      |__curr/
      |__prev/
"""

parms={'_TOTAL_NFRAMES_':2000, #total number of frames
       '_TOTAL_TIME_':200,     #total time, in picoseconds
       'TIME_PER_JOB':10,      #time for each job, in picoseconds
       '_DT_':0.001,           #time in between steps
       '_TEMP_':300,}

opts={'_DL_':0.05,    #increase in LAMBDA parameter
      '_PROD0_':'/home/jmborr/code/python/conductance_package/Prod',
      '_WALLTIME_':'01:00:00',
      '_NCPU_':'504', #total number of cores (a multiple of four)
      '_CPUPS_':'24', #cores per simulation
      '_TYP_':'dhfr'}

time_per_frame=float(parms['_TOTAL_TIME_'])/parms['_TOTAL_NFRAMES_']
frame_per_psec=parms['_TOTAL_NFRAMES_']/parms['_TOTAL_TIME_']
njobs=parms['_TOTAL_TIME_']/parms['TIME_PER_JOB']
nframes_per_job=parms['_TOTAL_NFRAMES_']/njobs
parms['_NSTLIM_']=parms['TIME_PER_JOB']/parms['_DT_'] #number steps per job
parms['_NTWX_']=int(time_per_frame/parms['_DT_']) #number steps per frame

#Series of questions to make sure all files are in place
ow=sys.stdout.write
ow('did you modify the initial parameters in this script?(y/n): ')
buf=sys.stdin.readline(); buf=buf[0].lower()
assert buf=='y', 'Sorry, do that first!'
ow('did you put the topology files under init/ ?(y/n): ')
buf=sys.stdin.readline(); buf=buf[0].lower()
assert buf=='y', 'Sorry, do that first!'
ow('did you modify initial parameters in script ener_redistr/doit.sh?(y/n): ')
buf=sys.stdin.readline(); buf=buf[0].lower()
assert buf=='y', 'Sorry, do that first!'
mssg="""Sorry, do that first!
How-to update 'mass' and 'sys':
1 Create a PDB (read rst with ptraj and save as pdb)
2 Create mass file with 'pdb_to_mass PDB-file mass'
3 In PDB file, set bfactor column to 100 for cold atoms
4 Run 'python gensys.py PDB-file' under ener_redistr/
"""
ow('did you updated files mass and sys under ener_redistr/ (y/n): ')
buf=sys.stdin.readline(); buf=buf[0].lower()
assert buf=='y', mssg
ow( 'Did you supply ref and rst files within each LAMBDA directory?(y/n): ' )
buf=sys.stdin.readline(); buf=buf[0].lower()
assert buf=='y', 'Sorry, do that first!'

#instantiate 'in' file
buf=open('../files/templates/in').read()
for (key,val) in parms.items(): buf=buf.replace(key,str(val))
open('../init/in','w').write(buf)

#helping globals
links={'../../init/r.top': 'r.top',
       '../../init/p.top': 'p.top',
       '../../init/in': 'in'}

#create Lambda directories, copy files, create soft links
L=0
while L<=1.01:
  print 'Preparing directory L-%4.2f'%L

  os.chdir( '%s/L-%4.2f/'%(opts['_PROD0_'],L) )

  #copy evb templates
  for s in ('r','p'):
    evbf='../files/templates/'+s+'.evb'
    buf=open(evbf).read().replace('_LAMBDA_','%4.2lf'%L)
    open('curr/%s.evb'%s,'w').write(buf)

  #copy "groupfile" files
  os.system('/bin/cp ../files/templates/groupfile curr/groupfile')

  #generate soft links
  for (target,link) in links.items():
    os.system('ln -s %s curr/%s'%(target,link))

  L+=opts['_DL_']
 
#template submission script
tpl0="""#!/bin/bash
#PBS -A BIP003
#PBS -N _TYP_
#PBS -j oe
#PBS -o _PROD0_/sub/${PBS_JOBID}.log
#PBS -l walltime=_WALLTIME_,size=_NCPU_

EFLOW=_PROD0_/ener_redistr/doit.sh
SANDER='aprun -n_CPUPS_ /ccs/home/jmborr/code/software/amber10/exe/sander'
#Energy(Kcal/mol) to remove from the hot-spot (24 atoms)
#deltaKs=(9.552 7.642 5.731 3.821 2.388 1.433 0.00)
deltaK=9.552

TYP=_TYP_
LAMBDAs=(0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00)
M=20 #number of LAMBDA values minus one

beg=_BEGIN_    #initial simulation frame
end=_END_      #final simulation frame
delta=_DELTA_  #number of frames in between two energy redistrib. events
let prev=beg
let next=beg+delta-1

while [ "$next" -le "$end" ];do
  range=`printf '%04d_%04d' $prev $next`
  
  for i in `seq 0 $M`;do #start jobs
    LAMBDA=${LAMBDAs[$i]}
    CURRD=_PROD0_/L-$LAMBDA/curr
    cd $CURRD
    /bin/cp rst rst.bak
    $EFLOW $CURRD/rst.bak $CURRD/rst $LAMBDA $deltaK
    $SANDER -ng 2 -groupfile groupfile &
    pids[$i]=$!
    echo "pids[$i]=${pids[$i]} TYP=$TYP $CURRD"
  done
  
  #check if simul finished. 
  for i in `seq 0 $M`;do #gather output
    wait ${pids[$i]} #wait until job finish
    LAMBDA=${LAMBDAs[$i]}
    CURRD=_PROD0_/L-$LAMBDA/curr
    PREVD=_PROD0_/L-$LAMBDA/prev
    cd $CURRD
    /bin/mv evbout ../prev/$range.evbout
    /bin/mv mdinfo ../prev/$range.mdinfo
    /bin/cp r.rst  ../prev/$range.r.rst
    /bin/mv r.out  ../prev/$range.r.out
    /bin/mv p.out  ../prev/$range.p.out
    /bin/mv r.crd  ../prev/$range.r.crd
    /bin/mv r.vel  ../prev/$range.r.vel
    /bin/mv r.rst rst
    cd $PREVD
    ln -s $range.r.rst $range.p.rst
    ln -s $range.r.rst $range.rst
    ln -s $range.r.crd $range.p.crd
    ln -s $range.r.crd $range.crd
    ln -s $range.r.vel $range.p.vel
    ln -s $range.r.vel $range.vel
    ln -s $range.r.out $range.out
  done

  let prev=next+1
  let next=next+delta
done

exit 0
"""

os.chdir( opts['_PROD0_']+'/sub' )
jn=os.path.join
subd=opts['_PROD0_']+'/sub'
iframe=1
eframe=nframes_per_job
range='%04d_%04d'%(iframe,eframe)
firstScript  = '%s.sub'%range
currscript = firstScript
scriptL=''

while eframe<=parms['_TOTAL_NFRAMES_']:
    #instantiate submission script from template
    tpl=tpl0
    for key,val in opts.items(): tpl=tpl.replace(key,str(val))
    for key,val in {'_BEGIN_':iframe,
                    '_END_':eframe,
                    '_DELTA_':frame_per_psec}.items():
      tpl=tpl.replace(key,str(val))
    #update range of simulation
    iframe=eframe+1
    eframe+=nframes_per_job
    range='%04d_%04d'%(iframe,eframe)
    #write script and update
    tpl=tpl.replace('_CURR_SCRIPT_',subd+'/'+currscript)
    nextscript='%s.sub'%range
    open(currscript,'w').write(tpl)
    scriptL+=' '+currscript
    currscript=nextscript

#template to create Bash script to submit the PBS jobs                                                                  
template="""#!/bin/bash
for script in _SCRIPT_LIST_;do
  if [ $script == _FIRST_SCRIPT_ ];then
    PBS_JOBID=$(qsub $script)
    echo \"qsub $script\"
  else
    PBS_JOBID=$(qsub -W depend=afterok:${PBS_JOBID} $script)
    echo \"qsub -W depend=afterok:${PBS_JOBID} $run.pbs\"
  fi
  sleep 3s
done
"""
tpl=template
tpl=tpl.replace('_SCRIPT_LIST_',scriptL)
tpl=tpl.replace('_FIRST_SCRIPT_',firstScript)
open('prod.sh','w').write(tpl)
os.system('chmod u+x prod.sh')
print '..done! Run ./prod.sh'
