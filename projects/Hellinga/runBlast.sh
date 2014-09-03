#!/bin/bash

#variables
strIDs="a0512 a0653"
groupJobs="-p 20" #number of blast jobs to group in same qsub submission
extraFlags="" #"-e -Wx=qos:blast" #personal has maximum 40 running jobs
numDays=0.249
memReserv=900
fake="-o yes" #put this as one of the arguments of generic_job_list.py for fake submission
rootd=/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/Cys_rich_turn/omega_loop/out2

#This block for debugging purposes only
#g=a0212
#outdir=$rootd/$g/xxxxx
#list=$rootd/$g/toy.list
#exeline="/library/jose/code/python/combo_jobs/predictor/blasts.py -a $outdir/xxxxx -b $outdir"
#prefix="-n ${g:0:4}"
#generic_job_list.py  $fake $prefix $groupJobs $extraFlags -d $numDays -k $memReserv -c $outdir -a $list -g "$exeline"

#loop to send jobs
for g in $strIDs;do
  outdir=$rootd/$g/xxxxx
  list=$rootd/$g/list
  exeline="/library/jose/code/python/combo_jobs/predictor/blasts.py -a $outdir/xxxxx -b $outdir"
  prefix="-n i${g:2:1}" #neccessary to distinguish group jobs for different $strIDs
  generic_job_list.py  $prefix $groupJobs $extraFlags -d $numDays -k $memReserv -c $outdir -a $list -g "$exeline"
done

exit 0
