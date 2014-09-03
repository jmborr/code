#!/bin/bash

#variables
groupJobs="-p 50" #number of blast jobs to group in same qsub submission
extraFlags="-e -Wx=qos:personal" #personal has maximum 40 running jobs
numDays=0.249
memReserv=900
fake="-o yes" #put this as one of the arguments of generic_job_list.py for fake submission
rootd=/gpfs1/scratch/jose/missing_coordinates
inpdir=$rootd/input

#loop to send jobs
outdir=$rootd/out/xxxxx
list=$rootd/preparing/filterout500.list
exeline="/library/jose/code/python/combo_jobs/predictor/blasts.py -a $inpdir/xxxxx -b $outdir"
prefix="-n z" #neccessary to distinguish group jobs for different $strIDs
generic_job_list.py  $prefix $groupJobs $extraFlags -d $numDays -k $memReserv -c $outdir -a $list -g "$exeline"

exit 0
