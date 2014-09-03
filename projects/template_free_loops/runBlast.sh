#!/bin/bash

#variables
groupJobs="-p 151" #number of blast jobs to group in same qsub submission
extraFlags="-e -Wx=qos:personal" #personal has maximum 40 running jobs
numDays=0.99
memReserv=900
fake="-o yes" #put this as one of the arguments of generic_job_list.py for fake submission
rootd=/gpfs1/scratch/jose/template_free_loops
fastadir=/gpfs1/archive/adrian/dat/pdb_mar07/input
genpy=/library/jose/code/python/combo_jobs/generic_job_list.py

#loop to send jobs
outdir=$rootd/out/xxxxx
list=$rootd/randomized.list
exeline="/library/jose/code/python/combo_jobs/predictor/blasts.py -a $fastadir/xxxxx -b $outdir"
prefix="-n x" 
$genpy  $prefix $groupJobs $extraFlags -d $numDays -k $memReserv -c $outdir -a $list -g "$exeline"

exit 0
