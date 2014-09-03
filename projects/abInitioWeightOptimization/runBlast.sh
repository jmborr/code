#!/bin/bash

#variables
groupJobs="-p 25"
extraFlags="-e -Wx=qos:blast"
numDays=0.08 #slightly below 2 hours
memReserv=900
fake="-o yes"
rootd=/gpfs1/scratch/jose/spliceVariants/benchmark_PDB300_disordered_regions

#loop to send jobs
outdir=$rootd/blastjobs
list=$rootd/pdb300.list #toy.list
exeline="/library/jose/code/python/combo_jobs/predictor/blasts.py -a $rootd/input/xxxxx -b $outdir"
prefix="-n dis"
generic_job_list.py  $prefix $groupJobs $extraFlags -d $numDays -k $memReserv -c $outdir -a $list -g "$exeline"

exit 0
