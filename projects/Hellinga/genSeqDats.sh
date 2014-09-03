#!/bin/bash

#variables
strIDs="a0651  a1172  a1213  c0130  c1037  c1154"
groupJobs="-p 200" 
extraFlags="-e -Wx=qos:blast"
numDays=0.5
memReserv=2040
fake="-o yes"

#loop to send jobs
for g in $strIDs;do
  rootd=/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/hongpred
  outdir=$rootd/$g/xxxxx
  list=$rootd/$g/list
  exeline="/library/jose/code/python/secondary/mkseq.py -a $outdir/xxxxx -o $outdir"
  prefix="-n ${g:0:4}"
  generic_job_list.py  $prefix $groupJobs $extraFlags -d $numDays -k $memReserv -c $outdir -a $list -g "$exeline"
done

exit 0
