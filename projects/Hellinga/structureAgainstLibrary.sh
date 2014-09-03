#!/bin/bash

groupJobs="-p 50" #50 jobs per qsub
extraFlags="-e -Wx=qos:critical"
numDays=0.5
memReserv=200
lib=/tmp/jose/filteredHard:/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/filteredHard.tbz2
list2=a.list #compare only to helical structures
for g in a0112 a0216 a0864;do
  rootd=/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/sequence_checking/$g
  outdir=$rootd/xxxxx/abinitio
  list=$rootd/abinitio.list
  for x in `seq 10`;do
    xx=`printf "%02d" $x`
    prefix="x${g:4:1}$x"

    inputline="/bin/cp /gpfs1/active/jose/code/f77/TM-align/TMalign $rootd/xxxxx/abinitio/combo0$x.pdb.rebuilt ."
    exeline="/library/jose/code/python/TMalign/findBestPDBhitByTMalign.py -a combo$xx.pdb.rebuilt -e /tmp/jose/filteredHard/$list2 -d /tmp/jose/filteredHard/yyyyy -i all > all.dat"
    outputline="/bin/mv all.dat $rootd/xxxxx/abinitio/combo0$x.instructures.dat"
  
    generic_job_list.py $groupJobs $extraFlags -n $prefix -d $numDays -k $memReserv -m $lib -c $outdir -a $list -f "$inputline" -g  "$exeline" -i "$outputline"

  done
done

exit 0
