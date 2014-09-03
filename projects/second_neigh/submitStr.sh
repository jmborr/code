#!/bin/bash

#outputs have to be marked with the $header and $arfavor parameters
pdbheader=$1
outspk=$2
outd=$3


pcd=/gpfs1/active/jose/code/projects/second_neigh
currd=`pwd` #current directory

#input parameter files and executables
/bin/tar jxf $outspk -C .
/bin/cp  $NCD/python/tertiary/ca2at.py .

echo "AMBER MINIMIZE SELECTED CLOSC STRUCTURES"
nc=`ls -1 closc*.pdb|wc -l`
for n in `seq $nc`;do
  echo "minimizing combo of cluster $n"
  ./ca2at.py -a combo${n}.pdb -b ${pdbheader}.${n}.pdb #dukka+pulchra+amber minimize
  done
  echo "moving output"
/bin/tar jcf out.structures.tbz2 ${pdbheader}.*.pdb
/bin/mv out.structures.tbz2 $outd #move generated structures

exit 0
