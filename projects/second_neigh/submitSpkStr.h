#!/bin/bash

#outputs have to be marked with the $header and $arfavor parameters
pdbheader=$1
outtas=$2
outd=$3
intas=$4

pcd=/gpfs1/active/jose/code/projects/second_neigh
currd=`pwd` #current directory

#input parameter files and executables
/bin/tar jxf $outtas -C .
/bin/cp  $NCD/python/tertiary/ca2at.py $CD/f77/spicker/spicker0.0/spicker .

if ! test -e combo1.pdb; then
  echo "run SPICKER"
  ls -1 rep*tra|wc -l > tra.in
  ls -1 rep*tra >> tra.in
  /bin/tar jxf $intas rmsinp seq.dat
  ./spicker
fi
echo "COLLECT SPICKER"
/bin/tar jcf out.spicker.tbz2 closc*.pdb combo*.pdb rst.dat
/bin/mv out.spicker.tbz2 $outd
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
