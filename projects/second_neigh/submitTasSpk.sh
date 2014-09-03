#!/bin/bash

#outputs have to be marked with the $header and $arfavor parameters
pdbheader=$1
inptas=$2
outd=$3

pcd=/gpfs1/active/jose/code/projects/second_neigh
indd=$pcd/indd
currd=`pwd` #current directory

#input parameter files and executables
/bin/tar jxf $inptas -C .
/bin/tar jxf $CD/f77/tasser/tasser.1.0/common.tbz2 -C $currd
/bin/cp $CD/f77/tasser/tasser.1.0/second_neigh/cas $CD/f77/tasser/tasser.1.0/second_neigh/in.dd /gpfs1/active/jose/code/f77/spicker/spicker0.0/spicker $NCD/python/tertiary/ca2at.py .

echo "TASSER RUN"
if ./cas; then
  
  echo `ls -1 rep*.tra|wc -l` > tra.in &&  ls -1 rep*.tra >> tra.in
  /bin/tar jcf out.tasser.tbz2 rep*.tra tra.in
  /bin/mv out.tasser.tbz2 $outd
  if ./spicker; then
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
  else
    exit 1
  fi
else
  exit 1
fi

exit 0
