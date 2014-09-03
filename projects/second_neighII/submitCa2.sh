#!/bin/bash

#input command-line arguments
pdbheader=$1
flavor=$2
#directories
scd=$SC/second_neighII/out
outd=$scd/$pdbheader/$flavor
currd=`pwd` #current directory
#executables and inputs
caslocal=cas.abinitio.x
cas=$CD/f77/tasser/tasser.1.0/abinitio/$flavor/$caslocal
indd=$scd/indds/$flavor/in.dd
spicker=/gpfs1/active/jose/code/f77/spicker/spicker0.0/spicker
ca2at=$NCD/python/tertiary/ca2at.py
#transfer inputs
/bin/tar jxf $scd/$pdbheader/in.tasser.tbz2 -C $currd
/bin/tar jxf $CD/f77/tasser/tasser.1.0/common.tbz2 -C $currd
/bin/cp $cas $indd $spicker $ca2at $currd

echo "TASSER RUN"
if ./$caslocal; then
  echo `ls -1 rep*.tra|wc -l` > tra.in &&  ls -1 rep*.tra >> tra.in
  /bin/tar jcf out.tasser.tbz2 rep*.tra tra.in
  /bin/mv out.tasser.tbz2 $outd
  if ./spicker; then
    /bin/tar jcf out.spicker.tbz2 closc*.pdb combo*.pdb rst.dat
    /bin/mv out.spicker.tbz2 $outd
    echo "AMBER MINIMIZE SELECTED CLOSC STRUCTURES"
    nc=`ls -1 combo*.pdb|wc -l`
    for n in `seq $nc`;do
      echo "minimizing combo of cluster $n"
      m=`printf %02d $n`
      ./ca2at.py -a combo${n}.pdb -b model.${m}.pdb #dukka+pulchra+amber minimize
    done
    echo "moving output"
    /bin/tar jcf out.structures.tbz2 model.*.pdb
    /bin/mv out.structures.tbz2 $outd #move generated structures
  else
    exit 1
  fi
else
  exit 1
fi

exit 0
