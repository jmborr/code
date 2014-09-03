#!/bin/bash

#outputs have to be marked with the $header and $arfavor parameters
pdbheader=$1
arfavor=$2   #ratio of strengths of long range to short range hydrogen bond
casversion=$3
outd0=$4 #/gpfs1/scratch/jose/tuneHBstrength/out/$pdbheader
outd=$outd0/$arfavor
pcd=/gpfs1/active/jose/code/projects/tuneHBstrength
indd=$pcd/indd
currd=`pwd` #current directory
avHBn=$5    #filter trajectory by average number of hydrogen bonds

#input from output directory
cd $outd0
/bin/cp rmsinp seq.dat chain.dat $currd

#input parameter files and executables
cd $pcd
/bin/cp indd/$arfavor/in.dd $currd
/bin/cp $casversion spicker_gene $currd
/bin/cp $NCD/python/tertiary/ca2at.py $currd #amber minimize
/bin/tar jxf $CD/f77/tasser/tasser.1.0/common.tbz2 -C $currd

#input more executables
/bin/cp $CD/cpp/pdb_utils/Rg.x $currd
/bin/cp $CD/f77/tertiary/filterByHBnumber.x $currd

#go to working directory
cd $currd
echo "TASSER RUN"
if ./$casversion; then
  echo "FILTERING TRAJECTORIES BY Rg"
  ls -1 rep*.tra|while read f;do  #filter trajectories by radius of gyration
    echo "working on trajectory $f"
    ./Rg.x -c $f > /dev/null; #filter by Rg
    ./filterByHBnumber.x -train $f.filtered -traout $f.hb -avHBn $avHBn #filter by avergae HB numb
    /bin/mv -f $f.hb $f #overwrite trajectory file
  done
  ls -1 rep*.tra|wc -l > tra.in && ls -1 rep*.tra >> tra.in
  /bin/tar jcf out.tasser.avHBn.tbz2 rep*.tra tra.in
  /bin/mv out.tasser.avHBn.tbz2 $outd
  if ./spicker_gene; then
    /bin/tar jcf out.spicker.avHBn.tbz2 closc*.pdb combo*.pdb rst.dat
    /bin/mv out.spicker.avHBn.tbz2 $outd
    echo "AMBER MINIMIZE SELECTED CLOSC STRUCTURES"
    nmax=`ls -1 closc*.pdb|wc -l`
    m=6
    for n in 1 2 3 4 5 6 7 8 9 10;do
      echo "minimizing closc of cluster $n"
      ./ca2at.py -a closc${n}.pdb -b ${pdbheader}.${n}.pdb #dukka+pulchra+amber minimize
    done
    echo "moving output"
    /bin/tar jcf out.structures.avHBn.tbz2 ${pdbheader}.*.pdb
    /bin/mv out.structures.avHBn.tbz2 $outd #move generated structures
  else
    exit 1
  fi
else
  exit 1
fi

exit 0
