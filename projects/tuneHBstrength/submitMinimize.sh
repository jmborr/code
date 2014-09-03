#!/bin/bash

#outputs have to be marked with the $header and $arfavor parameters
pdbheader=$1
arfavor=$2   #ratio of strengths of long range to short range hydrogen bond
outd0=/gpfs1/scratch/jose/tuneHBstrength/out/$pdbheader
outd=$outd0/$arfavor
pcd=/gpfs1/active/jose/code/projects/tuneHBstrength
currd=`pwd` #current directory

#input from output directory
cd $outd
/bin/cp out.spicker.tbz2 $currd

#input parameter files and executables
cd $pcd
/bin/cp $NCD/python/tertiary/ca2at.py $currd #amber minimize

#go to working directory
cd $currd
/bin/tar jxf out.spicker.tbz2
echo "AMBER MINIMIZE SELECTED CLOSC STRUCTURES"
nmax=`ls -1 closc*.pdb|wc -l`
m=6
for n in 1 2 3 4 5 6 7 8 9 10;do
echo "minimizing closc of cluster $n"
./ca2at.py -a closc${n}.pdb -b ${pdbheader}.${n}.pdb #dukka+pulchra+amber minimize
done
echo "moving output"
/bin/tar jcf out.structures.tbz2 ${pdbheader}.*.pdb
/bin/mv out.structures.tbz2 $outd #move generated structures

exit 0
