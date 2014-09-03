#!/bin/bash

pdbheader=$1
cas=$2    #basename of tasser executable
outd=$3   #output directory
switch=$4 #simulation protocol code

SC=/gpfs1/scratch/jose
CD=/gpfs1/active/jose/code
pcd=$CD/projects/structural_neighbors

currd=`pwd` #current directory

#discern input to use
SCP=$SC/structural_neighbors/out/$pdbheader
rmsinp=$SCP/rmsinp ; seqdat=$SCP/seq.dat ; chaindat=$SCP/chain.dat
if [ $switch = megamix ];then
  SCP=$SC/structural_neighbors/out/$pdbheader/megamix
  seqdat=$SCP/seq.dat
  chaindat=$SCP/chain.dat 
fi

#discern which in.dd file to use
if [ $cas = casB_gene4 ];then   #select in.dd appropriate to TASSER flavor
  indd=$pdc/indd/2.0/in.dd
elif [ $cas = casB_gene5 ];then
  indd=$pcd/indd/2.0/eh6b/in.dd
  if [ ${outd:(-9):9} = /1.0/eh6b ];then
    indd=$pcd/indd/1.0/eh6b/in.dd
  fi
fi

#cp input and executables
if ! /bin/cp $rmsinp $seqdat $chaindat $indd $pcd/$cas $pcd/spicker_gene $CD/python/tertiary/ca2at.py $CD/cpp/pdb_utils/Rg.x $currd;then
  echo "submitGenericTasser.sh could not copy input or executables"
  exit 1
fi
/bin/tar jxf $CD/f77/tasser/tasser.1.0/common.tbz2 -C $currd

#go to working directory
cd $currd
echo "TASSER RUN"
if ./$cas; then
  /bin/tar jcf out.tasser.unfiltered.tbz2 rep*.tra
  /bin/mv out.tasser.unfiltered.tbz2 $outd
  echo "FILTERING TRAJECTORIES BY Rg"
  ls -1 rep*.tra|while read f;do  #filter trajectories by radius of gyration
    echo "working on trajectory $f"
    ./Rg.x -c $f > /dev/null; 
  done 
  ls -1 rep*.tra.filtered|wc -l > tra.in && ls -1 rep*.tra.filtered >> tra.in
  /bin/tar jcf out.tasser.tbz2 rep*.tra.filtered tra.in
  /bin/mv out.tasser.tbz2 $outd
  if ./spicker_gene; then
    /bin/tar jcf out.spicker.tbz2 closc*.pdb combo*.pdb rst.dat
    /bin/mv out.spicker.tbz2 $outd
    echo "AMBER MINIMIZE SELECTED CLOSC STRUCTURES"
    nmax=`ls -1 closc*.pdb|wc -l`
    for n in 1 2 3 4 5 6 7 8 9 10;do
      echo "minimizing closc of cluster $n"
      ./ca2at.py -a closc${n}.pdb -b ${pdbheader}.${n}.pdb #dukka+pulchra+amber minimize
    done
    echo "moving output"
    /bin/tar jcf out.minimized.tbz2 ${pdbheader}.*.pdb
    /bin/mv out.minimized.tbz2 $outd #move generated structures
  else
    exit 1
  fi
else
  exit 1
fi

exit 0
