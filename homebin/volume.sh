#!/bin/bash
wd=`pwd`
vd="$CODED/software/voss_volume_voxelator"
/bin/cp $1 $vd/xyzr/junk.pdb
cd $vd/xyzr
./pdb_to_xyzr junk.pdb > junk.xyzr
x=`../bin/Volume.exe -i junk.xyzr -p 1.5 -g 0.5 2>/dev/null`
echo $x|cut -d ' ' -f 3
/bin/rm junk.pdb junk.xyzr
cd $wd
