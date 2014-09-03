#!/bin/bash

CODED=/home/jbq/code/projects/POSS/me8t8
PROJD=/projects/research/POSS/me8t8/charmm/nmolec_1

#for T in 10 50 75 100 125 150 175 200 225 250 300 350;do
 #cd $PROJD/T_$T
 #for file in me8t8.mdcrd  me8t8.pdb  me8t8.prmtop;do ln -s ../$file .;done
 #cp $PROJD/min0_template.conf $PROJD/T_$T/min0.conf
 #perl -p -i -e "s/_TEMPERATURE_/$T/g" $PROJD/T_$T/min0.conf
 #namd2-mc +idlepoll +p1 min0.conf> min0.log
 #cp $PROJD/run1_template.conf $PROJD/T_$T/run1.conf
 #perl -p -i -e "s/_TEMPERATURE_/$T/g" $PROJD/T_$T/run1.conf
 #namd2-mc +idlepoll +p1 run1.conf> run1.log &
 #sleep 2s
 #ptraj11 me8t8.prmtop < ../ptraj/rms2first_run1.in
 #python $CODED/me8t8.py single_molecule --kwargs="job=convert ptraj DCD to vmd DCD,prmtop=me8t8.prmtop,indcd=run1_rms2first.dcd,outdcd=junk.dcd"
 #/bin/mv  junk.dcd run1_rms2first.dcd
 #ln -s ../sassena_inc_run1_rms2first.xml .
 #ln -s ../sassena_coh_run1_rms2first.xml .
 #sassena --config=sassena_inc_run1_rms2first.xml
 #sassena --config=sassena_coh.xml
 #mv fqt_inc.hd5 fqt_inc.h5
 #python  ~/code/python/mysassena/version.py --filename=fqt_inc.h5
 #python  ~/code/python/mysassena/version.py --filename=fqt_coh.h5
 #python  ~/code/python/mysassena/orderByQmodulus.py --filename=fqt_inc_run1_rms2first.h5 --outfile=junk.h5
 #mv junk.h5 fqt_inc_run1_rms2first.h5
 #python  ~/code/python/mysassena/orderByQmodulus.py --filename=fqt_coh_run1_rms2first.h5 --outfile=junk.h5
 #mv junk.h5 fqt_coh_run1_rms2first.h5
 #mv fqt_inc.h5 fqt_inc.h5.BAK; mv fqt_Qordered_inc.h5 fqt_inc.h5
 #mv fqt_coh.h5 fqt_coh.h5.BAK; mv fqt_Qordered_coh.h5 fqt_coh.h5
 #mv fqt_inc.h5 fqt_inc_rms2first.h5
 #mv fqt_coh.h5 fqt_coh_rms2first.h5
 #ln -s ../sassena_inc_rms2first.xml .
 #ln -s ../sassena_coh_rms2first.xml .
 #mv fqt_coh_run1_rmsdfirst.h5 fqt_coh_run1_rms2first.h5
 #for dihedral in O11Si1C4H11 O12Si2C2H4 O5Si3C7H19 O4Si4C6H16 O9Si5C3H7 O7Si6C1H1 O12Si7C5H13 O8Si8C8H22;do
   #ptraj11 me8t8.prmtop < ../ptraj/dihedral_${dihedral}_run1_rms2first.in
 #  ls T_$T/dihedral_${dihedral}_run1_rms2first.in
 #done
 
#done

#cd $PROJD/
#for K in 0.02 0.04 0.06 0.07 0.08;do
#for K in 0.07;do
  #mkdir -p K$K/T_200
  #cp me8t8.prmtop K$K/
  #cd $PROJD/K$K/T_200/
  #ln -s ../me8t8.prmtop me8t8.prmtop
  #ln -s ../../me8t8.pdb 
  #ln -s ../../T_200/run0.restart.coor
  #ln -s ../../T_200/run1.conf
  #namd2-mc +idlepoll +p2 run1.conf > run1.log&
  #ptraj11 me8t8.prmtop < ../../ptraj/rms2first_run2.in
  #python $CODED/me8t8.py single_molecule --kwargs="job=convert ptraj DCD to vmd DCD,prmtop=me8t8.prmtop,indcd=run2_rms2first.dcd,outdcd=junk.dcd"
  #/bin/mv junk.dcd run2_rms2first.dcd
  #python ~/code/python/mysassena/version.py --filename fqt_inc_run1_rms2first.h5
  #python ~/code/python/mysassena/orderByQmodulus.py --filename fqt_inc_run1_rms2first.h5 --outfile junk.h5
  #mv junk.h5 fqt_inc_run1_rms2first.h5
#done

cd /projects/research/POSS/me8t8/amber/K0.07
sander -O -i min0.mdin -c me8t8.mdcrd -p me8t8.prmtop -o min0.mdout -x min0.mdcrd -r min0.restrt -e min0.mden -inf min0.mdinfo
nmode -O -i nmode0.nmdin -c min0.restrt -p me8t8.prmtop -o nmode0.nmdout -v nmode0.vecs
#python $CODED/me8t8.py mantidjobs --kwargs="job=create experimental Nexus files"

PROJD=/projects/research/POSS/me8t8
cd $PROJD/preparing/
python $CODED/me8t8.py preparing --kwargs "job=translate PDB from CSD to AMBER,inpdb=mercury_3x3x3.pdb,outpdb=me8t8_3x3x3.pdb"
python $CODED/me8t8.py preparing --kwargs="job=neutralize mol2 file,inmol2=me8t8.mol2,outmol2=me8t8_Qneutral.mol2"

