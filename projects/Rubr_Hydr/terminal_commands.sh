#!/bin/bash

CODED=/home/jbq/code/projects/Rubr_Hydr
PROJD=/projects/research/Rubr_Hydr


for type in H-Rub_D2O D-Rub_H2O;do
  for hydr in h_20 h_37;do
    for T in T_25 T_50 T_75 T_100 T_125 T_150 T_175 T_200 T_225 T_250 T_275 T_300 T_325 T_350;do
      cd $PROJD/out/$type/$hydr/$T/solv
      #ln -s dcd solv.dcd
      #ln -s top solv.prmtop
      #ln -s ../../solv/pdb solv.pdb
      #python $CODED/Rubr_Hydr.py sassena --kwargs "job=convert ptraj DCD to vmd DCD,prmtop=solv.prmtop,indcd=solv.dcd,outdcd=solv_vmd.dcd"
      #ln -s ../../solv/hd.pdb solv_hd.pdb
      #ln -s ../../../../sassena/sassena_inc.xml .
      #ln -s ../../../../sassena/sassena_coh.xml .
      #python ~/code/python/mysassena/version.py --filename fqt_inc.h5
      #python ~/code/python/mysassena/version.py --filename fqt_coh.h5
    done
  done
done

#cd $PROJD/out
#for typefile in solv_hd.pdb solv_vmd.dcd sassena_coh.xml sassena_inc.xml;do
#for typefile in solv_hd.pdb;do
#  for file in `find . -name $typefile`;do
#    tar ufh junk.tar $file 
#  done
#done

