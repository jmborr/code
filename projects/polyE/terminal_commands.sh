#!/bin/bash

CODED=/home/jbq/code/projects/polyE
PROJD=/projects/research/polyE/nchain_1/

cd $PROJD/charmm
#ptraj11 E496.psf ptraj/run1_rms2first.dcd2crd_H.in
#python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/run1_rms2first_H.crd,Nhydr=2978,storagedir=$PROJD/isf_run1_rms2drist_H,isffile=$PROJD/isf_run1_rms2drist_H.dat"

#ptraj11 E496.psf < ptraj/run1_rms2first.dcd2crd_HA.in
#python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/run1_rms2first_HA.crd,dt=0.1,tl=10.0,Nhydr=496,storagedir=$PROJD/isf_run1_rms2drist_HA,isffile=$PROJD/isf_run1_rms2drist_HA.dat"

#qsub sassena_inc_run2_rms2first.pbs #produces fqt_inc_run2_rms2first.h5
#ptraj11 E496.psf ptraj/rms2first_run2.in #produces run2_rms2first.dat and run2_rms2first.dcd
#ptraj11 E496.psf < ptraj/run2_rms2first.dcd2crd_HA.in
#python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/run2_rms2first_HA.crd,qlist=0.5 0.6 0.7 0.8 0.9 1.0,dt=0.01,tl=10.0,Nhydr=496,storagedir=$PROJD/isf_run2_rms2first_HA,isffile=$PROJD/isf_run2_rms2first_HA.dat"

#check if RMS of all heavy atoms remove diffussivity of side chains
#ptraj11 E496.psf < ptraj/rms2first_heavyatom_run2.in
#ptraj11 E496.psf < ptraj/run2_rms2first_heavyatom.dcd2crd_HA.in
#python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/run2_rms2first_heavyatom_HA.crd,qlist=0.5 0.6 0.7 0.8 0.9 1.0,dt=0.01,tl=10.0,Nhydr=496,storagedir=$PROJD/isf_run2_rms2first_heavyatom_HA,isffile=$PROJD/isf_run2_rms2first_heavyatom_HA.dat"

# check if Boson peak enhanced in the NVE ensemble
#ptraj11 E496.psf < ptraj/rms2first_run3.in
#ptraj11 E496.psf < ptraj/run3_rms2first.dcd2crd_HA.in
#time python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/run3_rms2first_HA.crd,qlist=0.5 0.6 0.7 0.8 0.9 1.0,dt=0.01,tl=40.0,Nhydr=496,isffile=$PROJD/isf_run3_rms2first_HA.dat,sffile=$PROJD/sf_run3_rms2first_HA.dat"

#ptraj11 E496.psf < ptraj/rms2first_run4.in
#ptraj11 E496.psf < ptraj/run4_rms2first.dcd2crd_HA.in
#time python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/run4_rms2first_HA.crd,qlist=0.5 0.6 0.7 0.8 0.9 1.0,dt=0.01,tl=90.0,Nhydr=496,isffile=$PROJD/isf_run4_rms2first_HA.dat,sffile=$PROJD/sf_run4_rms2first_HA.dat"

#ptraj11 E496.psf < ptraj/rms2first_run5.in
#ptraj11 E496.psf < ptraj/run5_rms2first.dcd2crd_HA.in
#time python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/charmm/run5_rms2first_HA.crd,qlist=0.5 0.6 0.7 0.8 0.9 1.0,dt=1.0,tl=2000.0,Nhydr=496,isffile=$PROJD/charmm/isf_run5_rms2first_HA.dat,sffile=$PROJD/charmm/sf_run5_rms2first_HA.dat"


cd $PROJD/hydr40
#packmol < packmol_h0.40.inp > packmol_h0.40.out &
#ptraj11 E496.psf < ptraj/rms2first_run1.in
#ptraj11 E496.psf < ptraj/run1_rms2first.dcd2crd_HA.in
#time python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/hydr40/run1_rms2first_HA.crd,qlist=0.5 0.6 0.7 0.8 0.9 1.0,dt=0.01,tl=90.0,Nhydr=496,isffile=$PROJD/hydr40/isf_run1_rms2first_HA.dat,sffile=$PROJD/hydr40/sf_run1_rms2first_HA.dat"
#ptraj11 E496.psf < ptraj/rms2first_run2.in
#ptraj11 E496.psf < ptraj/run2_rms2first.dcd2crd_HA.in
time python $CODED/polyE.py scattering --kwargs="job=isf,crdfile=$PROJD/hydr40/run3_rms2first_HA.crd,qlist=0.5 1.0,dt=0.01,tl=1000.0,Nhydr=496,isffile=$PROJD/hydr40/isf_run3_rms2first_HA.dat,sffile=$PROJD/hydr40/sf_run3_rms2first_HA.dat"
