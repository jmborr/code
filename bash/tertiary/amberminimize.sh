#!/bin/bash
#
#PBS -l mem=150mb
#PBS -l ddisk=100mb
#PBS -l walltime=2:00:00

modelf=$1
modeln=`basename $modelf`
outdir=$2

######## copy what is needed ####################################

annad=/gpfs1/u/anna/amber8
if ! cp $annad/exe/sander $annad/exe/tleap $annad/dat/leap/cmd/leaprc.ff03  /gpfs1/active/lila/optim/all/anal/ambpdb .;then
  echo "amberminimize.sh could not copy input and executables"
  exit 1
fi

######## execute pulchra here ####################################

#pulchra $camodel -b && /bin/mv $camodel.rebuilt $camodel.p2
#dukka.x $camodel.p2 $camodel
#pulchra best.pdb  && /bin/mv best.pdb.rebuilt $camodel.rebuilt

######## prepare inputs for sander ####################################

/bin/cat>nat03.5.cmd<<EOF
logfile nat.log
source ./leaprc.ff03
set default PBradii mbondi2
x = loadpdb $modelf
saveamberparm x ./p03.5.for ./m.for
quit
EOF

echo "./tleap -s -f nat03.5.cmd > nat03.5.out"
./tleap -s -f nat03.5.cmd > nat03.5.out 

######## prepare inputs min #####################################
/bin/cat>min4r.in<<EOF
minimization in vacuum, eps=4r, amb8
 &cntrl
  ibelly=0, imin=1, ntr=1, irest=0,
  restraint_wt=100,
  restraintmask='@CA',
  maxcyc=100, ncyc=100, ntmin=1, dx0=0.01, dxm=0.5, drms=0.0001,
  ntx=1,
   igb=0, dielc=2.0,
  ntc=1,
  ntb=0,
  ntf=1,
  cut=20.0, scee=1.2, nsnb=10, scnb=2.0, ipol=0,
  temp0=300, tempi=0.0, ig=71277, heat=0.0, ntt=1, tautp=0.5,
  ntpr=5, ntave=0, ntwr=500, ntwx=500, ntwe=500, ntwv=0,
 &end
 &ewald
 eedmeth=5
 &end
EOF

/bin/cat>min.in<<EOF
minimization in GB/SA, amb8
 &cntrl
  ibelly=0, imin=1, irest=0,
  maxcyc=1000, ncyc=1000, ntmin=1, dx0=0.01, dxm=0.5, drms=0.0001,
  ntpr=50, ntave=0, ntwr=500, ntwx=500, ntwv=0, ntwe=500,
  ntx=1,    ntxo=1,
  ntc=1,
  ntb=0,
  ntf=1,
  cut=25.0, scee=1.2, dielc=1.0, scnb=2.0, nsnb=10, ipol=0,
  ntpr=100, ntr=0,
  igb=5, saltcon=0.0, intdiel=1.0, extdiel=78.5, gbsa=1
  temp0=300, tempi=0.0, ig=71277, heat=0.0, ntt=1, tautp=0.5,
  ntp=0,
 &end
EOF

echo "sander pass 1"
./sander -O -i min4r.in \
       -o min4r.out   \
       -p p03.5.for  \
       -c m.for  \
       -r min4r.for  \
       -ref m.for  \
       -inf mdinfo  

echo "sander pass 2"
./sander -O -i min.in \
       -o min.out   \
       -p p03.5.for  \
       -c min4r.for  \
       -r min.for   \
       -inf mdinfo  

echo "ambdb"
./ambpdb -p p03.5.for < min.for | /bin/sed "s;HIE;HIS;" > $outdir/$modeln.min

exit 0