PROJD=/projects/research/LiCl
wd=/projects/research/LiCl/watBox30/equilibrate
wd=/global/homes/j/jmborr/m1503/jmborr/projects/LiCl/watBox30/production
for temp in 10 190 200 210 220 230 240 250 260 270 280 290;do
  #mkdir "T$temp"
  #cp -r ../annealing/toppar "T$temp/"
  #/bin/cp equilibrate.conf "T$temp/"
  #cd "$wd/T$temp/"
  #perl -p -i -e "s/_TEMP_/$temp/g" equilibrate.conf
  #python LiCl.py runSystem --kwargs "job=extract extended info from PDB header,inpdb=$PROJD/watBox30/equilibrate/T$temp/init.pdb,inxsc=$PROJD/watBox30/annealing/annealing.xsc,outxsc=$PROJD/watBox30/equilibrate/T$temp/init.xsc"
  #/bin/cp equilibrate/T$temp/equilibrate.restart.coor production/T$temp/production00.restart.coor
  #/bin/cp equilibrate/T$temp/equilibrate.restart.vel production/T$temp/production00.restart.vel
  #/bin/cp equilibrate/T$temp/equilibrate.restart.xsc production/T$temp/production00.restart.xsc
  #/bin/cp -r database db.xml sassenaInc.xml ../T$temp/
  #cd $PROJD/watBox30/production/T$temp; sleep 1s; sassena --config=sassenaInc.xml &
  #python $CODED/LiCl.py analyze --kwargs "job=parse out file,inOutFile=T$temp/equilibrate.out,outDatFile=T$temp/equilibrate.dat"
  #for irun in 0 1 2 3 4 5 6 7 8 ;do
  #  echo "T$temp $irun"; python $CODED/LiCl.py analyze --kwargs "job=parse out file,inOutFile=T$temp/production0$irun.out,outDatFile=T$temp/production0$irun.dat"
  #done
  #cd $PROJD/watBox30/production/T$temp; head -1 production00.dat > production.dat; grep -h -v TS production0?.dat >> production.dat
  #cd $PROJD/watBox30/production/T$temp; ptraj11 toppar/crd.md18_vmd_autopsf.psf < $PROJD/ptraj/rmsToFirst.ptraj
done

$CODED/LiCl.py runSystem --kwargs "job=extract initial conformations from annealing,anndir=$PROJD/watBox30/annealing,eqdir=$PROJD/watBox30/equilibrate”

for charge in 32 34 36 38 40 42 44 46 48 50;do
  #mkdir Q$charge
  #mkdir -p Q$charge/toppar
  #cp toppar/crd.md18_vmd_autopsf.pdb Q$charge/toppar/
  #cp toppar/par_all27_prot_lipid.inp Q$charge/toppar/
  #python LiCl.py runSystem --kwargs "job=change Hydrogen charge,Hq=0.$charge,inpsf=toppar/crd.md18_vmd_autopsf.psf,outpsf=Q$charge/toppar/Q$charge.psf"
  #/bin/cp /projects/research/LiCl/watBox30/equilibrate/T10/init* Q$charge/T10/equilibrate
  #/bin/cp /projects/research/LiCl/watBox30/equilibrate/T290/init* Q$charge/T290/equilibrate
  #/bin/cp -r database db.xml sassenaInc.xml ../Q$charge/T10/production/
  #/bin/cp -r database db.xml sassenaInc.xml ../Q$charge/T290/production/
  #for temp in 10 290;do
  #   echo "Q$charge/T$temp"
  #  /bin/cp Q$charge/T$temp/equilibrate/equilibrate.restart.coor Q$charge/T$temp/production/production00.restart.coor
  #  /bin/cp Q$charge/T$temp/equilibrate/equilibrate.restart.vel Q$charge/T$temp/production/production00.restart.vel
  #  /bin/cp Q$charge/T$temp/equilibrate/equilibrate.restart.xsc Q$charge/T$temp/production/production00.restart.xsc
  #  python $CODED/LiCl.py analyze --kwargs "job=parse out file,inOutFile=Q$charge/T$temp/equilibrate/equilibrate.out,outDatFile=Q$charge/T$temp/equilibrate/equilibrate.dat"
  #  cd $PROJD/watBox30/Hq/Q$charge/T$temp/production/; ptraj11 toppar/Q$charge.psf < $PROJD/ptraj/rmsToFirst.ptraj
  #  for irun in 01 02 03 04 05 06 07 08 09;do echo "Q$charge/T$temp/$irun"; python $CODED/LiCl.py analyze --kwargs "job=parse out file,inOutFile=Q$charge/T$temp/production/production$irun.out,outDatFile=Q$charge/T$temp/production/production$irun.dat";done
  #  python $CODED/LiCl.py analyze --kwargs "job=find average,inFiles=Q$charge/T$temp/production/production??.dat,column=15"
  #done
done

python $CODED/LiCl.py analyze --kwargs "job=parse out file,inOutFile=production00.out,outDatFile=production00.dat"

CODED=~/code/projects/LiCl/

for charge in 32 34 36 38 40 42 44 46 48 50;do
  cd /projects/research/LiCl/watBox30/Hq.II/Q$charge/T290/production
   python $CODED/LiCl.py prepareSystem --kwargs "job=create normalized convolution,resolution=$PROJD/EugeneData/resolution.nxs,simulated=simulated.nxs,convolved=convolved.nxs"
done

PROJD=/projects/research/LiCl/
CODED=~/code/projects/LiCl/
for charge in 421 422 423 424 425 426 427 428 429 430 431 432 433 434 435 436 437 438 439;do
 echo "charge=$charge"
 cd $PROJD/watBox30/Hq/NPT/Q$charge/T290/production
 # python /home/jmborr/code/python/mysassena/version.py --filename fqt_inc.hd5
 # python /home/jmborr/code/python/mysassena/version.py --filename fqt_coh.hd5

 # python /projects/development/CAMM/github/code/simulation/src/scattering/sassenatasks.py genSQE "fqt_inc.hd5 fqt_coh.hd5" simulated.nxs --wsname simulated --rebinQ "0.2 0.2 1.0" --LoadSassena "TimeUnit:1.0" --SassenaFFT "FFTonlyRealPart:1,DetailedBalance:1,Temp:290" --NormaliseToUnity "RangeLower:-0.1,RangeUpper:0.1"

 # python /projects/development/CAMM/github/code/simulation/src/beamline/convolve.py convolution --resolution=/projects/development/CAMM/github/code/simulation/test/LiCl/expdata/resolution.nxs --simulated=simulated.nxs --expdata=/projects/development/CAMM/github/code/simulation/test/LiCl/expdata/LiCl_290K.nxs --convolved=convolved.nxs

 # ptraj11 toppar/Q$charge.psf < $PROJD/ptraj/rmsToPrevious_equilibration_run.ptraj
 # ptraj11 toppar/Q$charge.psf < $PROJD/ptraj/rmsToPrevious.ptraj
 # ptraj11 toppar/Q$charge.psf < $PROJD/ptraj/rmsToFirst.ptraj
 # grep -P "ENERGY:" equilibrate.out|cut -d ':' -f 2 >> equilibrate.dat
 # echo -n "0.$charge "
 # python $CODED/LiCl.py analyze --kwargs "job=instant volume,xstfile=equilibrate.xst,volfile=volume.dat"
 # tar zxf production.xsc.tar.gz production00.restart.xsc
 # tail -1 production00.restart.xsc|tr -s ' '|cut -d ' ' -f 2,6,10; /bin/rm production00.restart.xsc
 # tar zxf production.out.tar.gz production01.out
 # echo "#      TS           BOND          ANGLE          DIHED          IMPRP               ELECT            VDW       BOUNDARY           MISC        KINETIC               TOTAL           TEMP      POTENTIAL         TOTAL3        TEMPAVG            PRESSURE      GPRESSURE         VOLUME       PRESSAVG      GPRESSAVG" > production01.dat
 # grep -P "ENERGY:" production01.out|cut -d ':' -f 2 >> production01.dat
 # /bin/rm production01.out
done


PROJD=/projects/research/LiCl/
for charge in 32 34 36 38 40 42 44 46 48 50;do
 cd $PROJD/watBox30/Hq.II/Q$charge/T290/production
 python /projects/development/CAMM/github/code/simulation/src/beamline/convolve.py convolution --resolution=/projects/development/CAMM/github/code/simulation/test/LiCl/expdata/resolution.nxs --simulated=simulated.nxs --expdata=/projects/development/CAMM/github/code/simulation/test/LiCl/expdata/LiCl_290K.nxs --convolved=/tmp/Q$charge/convolved.nxs
done

