#!/bin/bash

PROJD=/projects/research/BaSO4/andrew_simulation
CODED=/home/jbq/code/projects/BaSO4

cd $PROJD/
#for temp in 230 240 260 282 300;do
for temp in 230 240 260;do
  #python $CODED/BaSO4.py preparing --kwargs "job=translate to DCD,pdb=q300.pdb,lammpstrj=q${temp}_nopbc.lammpstrj,dcd=q${temp}_nopbc.dcd"
  #ln -s q${temp}_nopbc.dcd q${temp}.dcd
  #cp sassena_inc_template.xml sassena_inc_T${temp}.xml
  #perl -p -i -e "s/_DCDFILE_/q$temp.dcd/g" sassena_inc_T${temp}.xml
  #perl -p -i -e "s/_FQTFILE_/fqt_inc_T${temp}_water.h5/g" sassena_inc_T${temp}.xml
  #/bin/cp sassena_inc_template.pbs sassena_inc_T${temp}.pbs
  #perl -p -i -e "s/_XMLFILE_/sassena_inc_T${temp}.xml/g" sassena_inc_T${temp}.pbs
  #perl -p -i -e "s/_OUTFILE_/sassena_inc_T${temp}.out/g" sassena_inc_T${temp}.pbs
  #perl -p -i -e "s/_TEMP_/$temp/g" sassena_inc_T${temp}.pbs
  python ~/code/python/mysassena/version.py --filename fqt_inc_T${temp}_water.h5
  python ~/code/python/mysassena/orderByQmodulus.py --filename fqt_inc_T${temp}_water.h5 --outfile junk.h5
  /bin/mv junk.h5 fqt_inc_T${temp}_water.h5
done

popwindow.py "finished BaSO4 terminal_commands.sh"
