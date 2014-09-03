#!/bin/bash

CURRD=/home/jmborr/code/python/conductance_package/Prod/ener_redistr
INPUTFILE=$1   #input  rst file
OUTPUTFILE=$2  #output rst file
SUFFIX=$3      #label to uniquely identify the run
deltaK=$4      #remove this many Kcal/mol from the binding site (24atoms)
TOLD=300
NATOMS=26490

MASS=$CURRD/mass
SYS=$CURRD/sys
suffix=$SUFFIX                     #label to uniquely identify the run
inputmp=$CURRD/input.$suffix       #temporary input file
outtmp=$CURRD/output.temp.$suffix  #temporary output file

echo $NATOMS >  $inputmp
cat $SYS     >> $inputmp
cat $MASS    >> $inputmp

tail --lines=+3 $INPUTFILE | head -n $NATOMS >> $inputmp
$CURRD/energy.x $deltaK < $inputmp > $outtmp

head -n 2 $INPUTFILE >  $OUTPUTFILE
cat $outtmp          >> $OUTPUTFILE
tail -n 1 $INPUTFILE >> $OUTPUTFILE

/bin/rm $inputmp $outtmp #rm temporary files

exit 0
