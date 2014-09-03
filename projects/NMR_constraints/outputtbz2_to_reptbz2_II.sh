#!/bin/bash

dirname=`dirname $1`
filename=`basename $1`

cd $dirname 

if [ -e rep1.tra.bz2 ]; then
    cd -
    exit
elif [ -e rep1.tra ]; then
    cd -
    exit
else
    echo  "cd $dirname"
    tar jxf $filename
    bzip2 rep*.tra
    echo "cd -"
    cd -
fi