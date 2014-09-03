#!/bin/bash

dirname=`dirname $1`
filename=`basename $1`

if [ ! -e summary.txt ]; then
    echo $dirname
    cd $dirname 
    tar jxf $filename
    bzip2 -f rep*.tra
    cd -
fi