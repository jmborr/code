#!/bin/bash

dirname=`dirname $1`
filename=`basename $1`

echo "cd $dirname"
cd $dirname
echo "tar jxf $filename"
tar jxf $filename
echo "bzip2 rep*.tra"
bzip2 rep*.tra
echo "cd -"
cd -
#echo "/bin/rm -f $1"
#/bin/rm -f $1