#!/bin/sh
if [ $# -ne 3 ]
then
	echo "need three parameters:seq_file(name,seq in compact format), align_dir, output file." 
	exit 1
fi
#assumption: alignment_file = align_dir + name
/gpfs1/active/jose/code/bin/sspro4/script/predict_seq_sa.pl /gpfs1/active/jose/code/bin/sspro4/server/predict_seq_sa.sh $1 $2 $3 
