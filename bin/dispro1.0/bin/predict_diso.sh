#!/bin/sh
#predict disorder region for one sequence.
if [ $# -ne 2 ]
then
	echo "need 2 parameters:fasta seq_file, output prefix."
	exit 1
fi
/gpfs1/active/jose/code/bin/dispro1.0/script/predict_diso.pl /gpfs1/active/jose/code/bin/sspro4/bin/predict_ss_sa.sh /gpfs1/active/jose/code/bin/dispro1.0/server/predict_seq /gpfs1/active/jose/code/bin/dispro1.0/model/model.def $1 $2 
