#!/bin/sh
#predict sa for a single sequence from scratch.
if [ $# -ne 2 ]
then
	echo "need 2 parameters:seq_file(in fasta format), output_file" 
	exit 1
fi
#output: a file with predicted ss and sa.
/gpfs1/active/jose/code/bin/sspro4/script/predict_acc_ab.pl /gpfs1/active/jose/code/bin/sspro4/blast2.2.8/ /gpfs1/scratch/jose/db/sspro4data/big/big_98_X /gpfs1/scratch/jose/db/sspro4data/nr/nr /gpfs1/active/jose/code/bin/sspro4/server/predict_seq_sa.sh /gpfs1/active/jose/code/bin/sspro4/script/ $1 $2 
