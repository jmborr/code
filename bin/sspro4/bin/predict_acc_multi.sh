#!/bin/sh
#predict sa for a single sequence from scratch.
if [ $# -ne 3 ]
then
	echo "need 3 parameters:seq_file(in fasta format), output_file, threshold(0:0%,1:5%,...,5:25%,...19:95%)." 
	exit 1
fi
#output: a file with predicted ss and sa.
/gpfs1/active/jose/code/bin/sspro4/script/predict_acc_multi.pl /gpfs1/active/jose/code/bin/sspro4/blast2.2.8/ /gpfs1/scratch/jose/db/sspro4data/big/big_98_X /gpfs1/scratch/jose/db/sspro4data/nr/nr /gpfs1/active/jose/code/bin/sspro4/server/predict_seq_sa.sh /gpfs1/active/jose/code/bin/sspro4/script/ $1 $2 $3
