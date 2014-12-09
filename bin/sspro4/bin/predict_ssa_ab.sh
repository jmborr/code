#!/bin/sh
#predict ss for a single sequence from scratch using ab-initio neural network.
if [ $# -ne 2 ]
then
	echo "need 2 parameters:seq_file(in fasta format), output_file." 
	exit 1
fi
#output: a file with predicted ss and sa.
/gpfs1/active/jose/code/bin/sspro4/script/predict_ssa_ab.pl /gpfs1/active/jose/code/bin/sspro4/blast2.2.8/ /gpfs1/scratch/jose/db/sspro4data/big/big_98_X /gpfs1/scratch/jose/db/sspro4data/nr/nr /gpfs1/active/jose/code/bin/sspro4/server/predict_seq_ss.sh /gpfs1/active/jose/code/bin/sspro4/script/ $1 $2 