#!/bin/sh
if [ $# -ne 2 ]
then
	echo "need two parameters:seq_file, output_file." 
	exit 1
fi
/gpfs1/active/jose/code/bin/sspro4/script/generate_flatblast.pl /gpfs1/active/jose/code/bin/sspro4/blast2.2.8/ /gpfs1/active/jose/code/bin/sspro4/script/ /gpfs1/scratch/jose/db/sspro4data/big/big_98_X /gpfs1/scratch/jose/db/sspro4data/nr/nr $1 $2 
