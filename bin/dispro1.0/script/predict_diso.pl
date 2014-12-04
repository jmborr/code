#! /usr/bin/perl -w

##################################################################
# Predict disorder region for one singel sequence in fasta format
# Input: predict_ssa.sh, diso predictor, model def file,  fasta file, output file
# Output: two output files: output_file (raw output from disopro), output_file.diso (disorder prediction in casp format)

# Author: Jianlin Cheng, July 4, 2004
#Modified from dopro/script/predict_dom.pl
##################################################################

if (@ARGV != 5)
{
	die "need 5 parameters: ssa predictor, diso predictor, model def, fasta input file, output file.\n"; 
}
$ssa_predictor = shift @ARGV;
$diso_predictor = shift @ARGV; 
$model_def = shift @ARGV;
$fasta_file = shift @ARGV;
$output_file = shift @ARGV; 
if (! -f $ssa_predictor)
{
	die "can't find the ss, sa predictor.\n"; 
}
if (! -f $diso_predictor)
{
	die "can't find domain predictor.\n"; 
}
if (! -f $model_def)
{
	die "can't find the model definition file.\n"; 
}
if (! -f $fasta_file)
{
	die "can't find the fasta file.\n"; 
}
open(FASTA, "$fasta_file") || die "can't open fasta file.\n"; 
$target_name = <FASTA>; 
close FASTA;
$target_name = substr($target_name,1); 

#generate alignment, predict ss and sa
$pos = rindex($fasta_file, "/"); 
if ($pos < 0)
{
	$ssa_file = $fasta_file; 
}
else
{
	$ssa_file = substr($fasta_file, $pos+1); 
}
$ssa_file =~ s/\./_/g; 
$ssa_file .= "ssa"; 
$align_file = "${ssa_file}align"; 
`$ssa_predictor $fasta_file $ssa_file`; 
#notice: two files are generated from ssa predictor: one is ssa output, one is alignment file. 

#create tmp file for disorder prediction
open(TMP, ">$output_file.tmp") || die "can't create temporary file.\n"; 
open(INPUT, "$ssa_file") || die "can't open the ssa file.\n"; 
<INPUT>;
$seq = <INPUT>; 
chomp($seq); 
$ss = <INPUT>; 
chomp($ss); 
$sa = <INPUT>; 
chomp($sa); 
close INPUT; 
$length = length($seq); 
if ($length != length($ss) || $length != length($sa))
{
	die "sequence length doesn't match.\n"; 
}

#make a tmp file for diso prediction
print TMP "1 25 2\n"; 
print TMP "$align_file\n$length\n";
print TMP "$seq\n$ss\n$sa\n";
$target = ""; 
for ($i = 0; $i < $length; $i++)
{
	$target .= "N"; 
}
print TMP "$target\n"; 
close TMP; 

#predict disorder  
system("$diso_predictor $model_def $output_file.tmp ./ > $output_file.res"); 
open(RES, "$output_file.res") || die "can't open result file.\n"; 
open(OUT, ">$output_file") || die "can't create output file.\n"; 
<RES>; 
$pre = <RES>; 
$pre =~ s/T/D/g;
$pre =~ s/N/O/g;
$probs = <RES>;
print OUT "$seq\n$pre$probs"; 
close RES;
close OUT; 

#remove temporary files
`rm $output_file.tmp $output_file.res`; 
`rm $ssa_file`; 
`rm $align_file`; 






