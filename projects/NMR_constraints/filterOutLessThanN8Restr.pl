#!/usr/bin/perl

$rootd="/gpfs1/active/jose/projects/NMR_constraints/mkinput_NMR/run/input_output"; 
$listf="/gpfs1/active/jose/projects/NMR_constraints/tasser1.0_NMR_II/out/pdb1484_nonzeronoe.list"; #list of headers
open(IN,$listf); @list=<IN>; close(IN);


for $header (@list){
    chomp($header);
    $subd=substr($header,1,1);
    $file="$rootd/$subd/$header/noe.dat";
    open(IN,$file); chomp($nnoe=<IN>); close(IN);
    $file="$rootd/$subd/$header/CA";
    open(IN,$file); @list2=<IN>; close(IN);
    $nres=$#list2;
    #print "header=$header nres=$nres nnoe=$nnoe\n";
    if($nres<=$nnoe*8){ print "$header $nres $nnoe\n"; }
}
exit(0);
