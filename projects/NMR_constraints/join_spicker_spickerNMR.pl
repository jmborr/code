#!/usr/bin/perl

$fileSpicker1="/gpfs1/active/jose/projects/nonhomol_pdbJul05/out/summary.all";
$fileSpicker2="/gpfs1/active/jose/projects/nonhomol_pdbJul05/out_atlanta_cluster/summary.all";

open(IN,"$fileSpicker1"); @list =<IN>; close(IN);
open(IN,"$fileSpicker2"); @list2=<IN>; close(IN);

shift @list; shift @list2; #remove first line, which is a comment

for $item (@list){
    chomp($item);
    ($header,$rest)=split(' ',$item);
    for $item2 (@list2){
	($header2,$rest2)=split(' ',$item2);
	if($header eq $header2){
	    print "$item $item2";
	    last;
	}
    }
}
