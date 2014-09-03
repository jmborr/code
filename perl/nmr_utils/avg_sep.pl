#!/usr/bin/perl
$data="/home/jmborr/Projects/NMR_constraints/data/";
$sp=$n0;
while($d=<>){
    chomp($d);
    $header=substr($d,2,5); print "$header\n";
    $f=$data.$d."/".$header.".NOE"; #print "$f\n";
    open(IN,$f);  $gbg=<IN>; @l=<IN>;   close(IN);
    foreach $pair (@l){
	$i1=substr($pair,0,10);
	$i2=substr($pair,10,10); 
	$sp+=abs($i2-$i1); #print"sp=$sp\n";
	$n++;
    }
}
print $sp/$n;
