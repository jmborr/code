#!/usr/bin/perl

$root="/home/jmborr/Projects/NMR_constraints";
$data="$root/data";
$list="$root/pdb1484.list";
open(in,"$root/analysis/lengths.dat"); @lengths=<in>;close(in);
open(in,"$list");
while($d=<in>){
    print "$d";
    chomp($d);
    $subdir=substr($d,0,1);
    $header=substr($d,2,5);
    $d_i="$data/$subdir/$header";
    foreach $x (@lengths){ #find length of protein
	if($x=~/$header/){
	    ($rest,$l)=split ' ', $x;
	    chomp($l);
	    last;
	}
    }
    $n=1+int($l/8);
    chomp($m=`head -1 $d_i/noe.dat`);
    if($n!=$m){printf "ERROR: n=$n m=$m in noe,  $header\n"; exit; }
    chomp($m=`head -1 $d_i/mcnoe.dat`);
    if($m!=0){printf "ERROR: m=$m in mcnoe, $header\n"; exit; }
    chomp($m=`head -1 $d_i/smnoe.dat`);
    if($m!=0){printf "ERROR: m=$m in smnoe, $header\n"; exit; }    
}
close(in);
