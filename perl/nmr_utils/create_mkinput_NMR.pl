#!/usr/bin/perl

$root="/home/jmborr/Projects/NMR_constraints";
$data="$root/data";
$image="$root/image";
`mkdir -p $image`;
$list="$root/pdb1484.list";
open(in,"$list");
while($d=<in>){
    print $d;
    chomp($d);
    $subdir=substr($d,0,1);
    `mkdir -p $image/$subdir`;
    $header=substr($d,2,5);
    $d_o="$image/$subdir/$header";
    `mkdir -p $d_o`;
    $d_i="$data/$subdir/$header";
    `cp $d_i/noe.dat $d_o/noe.dat`;
    `cp $d_i/mcnoe.dat $d_o/mcnoe.dat`;
    `cp $d_i/smnoe.dat $d_o/smnoe.dat`;
    `cp $d_i/CA $d_o/CA`;
}
close(in);
