#!/usr/bin/perl

$rootd=$ARGV[0]; 
$listf=$ARGV[1]; #list of headers
open(IN,$listf); @list=<IN>; close(IN);


for $header (@list){
    chomp($header);
    $subd=substr($header,1,1);
    $file="$rootd/$subd/$header/noe.dat";
    open(IN,$file); chomp($line=<IN>); close(IN);
    if($line>0){print "$header\n";}
}
exit(0);
