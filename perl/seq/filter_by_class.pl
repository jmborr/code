#!/usr/bin/perl

@class=('alph','beta',' a/b',' a+b','mult','memb','smal','coil','lowr',
	'pept','desg');

foreach $c (@class){ open($c,">pdb1484.${c}.list"); }
open(IN,"/home/jmborr/Databases/scop/pdb1484.class.list");
@cl=<IN>; close(IN);
open(IN,"/home/jmborr/Projects/NMR_constraints/pdb1484.list");
while($d=<IN>){
    print $d;
    $header=substr($d,2,5); #print "header=$header\n";exit;
    foreach $l (@cl){
	if($l=~$header){
	    $c=substr($l,7,4);  #print "class=|$class|\n"; exit;
	    last;
	}
    }
    print $c "$d";
}
foreach $c (@class){ close($c); }
close(IN);
