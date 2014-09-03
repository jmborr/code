#!/usr/bin/perl
$execut="/home/jmborr/Projects/NMR_constraints/bin/outputNOE.x";
$data="/home/jmborr/Projects/NMR_constraints/data/";
$list="/home/jmborr/Projects/NMR_constraints/pdb1484.list";
open(in,$list);
while($d=<in>){
    chomp($d);
    $header=substr($d,2,5); 
    print "header=$header\n";
    $d=$data.$d."/";
    $pdb=$d.$header.".pdb";
    %outf=(
	'HVSC_HVSC', $d.$header.".HVSC_HVSC",
	'HVBK_HVBK', $d.$header.".HVBK_HVBK",
	'F_HVSC_S_HVBK', $d.$header."F_HVSC_S_HVBK",
	   );
    foreach $type (keys %outf){
	#print "$execut $pdb $type $outf{$type}\n";
	print "$header\n";
	system("$execut $pdb $type $outf{$type}");
    }

}
close(in);
