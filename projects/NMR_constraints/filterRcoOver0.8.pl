#!/usr/bin/perl

$listf="/gpfs1/active/jose/projects/NMR_constraints/tasser1.0_NMR_II/out/pdb_overN8.list";
$rootd="/gpfs1/active/jose/projects/NMR_constraints/tasser1.0_NMR_II/out";

open(IN,"$listf"); @list=<IN>; close(IN); 
shift @list; #remove firs line, which is a comment

for $line (@list){
    $header=substr($line,0,5); 
    $subd=substr($header,1,1);
    #print "$rootd/$subd/$header/summary.txt\n";
    #print "$header\n";
    $summf="$rootd/$subd/$header/summary.txt";
    if( -e $summf){
	open(IN,$summf);
	do{
	    $summline=<IN>; 	#print $summline;
	}until($summline=~/select structures with rco >\s+([\d\.]+)/);    
	$rco=$1;    
	if($rco>=0.8){
	    #print "$header $rco $rootd/$subd/$header/summary.txt\n";
	    print "$rco $line";
	}
	close(IN);
    }
}
 
