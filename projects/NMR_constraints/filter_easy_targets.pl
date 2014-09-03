#!/usr/bin/perl

$listf="/gpfs1/active/jose/projects/NMR_constraints/tasser1.0_NMR_II/out/pdb_overN8_rbelow0.8.list";
$rootd="/gpfs1/active/jose/projects/NMR_constraints/tasser1.0_II/out";

open(IN,"$listf"); @list=<IN>; close(IN); 
shift @list; #remove firs line, which is a comment

for $line (@list){
    $header=substr($line,0,5); 
    $subd=substr($header,1,1);
    open(IN,"$rootd/$subd/$header/chain.dat");
    $iseasy=<IN>; 
    #print "$header $iseasy $rootd/$subd/$header/chain.dat \n";
    if($iseasy=~'easy'){ print $line; }
    close(IN);
}
