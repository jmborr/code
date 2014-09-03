#!/usr/bin/perl

$rootd=$ARGV[0];
open(IN,"$ARGV[1]"); 
@list=<IN>; 
close(IN);

print "#header first_rank first_rcombo best_rank best_rcombo\n";
for $item (@list){
    chomp($item);
    $header = substr($item,3,5);
    $file="$rootd/$item"; #print "file=$file\n";
    open(IN,"$file"); @list2=<IN>;  close(IN);
    $i=0;
    for $line (@list2){
	if($line=~/A---/){last;}
	$i++;
    }
    $i++;
    ($rank,$first_rcombo,$junk)=split(' ',$list2[$i]);
    $best_rank=$rank;
    $best_rcombo=$first_rcombo;
    $i++;
    while(1){
	($rank,$rcombo,$junk)=split(' ',$list2[$i]);
	if($rank eq ''){last;}
	if($rcombo<$best_rcombo){
	    $best_rcombo=$rcombo;
	    $best_rank=$rank;
	}
	$i++;
    }
    printf "$header %2d %5.2lf %2d %5.2lf\n",
    1,$first_rcombo,$best_rank,$best_rcombo ;
}


