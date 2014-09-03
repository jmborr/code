#!/usr/bin/perl

$data="/home/jmborr/Projects/NMR_constraints/data";
$list="/home/jmborr/Projects/NMR_constraints/analysis/lengths.dat";

@noe=("mcnoe.dat","noe.dat");
@NOE=("HN","MCHN");
open(out,">relCO.dat");
print out "header length resid/restr   COall    CO   CO/COall\n";  
open(in,$list);
while($l=<in>){
    chomp($l);
    ($header,$L)=split(' ',$l); #print "l=$l header=$header L=$L\n";exit;
    $subdir=substr($l,1,1);

    $co_all=$n_all=0;
    foreach $suffix (@NOE){
	$f="$data/$subdir/$header/$header.$suffix"; #print "f=$f\n";
	open(in2,$f);
	$l2=<in2>;  ($i,$rest)=split(' ',$l2); #print "i=$i\n";
	if($i>0){
	    while($l2=<in2>){
		chomp($l2);
		($r1,$r2)=split(' ',$l2);
		$co_all+=abs($r1-$r2);
		$n_all++;
	    }
	}
	close(in2);
    }

    $co=$n=0;
    foreach $x (@noe){
	$f="$data/$subdir/$header/$x";
	open(in2,$f);
	$l2=<in2>; ($i,$rest)=split(' ',$l2);
	if($i>0){
	    while($l2=<in2>){
		chomp($l2);
		($r1,$r2)=split(' ',$l2);
		$co+=abs($r1-$r2);
		$n++;
	    }
	}
	close(in2);
    }
    
    print "$header $n_all\n";
    if($n_all>0 && $n>0){
#	printf out "$header %3d %3d %6.2lf %6.2lf\n",
#	$L, $n_all, $co_all, $co;
	printf out "$header    %3d     %3d       %6.2lf %6.2lf %6.2lf\n",
	$L,1+$L/$n,$co_all/($L*$n_all),$co/($L*$n),$co*$n_all/($co_all*$n); 
    }
    elsif($n_all>0 && $n==0){ #this is a bug !
	printf out "$header    %3d      INF      %6.2lf INF INF\n",
	$L,$co_all/($L*$n_all); 
    }
    else{
	printf out "$header    %3d      INF         INF    INF    INF\n",
	$L; 
    }
}
close(in);
close(out);
