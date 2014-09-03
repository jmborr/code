#!/usr/bin/perl
$execut="/home/jmborr/Projects/NMR_constraints/bin/outputNOE.x";
$data="/home/jmborr/Projects/NMR_constraints/data";
$list="/home/jmborr/Projects/NMR_constraints/pdb1484.list";
open(in,$list);
@xx=(8); #alternatively, (4,8)
while($d=<in>){
    chomp($d);
    $header=substr($d,2,5); #print "header=$header\n";
    $d="$data/$d";
    %inf=(
	  "$d/101m_.HVSC_HVSC" => "$d/noe.dat"
	  );
    foreach $type (keys %inf){ #for each of the three types of NOEs.
	#print "$execut $pdb $type $outf{$type}\n";
	&permute($type);
    }
    open(OUT,">$d/mcnoe.dat"); print OUT "0\n"; close(OUT);
    open(OUT,">$d/smnoe.dat"); print OUT "0\n"; close(OUT);
}
close(in);

sub permute{
    $f=shift;    open(IN,$f);
    chomp($N=<IN>); $I=10*$N; #print "$N\n";
    @l=<IN>;     close(IN);
    for($i=0;$i<$I;$i++){ 
	$i1=int(rand($N)); $i2=int(rand($N));
	$t=$l[$i1]; $l[$i1]=$l[$i2]; $l[$i2]=$t; #swap order of these two lines
    }
    foreach $x (@xx){
	$out=$inf{$f}; open(OUT,">$out");    print "$out\n";
	$n=1+int($N/$x); $m=@l; if($m<$n){ $n=$m; } print OUT "$n\n";
	for($i=0;$i<$n;$i++){ print OUT $l[$i]; } #print only first $n lines
	close(OUT);
    }
}
