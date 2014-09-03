 #!/usr/bin/perl
#===========================================================================
#Look at files xxxxx.HN and xxxxx.MCMC and obtain in total N/8 restrains, 
#where N is lenght of the protein, as follows
#(1) Shuffle restrains of xxxxx.HN
#(2) Pick as much as N/8 restraints (if possible) and save them as mcnoe.dat
#(3) If there were N/8 restraints in xxxxx.HN, then create noe.dat with zero
#    entries. If not, then suffle xxxxx.MCMC and pick as many restrains so 
#    that total number of restrains if N/8. If the number of restraints in
#    xxxxx.HN and xxxxx.MCMC is less than N/8, then output the protein name
#(4) create smnoe.dat with zero entries
#===========================================================================

$root="/home/jmborr/Projects/NMR_constraints";

#read protein lengths
open(in,"$root/analysis/lengths.dat");#length
while($line=<in>){ $prot2clength{substr($line,0,5)}=substr($line,7,3); }
close(in);

$list="$root/pdb1484.list";
#$list="$root/toylist"; #print "list=$list\n"; 
open(in,"$list");
while($d=<in>){
    $t=0; #current number of stored restraints
    $y=0; #number of restraints in the current file
    chomp($d); #print "d=$d\n";exit(0);
    $header=substr($d,2,5); #$print "header=$header\n";exit(0);

    $x=1+int( $prot2clength{$header} / 8.0 ); #number of needed restraints
    #print "x=$x\n";
    $d="$root/data/$d"; #print "d=$d\n";exit(0);
    @inf=("$d/$header.HN","$d/$header.MCMC"); #input files
    @outf=("$d/mcnoe.dat","$d/noe.dat","$d/smnoe.dat"); #output files
    
    #save mcnoe.dat
    @restraints=&permute($inf[0]);  #permute also sets up $x
    $y=@restraints;
    $t+=$y;
    open(out,">$outf[0]");  #print "$outf[0]\n";
    print out "$y\n";
    print out @restraints;     
    close(out);

    #save noe.dat
    open(out,">$outf[1]"); #print $outf[1];
    if($t<$x){ #more restraints needed
	@restraints=&permute($inf[1]);
	$y=@restraints;
	$t+=$y; #print "y=$y\n";exit(0);
	print out "$y\n" ;
	print out @restraints;    	
    }
    else{
	print out "0\n";
    }
    close(out);

     #save smnoe.dat
#    open(out,">$outf[2]");
#    if($t<$x){ #more restraints needed
#	@restraints=&permute($inf[2]);
#	$y=@restraints;
#	$t+=$y;
#	print out "$y\n" ;
#	print out @restraints;    	
#    }
#    else{
#	print out "0\n";
#    }
#    close(out);
    open(out,">$outf[2]");
    print out "0\n";
    close(out);

    if($t<$x){ print "$header insufficient number of restraints\n"; }
}
close(in);
#============================================================================
sub permute{
    my ($f,$N,$I,$i,$i1,$i2,$sw);
    my @l;
    $f=shift;    #print "f=$f\n";
    open(IN,$f);
    chomp($N=<IN>); $I=10*$N; #print "I=$I\n"; exit; #permute $I times 
    @l=<IN>;     #print @l; exit;
    close(IN);
    for($i=0;$i<$I;$i++){ 
	$i1=int(rand($N)); $i2=int(rand($N));
	$sw=$l[$i1]; $l[$i1]=$l[$i2]; $l[$i2]=$sw;#swap order of two lines
    }
    $z=$x-$t; #print "z=$z\n"; exit; #$z is number of restraints still needed
    if(@l<=$z){ return @l; }
    else{ return @l[0..$z-1]; } #return at most $N/8 elements
}
#============================================================================
