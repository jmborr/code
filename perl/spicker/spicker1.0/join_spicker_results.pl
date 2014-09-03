#!/usr/bin/perl
use Getopt::Std;
use small_util qw(addAbsPath $codedir);

$nc=5;
$listf="none";
&parse_input();

#store list of headers in memory
if($listf eq "none"){
    @list=`find $root -name "summary.txt"`;
    for $ref (@list){ $ref=substr($ref,-20,7); } #reference to list elements
}
else{  
    open(in,$listf); @list=<in>; close(in); 
}
open(out,">$outf");
print out "#header toprank toprmsd firstrmsd \n";
foreach $l (@list){
    $subdir=substr($l,0,1);
    $header=substr($l,2,5); #print "$subdir $header\n"; exit;
    print "$header\n";
    #compose absolute path of "summary.txt" file
    $f="$root/$subdir/$header/summary.txt"; #print "f=$f\n";
    #output rmsd of first and top centroid. Output rand of top centroid
    open(in,"$f") or next;
    do($l=<in>)until($l=~'#cluster   rmsd');
    chomp($l=<in>);
    @x=split(' ',$l); $tr=$x[0]; $toprmsd=$firstrmsd=$x[1]; $n=1;
    if($firstrmsd=~'nan'){ print "ERROR: firstrmsd=nan for $header\n";exit;}
    chomp($l=<in>);
    while(!($l=~"rmsd of closests") && $n<$nc){#look first nc centroids
	@x=split(' ',$l); 
	$rank=$x[0]; $rmsd=$x[1];
	if($rmsd=~'nan'){ print "ERROR: rmsd=nan for $header\n";exit;}
	if($rmsd<$toprmsd){
	    $toprmsd=$rmsd;
	    $tr=$rank;
	}
	chomp($l=<in>);
	$n++;
    }
    printf out "$header     %2d     %5.2f     %5.2f\n",$tr,$toprmsd,$firstrmsd;
    close(in);
}
close(out);

sub parse_input(){
    my $f=0;
    if(!getopts('l:o:r:n:h')){ &message(); }
    if($opt_l){ $listf=$opt_l;                   }
    if($opt_n){ $nc=$opt_n;                      }
    if($opt_o){ $outf=$opt_o;              $f++; }
    if($opt_r){ $root=addAbsPath($opt_r);  $f++; }
    if($f<2){ &message(); }
    if($opt_h){ &message(); }
}

sub message(){
    system("clear");
    print "Usage: join_spicker_results.pl -l list -r root -o out\n";
    print " Mandatory:\n";
    print " -r root: root directory (root of outputs)\n";
    print " -o out: output file to store results\n";
    print " Optional:\n";
    print " -l list: list of type \"subdir/header\". If not specify, will use instead Unix \"find\" command to find all occurences of file \"summary.txt\"\n";
    print " -n nc : number of top clusters to consider(def=5)\n";
    print "\n";
    exit(0);
}

