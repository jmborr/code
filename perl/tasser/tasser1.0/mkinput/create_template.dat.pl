#!/usr/bin/perl


use Getopt::Std;
use small_util qw(addAbsPath codedir db);

$currd=`pwd`; chomp($currd);
$outd=$currd;

&parse_input();
if(!defined($predictedf)){ $predictedf="${header}rap3orienrev1.pdb" ;}
if(! -e $outd){ system("mkdir -p $outd"); } #creates non-existent directory
if(!defined($outf)){ $outf="template.dat"; }
$outf="$outd/$outf"; #add absolute path

open(RAW,"$predictedf");  open(NEAT,">$outf");
$l=`grep "TER" $raw | wc`;
split(' ',$l);
print NEAT "$_[0]\n"; #number templates
# print "number templates=$_[0]\n";exit;
$c=1; #template number
while($l=<RAW>){
#print $l;
    if($l=~/TER/){ print NEAT "TER\n"; $c++; }
    elsif($l=~/ATOM/){ #print $l;
	@it=split(' ',$l); #items in $l
	printf NEAT "ATOM %6d  CA  $it[3]%6d    %8.3f%8.3f%8.3f%5d\n",
	$it[4],$it[4],$it[5],$it[6],$it[7],$it[8] ;
    }
    else{
	@it=split(' ',$l); #items in $l
	printf NEAT "%6d%12.5f%7d     $it[0]\n", $it[1],$it[2],$c,
    }
}
close(RAW); close(NEAT);

sub parse_input(){
    my $f=0;
    if(!getopts('i:f:o:p:h')){ &message(); }
    if($opt_i){ $header =$opt_i; $f++; }
    if($opt_f){ $predictedf =addAbsPath($opt_f); }
    if($opt_o){ $outd =addAbsPath($opt_o); }
    if($opt_p){ $outf =$opt_p; }
    if($opt_h){ &message(); }
    if($f<1){ $message(); }
}

sub message(){
    system("clear");
    print "Usage:create_template.dat.pl \n";
    print " Required:\n";
    print " -i header : file-letter (xxxxx) pdb code\n";
    print " Optional:\n";
    print " -f predictedf : ${header}rap3orienrev1.pdb type file\n";
    print " -o outd   : output directory (default=current directory)\n";
    print " -p outf   : output filename (default=template.dat)\n";
    print "\n";
    exit(0);
}
