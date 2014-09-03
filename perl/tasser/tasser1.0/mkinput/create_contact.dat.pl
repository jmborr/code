#!/usr/bin/perl

use Getopt::Std;
use small_util qw(addAbsPath codedir db);

$currd=`pwd`; chomp($currd);
$outd=$currd;

&parse_input();
if(!defined($predictedf)){ $predictedf="$header.predictedrap3orienrev" ;}
if(! -e $outd){ system("mkdir -p $outd"); } #creates non-existent directory
if(!defined($outf)){ $outf="contact.dat"; }
$outf="$outd/$outf"; #add absolute path

open(RAW,"$predictedf");  open(NEAT,">$outf");
$l=<RAW>; print NEAT "$l"; #first line is different than rest
while($l=<RAW>){ print NEAT substr($l,0,25)."\n"; }
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
    print "Usage:create_contact.dat.pl \n";
    print " Required:\n";
    print " -i header : file-letter (xxxxx) pdb code\n";
    print " Optional:\n";
    print " -f predictedf : ${header}.predictedrap3orienrev type file\n";
    print " -o outd   : output directory (default=current directory)\n";
    print " -p outf   : output filename (default=contact.dat)\n";
    print "\n";
    exit(0);
}
