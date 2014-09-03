#!/usr/bin/perl

use small_util qw(addAbsPath $codedir);
use Getopt::Std;

&parse_input();

open(in,$list); @headers=<in>; close(in);
open(out,">$outf");
for $header (@headers){
    print $header;
    chomp($header);
    $branch=substr($header,1,1)."/".$header;
    $noef="$root/$branch/noe.dat";
    @count=`cat $noef`;
    if(@count>1){ print out "$header\n"; }
}
sub parse_input(){
    my $f=0;
    if(!getopts('t:r:o:h')){ &message(); }
    if($opt_t){ $list=$opt_t; $f++; }
    if($opt_r){ $root=addAbsPath($opt_r); $f++; }
    if($opt_o){ $outf=addAbsPath($opt_o); $f++; }
    if($f<3){ &message(); } 
    if($opt_h){ &message(); }
}

sub message(){
    system("clear");
    print "Usage: filterNoeFiles.pl -t -r -o\n";
    print " Mandatory:\n";
    print " -t list: list of headers (five letter codes)\n";
    print " -r root: root directory where the find the NOE files (will append x/xxxxx to each of them)\n";
    print " -o out: outputfile\n";
    print "Will scan each \"noe.dat\" file and check if it has one or more restraints. In this case, its header will be added to the outputfile\n";
    print "\n";
    exit(0);
}
