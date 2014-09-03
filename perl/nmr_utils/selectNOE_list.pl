#!/usr/bin/perl

#apply selectNOE.pl to a list of headers

use small_util qw(addAbsPath $codedir);
use Getopt::Std;

$e=10;
$noefile='xxxxx.HVSC_HVSC';
$m=8;
&parse_input();

open(in,$list); @headers=<in>; close(in);

for $header (@headers){
    print $header;
    chomp($header);
    $branch=substr($header,1,1)."/".$header;
    $noefile="$header.HVSC_HVSC";
    $noef="$root/$branch/$noefile";
    `/bin/mkdir -p $root/$branch`;
    $outfile="$outd/$branch/noe.dat";
    chomp($l=`$codedir/python/seq/length.py -i $header`);
    ($junk,$l)=split(' ',$l);
    #print "$codedir/perl/nmr_utils/selectNOE.pl -i $noef -o $outfile -m $m -l $l -e $e -a $header\n";
    system("$codedir/perl/nmr_utils/selectNOE.pl -i $noef -o $outfile -m $m -l $l -e $e -a $header");
}


sub parse_input(){
    my $f=0;
    if(!getopts('t:i:r:o:m:e:h')){ &message(); }
    if($opt_t){ $list=$opt_t; $f++; }
    if($opt_i){ $inf=$opt_i; }
    if($opt_r){ $root=addAbsPath($opt_r); $f++; }
    if($opt_o){ $outd=addAbsPath($opt_o); $f++; }
    if($opt_m){ $m=$opt_m; }
    if($opt_e){ $e=$opt_e; }
    if($f<3){ &message(); } 
    if($opt_h){ &message(); }
}

sub message(){
    system("clear");
    print "Usage: selectNOE_list.pl -t -i -r -o  -m -l-e\n";
    print " Mandatory:\n";
    print " -t list: list of headers (five letter codes)\n";
    print " -r root: root directory where the find the NOE files (will append x/xxxxx to each of them)\n";
    print " -o outd: root output directory (will append x/xxxxx)\n";
    print " Optional:\n";
    print " -i noefile: list of NOE's, and first line is N, number of NOE's (def='101m_.HVSC_HVSC')\n";
    print " -m m: select 1+int(N/m) restraints (def=8)\n";
    print " -e e: exclude contacts with contact order smaller than e (def=10)";
    print "\n";
    exit(0);
}
