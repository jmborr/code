#!/usr/bin/perl

use Getopt::Std;
use small_util qw(addAbsPath codedir db);

$currd=`pwd`; chomp($currd);
$outd=$currd;
$datab=$db/pbdNonHom.fasta;
$E=10;
$seqIdCutOff=0.35;

&parse_input();

if(! -e $outd){ system("mkdir -p $outd"); } #creates non-existent directory
if(!defined($outf)){ $outf="$header.homol"; }
$outf="$outd/$outf"; #add absolute path

#create temporary working directory and work there
$workd="$currd/tmp_homol".int(rand(1000000));
system("mkdir -p $workd");
chdir($workd);

#get/create fasta file
if(!defined(fastaf)){
    $cmd="$codedir/python/seq/mkfasta.py -i $header -o $header.fasta";
    if(defined(rawf)){$cmd="$cmd  -s $raw"; }
    system(cmd);
}

#generate alignment
`$codedir/bin/msa/fasta $fastaf $datab -b 2000 -E $E -Q -m 10 -p -w 50 -O $header.out >/dev/null`;

@l=`egrep ">>|sw_ident" $header.out`; 
shift @l; #first line is useless
$buff="";
$n=0;
for($i=0;$i<10;$i++){
    $ident=substr($l[$i+1],12,4);
    if($ident>seqIDCutOff){ $n++; $buff= $buff.substr($l[$i],2,5)."\n"; }
}
open(out,">$outf");   print out "$n\n$buff";   close(out);

chdir($currd); system("/bin/rm -rf $workd");

sub parse_input(){
    my $f=0;
    if(!getopts('i:f:s:o:p:d:e:d:h')){ &message(); }
    if($opt_i){ $header =$opt_i; $f++; }
    if($opt_f){ $fastaf =addAbsPath($opt_f); }
    if($opt_s){ $rawf =addAbsPath($opt_f); }
    if($opt_o){ $outd =addAbsPath($opt_o); }
    if($opt_p){ $outf =$opt_p; }
    if($opt_d){ $datab =addAbsPath($opt_d); }
    if($opt_e){ $E =$opt_e; }
    if($opt_d){ $ =$opt_d; }
    if($opt_h){ &message(); }
    if($f<1){ $message(); }
}

sub message(){
    system("clear");
    print "Usage:find_homologs.pl \n";
    print " Required:\n";
    print " -i header : file-letter (xxxxx) pdb code\n";
    print " Optional:\n";
    print " -f fastaf : fasta file\n";
    print " -s rawf   : raw sequence file\n";
    print " -o outd   : output directory (default=current directory)\n";
    print " -p outf   : output filename (default=${header}.homol)\n";
    print " -d datab  : database of sequences to compare (default=$db/pdbNonHom.fasta)\n";
    print " -e E      : Evalue (default=10)\n";
    print " -d ident  : sequence identity threshold (default=0.35)\n";
    print "\n";
    exit(0);
}
