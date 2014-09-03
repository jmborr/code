#!/usr/bin/perl

use small_util qw(addAbsPath $codedir);
use Getopt::Std;

undef($header);
undef($outf);
undef($classf);
&parse_input(); #may define some of the previous undefs

%class=(
	'a'=>'alph',  #all alpha                  protein
	'b'=>'beta',  #all beta                     "
	'c'=>' a/b',  #alpha/beta                   "
	'd'=>' a+b',  #alpha and beta               "
	'e'=>'mult',  #multidomainn                 "
	'f'=>'memb',  #membrane and cell surface    "
	'g'=>'smal',  #samll                        "
	'h'=>'coil',  #coiled-coil                  "
	'i'=>'lowr',  #low-resolution structure     "
	'j'=>'pept',  #peptide                      "
	'k'=>'desg'   #designed                     "
 # if undefined, then 'modl', for theoretical model "
	);

#store in memory the protein class for each header
#typical line:
#d1vwtc_	1vwt	C:	a.1.1.2	15283	cl=46456,...some more stuff
open(in, "$codedir/perl/pdb_utils/dir.cla.scop.txt_1.67");
$n=0;
while($l=<in>){
    unless(substr($l,0,1) eq '#'){
	$n++;
	@list=split(' ',$l); 
	$h4=$list[1]; #four letter pdb code
	if($list[1]=~/__/){ $ids='_'; }
	else{ $list[2]=~s/[^A-Z]//g; $ids=$list[2];} #$ids contain chain ID's
	$type=substr($list[3],0,1); #in the example, 'a' of 'a.1.1.2'
	for $id (split(//,$ids)){
	    $h=$h4.$id; #five-letter code
            #dictonary maps headers to one-letter code protein-class types
	    unless(defined($x{$h})){ $x{$h}=$type; }
	}
    }
}
close(in);
#for $key (keys %x){ print "$key $x{$key} $class{$x{$key}}\n";} exit(0);

#did we pass only one header, or a list of them?
if(defined($header)){ @list=($header); }
else{ @list=`cat $listh`; for $h (@list){ chomp($h); } }
#$,="\n";print @list;

#if we didn't pass output file, or if we passed only one header, print results to STDOUT
if(defined($outf)){ open(out,">$outf"); }
else{ open(out,">&STDOUT"); } #print to STDOUT

#go through the list of headers
for $h (@list){
    unless( defined($class{$x{$h}}) ){ $type="modl"; }
    else{ $type=$class{$x{$h}}; }
    if( !defined($classf) ||
	( defined($classf) && ($classf eq $type) ) ){
	print out "$h  $type\n";
    }
}
sub parse_input(){
    my $f=0;
    if(!getopts('o:i:l:f:a:h')){ &message(); }
    if($opt_o){ $outf=addAbsPath($opt_o); }
    if($opt_i){ $header=$opt_i; }
    if($opt_l){ $listh=addAbsPath($opt_l); }
    if($opt_f){ $classf=$opt_f; }
    if($opt_a){ $a=1; }
    #if($f<1){ &message(); } 
    if($opt_h){ &message(); }
}

sub message(){
    system("clear");
    print "Usage: class_find.pl -i -l -o\n";
    print " -o outf: output file where to write the results. If not provided, or if -i flag is used, will output instead to STDOUT\n";
    print " -i header: a file-letter pdb-code\n";
    print " -l listh: list of headers. Use this option if you want to run more than one header\n";
    print " -f classf: filter the list of headers, outputing only those headers of four-letter code for protein class (available classes are 'alph','beta',' a/b',' a+b','mult','memb','smal','coil','lowr','pept','desg'). Don't forget to enclose the code in quotes!\n";
    print " -a will output the list of headers plus and extra column with the class of each header. Use in place of -f flag\n";
    print "\n";
    exit(0);
}
