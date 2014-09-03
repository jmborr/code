#!/usr/bin/perl
use Getopt::Std;

&parse_input();
$sd=substr($h,1,1);
#print "h=$h\n";print "c=$c\n";print "out=$out\n"; 
chomp($currd=`pwd`);
$d="/tmp/jose/tasser_" . rand; `mkdir -p $d`; chdir "$d";
`cp /library/pdb_jul05/input/$h seq.txt`;
`cp /library/pdb_jul05/seq/$h.SEQ seq.dat`;
`cp /library/pdb_jul05/CA/$h.pdb CA`;
`cp $c common.tbz2`; `tar jxf common.tbz2`;
exit(0);
`./create_rmsinp.pl -i $h -o . -s seq.txt`;
print "INVOKING MKDAT\n"; `./mkdat.pl -d .`;
print "INVOKING MKPAIR\n"; `./mkpair.pl -d .`;
`echo 0 > comb.dat`; `echo 0 > combCA.dat`; #remove contact info
`echo 0 > dist.dat`; `echo 0 > distL.dat`;  #remove distance-restraints info
`echo "  0 hard" > chain.dat`;              #remove template info
`./cas_rd`;
`echo 16 > tra.in`; `ls -1 rep*.tra >> tra.in`;
`./spickerZZ.x -l tra.in -s seq.dat -a CA -m $h`;

`tar jcf casout.tbz2 out.d rep*.tra tra.in`;
`tar jcf spickout.tbz2 centroids.pdb closests.pdb summary.txt`;
`tar jcf inpout.tbz2 CA chain.dat *.dat pair.* rmsinp`;
`mkdir -p $out`;
`mv casout.tbz2 spickout.tbz2 inpout.tbz2 $out`;
`rm -rff $d`;

sub parse_input(){
    my $f=0;
    my ($currd,$dir,$fil);
    chomp($currd=`pwd`);
    if(!getopts('i:c:o:h')){ &message(); }
    if($opt_i){ $h=$opt_i; $f++; }
    if($opt_c){ $c =addAbsPath($opt_c); $f++; }
    if($opt_o){ $out=addAbsPath($opt_o); $f++; }
    if($f<3){ &message(); }
    if($opt_h){ &message(); }
}

#add absolute path to files or directories (even if inexistent)
sub addAbsPath(){
    my $c=shift;
    $c=~s/^\.\///; #remove leading './' if present
    unless($c=~/^\//){ chomp($dir=`pwd`); $c="$dir/$c"; } #add abs path
    while($c=~s/\/\w+\/\.\.//){;} #substitute '/*/..' with '/' one by one
    return $c;
}

sub message(){
    system("clear");
    print "Usage: mkjob.pl -i -c -o \n";
    print " Mandatory:\n";
    print " -i header (five letter header)\n";
    print " -c common.tbz2 (executables and common input files)\n";
    print " -o outdir (directory where to store the results)\n";
    print "\n";
    exit(0);
}

