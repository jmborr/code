#!/usr/bin/perl
use Getopt::Std;
use small_util qw(addAbsPath $codedir);

print "...running mkseq.pl...\n";

$rootdir="$codedir/perl/tasser/tasser1.0/mkinput";
$seqraw='seq.raw';

&parse_input(); #print "rootdir=$rootdir seqraw=$seqraw\n";
`/bin/cp $seqraw seq.txt 2> /dev/null`;
########### setup  the environment and Working DIRectory ###
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/usr/pgi/linux86/bin";
$ENV{'LD_LIBRARY_PATH'}="/usr/local/lib:/usr/lib:/lib";

%ts=(
     'GLY'=>'G',
     'ALA'=>'A',
     'VAL'=>'V',
     'LEU'=>'L',
     'ILE'=>'I',
     'SER'=>'S',
     'THR'=>'T',
     'CYS'=>'C',
     'MET'=>'M',
     'PRO'=>'P',
     'ASP'=>'D',
     'ASN'=>'N',
     'GLU'=>'E',
     'GLN'=>'Q',
     'LYS'=>'K',
     'ARG'=>'R',
     'HIS'=>'H',
     'PHE'=>'F',
     'TYR'=>'Y',
     'TRP'=>'W',
     'G'=>'GLY',
     'A'=>'ALA',
     'V'=>'VAL',
     'L'=>'LEU',
     'I'=>'ILE',
     'S'=>'SER',
     'T'=>'THR',
     'C'=>'CYS',
     'M'=>'MET',
     'P'=>'PRO',
     'D'=>'ASP',
     'N'=>'ASN',
     'E'=>'GLU',
     'Q'=>'GLN',
     'K'=>'LYS',
     'R'=>'ARG',
     'H'=>'HIS',
     'F'=>'PHE',
     'Y'=>'TYR',
     'W'=>'TRP',
     'a'=>'CYS',
     'b'=>'CYS',
     'c'=>'CYS',
     'd'=>'CYS',
     'e'=>'CYS',
     'f'=>'CYS',
     'g'=>'CYS',
     'h'=>'CYS',
     'i'=>'CYS',
     'j'=>'CYS',
     'k'=>'CYS',
     'l'=>'CYS',
     'm'=>'CYS',
     'n'=>'CYS',
     'o'=>'CYS',
     'p'=>'CYS',
     'q'=>'CYS',
     'r'=>'CYS',
     's'=>'CYS',
     't'=>'CYS',
     'u'=>'CYS',
     'v'=>'CYS',
     'w'=>'CYS',
     'x'=>'CYS',
     'y'=>'CYS',
     'z'=>'CYS',
     'X'=>'GLY',
    );

$bindir="$rootdir/bin";
if( ! -e $bindir ){
    print "ERROR: no $bindir directory found !\n";
    exit(0);
}

#########make fasta input from 'seq.raw' #########################
open(fasta,">protein.fasta");
printf fasta "> protein\n";
open(seqraw,"seq.raw");
$i=0;
$j=0;
while($line=<seqraw>){
    goto pos1 if($line=~/^>/);
    if($line=~/(\S+)/){
	$sequence=$1;
	$Lch=length $sequence;
	for($k=1;$k<=$Lch;$k++){
	    $i++;
	    $j++;
	    $seq1=substr($sequence,$k-1,1);
	    $seq{$j}=$ts{$seq1};
	    printf fasta "$seq1";
	    if($i==60){
		printf fasta "\n";
		$i=0;
	    }
	}
    }
  pos1:;
}
if($i != 60){
    printf fasta "\n";
}
close(seqraw);
close(fasta);

########### run PSIPRED ####################
#print "$bindir/psipred/runpsipred protein.fasta $bindir/psipred";
`$bindir/psipred/runpsipred protein.fasta $bindir/psipred`;
sleep(1);
########### make 'seq.dat' ########################
open(psipred,"protein.horiz");
open(yan,">seq.dat");
$j=0;
while($line=<psipred>){
    if($line=~/Conf:\s+(\d+)/){
	$conf=$1;
	<psipred>=~/Pred:\s+(\S+)/;
	$pred=$1;
	<psipred>=~/AA:\s+(\S+)/;
	$aa=$1;
	$num=length $aa;
	for($i=1;$i<=$num;$i++){
	    $j++;
	    $conf1=substr($conf,$i-1,1);
	    $pred1=substr($pred,$i-1,1);
	    $aa1=substr($aa,$i-1,1);
	    $sec{$j}=1;
	    $sec{$j}=2 if($conf1 >=1 && $pred1 eq 'H');
	    $sec{$j}=4 if($conf1 >=1 && $pred1 eq 'E');
	    printf yan "%5d   %3s%5d%5d\n",$j,$seq{$j},$sec{$j},$conf1;
	}
    }
}
close(yan);
close(psipred);

print "...finished mkseq.pl...\n";

exit();

sub parse_input(){
    my $f=0;
    if(!getopts('r:s:h')){ &message(); }
    if($opt_r){ $rootdir =addAbsPath($opt_r);       }
    if($opt_s){ $seqtxt  =addAbsPath($opt_s);       }
    if($f<0){   &message(); }
    if($opt_h){ &message(); }
}

sub message(){
    system("clear");
    print "Usage: mkseq.pl [-r -s]\n";
    print " Optional:\n";
    print " -r rootdir (def=$codedir/perl/tasser/tasser1.0/mkinput)\n";
    print " -s raw seq file (def='seq.raw')\n";
    print "\n";
    exit(0);
}
