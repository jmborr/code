#!/usr/bin/perl

use Getopt::Std;
use small_util qw(addAbsPath pastry $codedir);

sub message(){
    #system("clear");
    print "Usage: mkdat.pl -d -o [-b -c -l -t -w]\n";
    print " Mandatory:\n";
    print " -d data_dir (where to find seq.fasta, seq.dat, template.dat and contact.dat. Default=current dir)\n";
    print " -o output directory (def='.')\n";
    print " Optional:\n";
    print " -b bin_dir  (def=$codedir/perl/tasser/tasser1.0/mkinput/bin)\n";
    print " -c comm_dir (def=$codedir/perl/tasser/tasser1.0/mkinput/data)\n";
    print " -l libdir   (default=/library/yzhang/PDB)\n";
    print " -t database (default=/library/yzhang/nr/nr.filter)\n";
    print " -w work_dir (default=./tmp_mkdat)\n";
    print " -v sequence file in FASTA format, located in data_dir (def='seq.fasta')\n";
    print " -x secondary structure propensity, located in data_dir (def='seq.dat')\n";
    print " -y template file, located in data_dir (def='template.dat')\n";
    print " -z predicted contacts, located in data_dir (def='contact.dat')\n";
    print "\n";
    exit(0);
}

sub parse_input(){
    my $f=0;
    if(!getopts('b:c:d:l:w:t:v:x:y:z:o:h')){ &message(); }
    if($opt_d){ $data_dir=addAbsPath($opt_d); $f++; }
    if($opt_b){ $bin_dir =addAbsPath($opt_b);       }
    if($opt_c){ $comm_dir=addAbsPath($opt_c);       }
    if($opt_l){ $libdir  =addAbsPath($opt_l);       }
    if($opt_w){ $work_dir=addAbsPath($opt_w);       }
    if($opt_t){ $db      =addAbsPath($opt_t);       }
    if($opt_v){ $seqfasta=$opt_v;                   }
    if($opt_x){ $seqdat  =$opt_x;                   }
    if($opt_y){ $tpldat  =$opt_y;                   }
    if($opt_z){ $cntdat  =$opt_z;                   }
    if($opt_o){ $outd    =addAbsPath($opt_o); $f++; }
    if($f<2){ print "Hello\n";  &message(); }
    if($opt_h){ &message(); }
}

print "...running mkdat.pl...\n";

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

     'ASX'=>'B',
     'GLX'=>'Z',
     'UNK'=>'X',

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

     'B'=>'ASX',
     'Z'=>'GLX',
     'X'=>'CYS',
     );

@ww=qw(
       0.05
       0.10
       0.15
       0.20
       0.25
       0.30
       0.35
       0.40
       0.45
       0.50
       0.55
       0.60
       0.65
       0.70
       0.75
       0.80
       0.85
       );

@AA=qw(
       C
       M
       F
       I
       L
       V
       W
       Y
       A
       G
       T
       S
       Q
       N
       E
       D
       H
       R
       K
       P
       );

################# default directories and files #############################
$data_dir=addAbsPath(".");
$bin_dir="$codedir/perl/tasser/tasser1.0/mkinput/bin";
# comm_dir for *.comm
$comm_dir="$codedir/perl/tasser/tasser1.0/mkinput/data";
$work_dir=addAbsPath("./tmp_mkdat");
$libdir="/library/yzhang/PDB/";
$blastdir="/library/yzhang/blast";
$db="/library/yzhang/nr/nr.filter";
$seqfasta='seq.fasta';
$seqdat='seq.dat';
$tpldat='template.dat';
$cntdat='contact.dat';
# customary files and directories
&parse_input();
$libdir.='/' unless $libdir=~/\/$/ ; #libdir must end in '/'

if($seqfasta=~/\// || $seqdat=~/\// || $tpldat=~/\// || $cntdat=~/\//){
    print "ERROR flags v,x,y,z don't support '/'\n";
}
################ working directory ########################
`/bin/mkdir -p $work_dir`;
chdir "$work_dir"; 
unless($work_dir=~/^\//){ chomp($work_dir=`pwd`); }
#`/bin/rm -f $work_dir/*`;
`/bin/cp $data_dir/$seqfasta ./seq.fasta`;
`/bin/cp $data_dir/$seqdat ./seq.dat`;
`/bin/cp $data_dir/$tpldat ./template.dat`;
`/bin/cp $data_dir/$cntdat ./contact.dat`;
`/bin/cp $comm_dir/wgt.tar.bz2 .`;
`/usr/bin/bunzip2 wgt.tar.bz2`;
`/bin/tar -xvf wgt.tar`;

`/bin/cp $bin_dir/zal ./zal`;
`/bin/cp $bin_dir/dat ./dat`;
`/bin/cp $bin_dir/solve ./solve`;

################ make 'protein.seq' #################
@seqfastas=`cat seq.fasta`;
$sequence="";
foreach $seqfasta(@seqfastas){
    goto pos6 if($seqfasta=~/\>/);
    $seqfasta=~s/\s//mg;
    $seqfasta=~s/\n//mg;
    $sequence=$sequence.$seqfasta;
  pos6:;
}
$Lch=length $sequence;
open(seq,">protein.seq");
printf seq ">protein\n";
for($i=1;$i<=$Lch;$i++){
    $a=substr($sequence,$i-1,1);
    printf seq "$a";
    $seqQ{$i}=$a;
    $log{$i,$seqQ{$i}}++;
    if($i==int($i/60)*60){
	printf seq "\n";
    }
}
printf seq "\n";
close(seq);

############## make exp.dat #########################
print "...mkdat.pl calling blastpgp...\n...running blastpgp...\n"; 
system("$blastdir/blastpgp -i protein.seq -Q protein.mat3 -d $db -o protein.out -e 1e-3 -h 1e-3 -M BLOSUM62 -j 3 -m 6");
print "...finished blastpgp...\n";

print "...mkdat.pl calling solve...\n...runing solve...\n";
foreach $w(@ww){
    `./solve weight.$w protein`; #protein.neu
    open(neu,"protein.neu");
    <neu>=~/(\d+)/;
    $Lch=$1;
    for($i=1;$i<=$Lch;$i++){
	<neu>=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	$exp{$i,$w}=$5;
    }
    close(neu);
}
print "...finished solve...\n";
open(exp1,">exp.dat");
printf exp1 "$Lch  ";
foreach $w(@ww){
    printf exp1 " %5s",$w;
}
printf exp1 "\n";
for($i=1;$i<=$Lch;$i++){
    printf exp1 "%5d",$i;
    foreach $w(@ww){
	printf exp1 " %5d",$exp{$i,$w};
    }
    printf exp1 "\n";
}
close(exp1);

########### run psi-blast #######################
print "...mkdat.pl calling blastpgp...\n...running blastpgp...\n"; 
system("$blastdir/blastpgp  -b 1000 -j 3 -h 0.001 -d $db -i protein.seq -C psitmp.chk > blast.out");
print "...finished blastpgp...\n";

########### calculate 'protein.frq' ###################
open(blast,"blast.out");
while($line=<blast>){
    if($line=~/Results from round\s+(\d+)/){
	$ROUND=$1;
    }
}
seek(blast,0,0);
$it=0;
while($line=<blast>){
    if($line=~/round\s+$ROUND/){
	while($line=<blast>){
	    if($line=~/Expect =\s*(\S+)/){
		$ev=$1;
		<blast>=~/Identities =\s*\S+\s+\((\S+)\%/;
		<blast>;
		$id=$1;
		$it++;
		if($ev<0.001 && $id < 98){
		    while($line=<blast>){
			if($line=~/Query\:\s*(\d+)\s+(\S+)\s+(\S+)/){
			    $i1=$1;
			    $seq1=$2;
			    <blast>;
			    <blast>=~/Sbjct\:\s*(\S+)\s+(\S+)\s+(\S+)/;
			    $seq2=$2;
			    <blast>;
			    ###
			    $L=length $seq1;
			    $m1=$i1-1;
			    for($i=1;$i<=$L;$i++){
				$q1=substr($seq1,$i-1,1);
				$q2=substr($seq2,$i-1,1);
				$m1++ if($q1 ne '-');
				if($q1 ne '-' && $q2 ne '-'){
				    $log{$m1,$q2}++;
				}
			    }
			    ###
			}else{
			    goto pos1;
			}
		    }
		}
	      pos1:;
	    }
	}
    }
}
close(blast);
open(freq,">protein.frq");
printf freq "$Lch\n";
for($i=1;$i<=$Lch;$i++){
    printf freq "%3d",$i;
    $norm=0;
    foreach $A(@AA){
	$norm+=$log{$i,$A};
    }
    foreach $A(@AA){
	printf freq "%10.7f",$log{$i,$A}/$norm;
    }
    printf freq "\n";
}
close(freq);

########### make in.dd #############
open(in,">in.dd");
printf in "protein.seq\n";
printf in "protein.frq\n";
printf in "$libdir\n";
printf in "seq.pdb\n";
printf in "seq.dat\n";
printf in "0.3\n"; #sequence cutoff
printf in "template.dat\n";
printf in "contact.dat\n";
printf in "8   15.5\n"; #zscore cutoff for medium and easy classes for PROSPECTOR_3 templates. This two numbers are adjustable when PROSPECTOR_3 is updated.
close(in);

########### run dat ################
print "...mkdat.pl calling zal...\n...running zal...\n"; 
`./zal`;
print "...finished zal...\n";

print "...mkdat.pl calling dat...\n...running dat...\n"; 
`./dat`;
print "...finished dat...\n";

`/bin/mkdir -p $outd`;

@listOut=("comb.dat","combCA.dat","dist.dat","distL.dat","chain.dat","par.dat","exp.dat");
for $file (@listOut){
    if(system("ls $file >/dev/null 2>/dev/null")){
	print "Error: file $file does not exists\n";
	exit(1);
    }
    `/bin/cp $file  $outd`;
}

#`/bin/rm -rf $work_dir`;
`sync`;
`sync`;
sleep(2);
print "...finished mkdat.pl...\n";

exit(0);
