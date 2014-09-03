#!/usr/local/bin/perl
# usage: >mkseq.pl

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

$bindir="/net/dell/02/users/yzhang/public/mkinput/bin";

#########make fasta input from 'seq.txt' #########################
open(fasta,">protein.fasta");
printf fasta "> protein\n";
open(seqtxt,"seq.txt");
$i=0;
$j=0;
while($line=<seqtxt>){
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
close(seqtxt);
close(fasta);

########### run PSIPRED ####################
`$bindir/psipred/runpsipred protein.fasta`;
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


exit();


