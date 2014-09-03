#!/usr/bin/perl

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

################# directories #############################
$data_dir="/net/dell/02/users/yzhang/public/mkinput/examples/example_mkdat";
$bin_dir="/net/dell/02/users/yzhang/public/mkinput/bin";
$comm_dir="/net/dell/02/users/yzhang/public/mkinput/data"; #for *.comm
$work_dir="/tmp/yzhang/MKINPUT";
$libdir="/library/yzhang/PDB/";
$blastdir="/library/yzhang/blast";
$db="/library/yzhang/nr/nr.filter";

################ working directory ########################
`/bin/mkdir -p $work_dir`;
chdir "$work_dir";
`/bin/rm -f $work_dir/*`;
`/bin/cp $data_dir/seq.txt ./seq.txt`;
`/bin/cp $data_dir/seq.dat ./seq.dat`;
`/bin/cp $data_dir/template.dat ./template.dat`;
`/bin/cp $data_dir/contact.dat ./contact.dat`;
`/bin/cp $comm_dir/wgt.tar.bz2 .`;
`/usr/bin/bunzip2 wgt.tar.bz2`;
`/bin/tar -xvf wgt.tar`;

`/bin/cp $bin_dir/zal ./zal`;
`/bin/cp $bin_dir/dat ./dat`;
`/bin/cp $bin_dir/solve ./solve`;

################ make 'protein.seq' #################
@seqtxts=`cat seq.txt`;
$sequence="";
foreach $seqtxt(@seqtxts){
    goto pos6 if($seqtxt=~/\>/);
    $seqtxt=~s/\s//mg;
    $seqtxt=~s/\n//mg;
    $sequence=$sequence.$seqtxt;
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
system("$blastdir/blastpgp -i protein.seq -Q protein.mat3 -d $db -o protein.out -e 1e-3 -h 1e-3 -M BLOSUM62 -j 3 -m 6");
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
`$blastdir/blastpgp  -b 1000 -j 3 -h 0.001 -d $db -i protein.seq -C psitmp.chk > blast.out`;

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
`./zal`;
`./dat`;

`/bin/cp comb.dat    $data_dir`;
`/bin/cp combCA.dat  $data_dir`;
`/bin/cp dist.dat    $data_dir`;
`/bin/cp distL.dat   $data_dir`;
`/bin/cp chain.dat   $data_dir`;
`/bin/cp par.dat     $data_dir`;
`/bin/cp exp.dat     $data_dir`;

chdir "/tmp/yzhang";
`rm -fr $work_dir`;

`sync`;
`sync`;
sleep(2);

exit();
