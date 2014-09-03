#!/usr/bin/perl

########################################################
# Given a directory ($lib) with pdb's, mtx's and cnt2's files, the program
# will create file seq.pdb with all sequences in fasta-like format
########################################################

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

$lib="/gpfs1/active/jose/projects/updateMkLib/PDBextras";

############## mkdir list ########################
@lines=`cd $lib; /bin/ls |more`;
$n=0;
foreach $line(@lines){
    if($line=~/(\S+)\.mtx/){
	$n++;
	$s{$n}=$1;
    }
}

open(list1,">list");
printf list1 "$n\n";
for($i=1;$i<=$n;$i++){
    printf list1 "$s{$i}\n";
}
close(list1);

############# make seq.pdb #########################
open(seqpdb,">seq.pdb");
for($j=1;$j<=$n;$j++){
    $pdb="$lib/$s{$j}\.pdb";
    if(!-s "$pdb"){
	printf "$s{$j} without $lib/$s{$j}\.pdb!!!\n";
	exit();
    }
    if(!-s "$lib/$s{$j}\.cnt2"){
	printf "$s{$j} without $lib/$s{$j}\.cnt2!!!\n";
	exit();
    }
    printf seqpdb "> $s{$j}\n";
    $i=0;
    open(str,"$pdb");
    while($line=<str>){
	if(substr($line,12,4)=~/CA/){
	    $i++;
	    $seq=$ts{substr($line,17,3)};
	    printf seqpdb "$seq";
	    $mk=1;
	    if(int($i/60)*60 == $i){
		printf seqpdb "\n";
		$mk=-1;
	    }
	}
    }
    close(str);
    if($mk == 1){
	printf seqpdb "\n";
    }
    printf seqpdb "\n";
}
close(seqpdb);

exit();
