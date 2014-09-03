#!/usr/bin/perl
use Math::Trig;

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

############# something you need to change #######################
#$run="real";
$run="benchmark";
$data_dir="/net/dell/02/users/yzhang/public/mkinput/examples/example_mkpair"; #for I/O
$bin_dir="/net/dell/02/users/yzhang/public/mkinput/bin"; #for bin
$comm_dir="/net/dell/02/users/yzhang/public/mkinput/data"; #for *.comm
$libdir="/library/yzhang/PDB";
$rundir="/tmp/yzhang/MKPAIR";

########### Working directory ##############################
`mkdir -p $rundir`;
`rm -fr $rundir/*`;
chdir "$rundir";

`/bin/cp $data_dir/seq.txt ./seq.txt`;
`/bin/cp $data_dir/seq.dat ./seq.dat`;
`/bin/cp $bin_dir/align ./align`;
`/bin/cp $bin_dir/pair65 $rundir/pair`;
`/bin/cp $comm_dir/matrix1.comm $rundir`;
`/bin/cp $comm_dir/matrix3.comm $rundir`;
`/bin/cp $comm_dir/blosum.comm  $rundir`;

############################################################
########### make MSA by PSIBLAST ###########################
############################################################
# make fasta sequence file ------------------>
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
    if($i==int($i/60)*60){
	printf seq "\n";
    }
}
printf seq "\n";
close(seq);

#run PSIBLAST ------------------------------------>
printf "running PSI-blast ............\n";
system("/library/yzhang/blast/blastpgp -b 5000 -j 3 -h 0.001 -d /library/yzhang/nr/nr.filter -i protein.seq -m 0 -o blast.out");
sleep(1);
`sync`;

#generate MSA file ------------------------------->
# read PSIBALST output ==>
printf "reading blast output file .........\n";
open(blast,"blast.out");
$i=0;
while($line=<blast>){
    ### decided starting point:--->
    if($line=~/Results from round 3/ || $line=~/^CONVERGED\!/){
	while($line=<blast>){
	    goto pos2 if($line=~/Results from round/);
	    if($line=~/Expect =\s*(\S+)/){
		$Evalue=$1;
		if($Evalue<1.0){
		    <blast>=~/Identities =\s*(\d+)\/(\d+)/;
		    if($2>0){
			$id=$1/$2;
		    }else{
			$id=0;
		    }
		    if($id<=0.9){
			$i++;
			for($j=1;$j<=$Lch;$j++){
			    $align{$i,$j}="-";
			}
			<blast>;
			while($line=<blast>){
			    if($line=~/^Query:\s+(\d+)\s+(\S+)\s+(\d+)/){
				$i1=$1;
				$se0=$2;
				$i2=$3;
				<blast>;
				<blast>=~/^Sbjct:\s+\d+\s+(\S+)/;
				$se=$1;
				$L_se=length $se0;
				$k=$i1-1;
				for($j=1;$j<=$L_se;$j++){
				    if(substr($se0,$j-1,1) ne "-"){
					$k++;
					$align{$i,$k}=substr($se,$j-1,1);
				    }
				}
				<blast>;
			    }else{
				goto pos1;
			    }
			}
		    }
		}
	    }
	    pos1:;
	}
    }
}
pos2:;
$N_MSA=$i;
# output MSA ==>
printf "output MSA files .....................\n";
open(msa,">msa.aln");
printf msa "%5d %5d protein\n",$N_MSA+1,$Lch;
for($j=1;$j<=$Lch;$j++){  # <---------- Query
    $a=substr($sequence,$j-1,1);
    printf msa "$a";
    $m=0;
    if($j==int($j/50)*50){
	printf msa "\*\n";
	$m=1;
    }
}
if($m==0){
    printf msa "\n";
}
for($i=1;$i<=$N_MSA;$i++){ # <-------- alignments
    for($j=1;$j<=$Lch;$j++){
	printf msa "$align{$i,$j}";
	$m=0;
	if($j==int($j/50)*50){
	    printf msa "\n";
	    $m=1;
	}
    }
    if($m==0){
	printf msa "\n";
    }
}
close(msa);

###########################################################
######### create pair3.dat ################################
###########################################################
### decided library list -------------------------------->
printf "removing homology templates >30% id ...........\n";
if($run eq "real"){
    `/bin/cp $libdir/list ./list`;
}else{
    open(listall,"$libdir/list");
    <listall>=~/(\d+)/;
    $n=$1; #total number of templates
    $k=0;
    for($i=1;$i<=$n;$i++){
	<listall>=~/(\S+)/;
	$p=$1;
	$ali=`./align protein.seq $libdir/$p\.pdb 2`;
	$ali=~/Identical length\:\s*(\d+)/;
	$L_id=$1;
	$id=$L_id/$Lch;
	printf "$p $id $i\n";
	if($id<0.35){
	    #printf "$p $id $k ---------------\n";
	    $k++;
	    $pp{$k}=$p;
	}
    }
    close(listall);
    $K=$k;
    open(list,">list");
    printf list "$K\n";
    for($i=1;$i<=$K;$i++){
	printf list "$pp{$i}\n";
    }
    close(list);
}

sleep(2);

# run pair ------------------------------------------->
printf "Running pair ................\n";
`./pair 0.3 $libdir/`; #please do not change 0.3
sleep(5);

# copy back the output files -------------------------->
`/bin/cp -f pair.3  $data_dir/pair.3`;
`/bin/cp -f pair.1  $data_dir/pair.1`;

chdir "/tmp/yzhang";
`rm -fr $rundir`;

sleep(5);
`sync`;
`sync`;

exit();
