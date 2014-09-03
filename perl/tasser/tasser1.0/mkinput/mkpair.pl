#!/usr/bin/perl
use Math::Trig;
use Getopt::Std;
use small_util qw(addAbsPath $codedir $pdbnonhom);

print "...running mkpair.pl...\n";

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

$run="benchmark"; #$run="real";
$bin_dir="$codedir/perl/tasser/tasser1.0/mkinput/bin"; #for bin
$comm_dir="$codedir/perl/tasser/tasser1.0/mkinput/data"; #for *.comm
$libdir="/library/yzhang/PDB";
$rundir="./tmp_mkpair".int(rand(1000000));

&parse_input();

########### Working directory ##############################
`mkdir -p $rundir`;
`rm -fr $rundir/*`;
chdir "$rundir";
if(defined($data_dir)){
    `/bin/cp $data_dir/seq.fasta seq.fasta`;
    `/bin/cp $data_dir/seq.dat seq.dat`;
}
else{
    if(defined($fastaf)){
	`/bin/cp $fastaf seq.fasta`; 
    }  
    elsif(defined($header)){
	$fastaf="$pdbnonhom/input/$header"; #print "$fastaf\n";exit(0);	
	if( -e $fastaf){ 
	    `/bin/cp $fastaf seq.fasta`; #print "fastaf done!\n";
	}
	else{
	    `$codedir/python/seq/mkfasta.py -i $header`;
	}
    }
    else{
	print "ERROR: no means to obtain the fasta file\n"; exit(1);
    }

    if(defined($datf)){
	`/bin/cp $datf seq.dat`; 
    }  
    elsif(defined($header)){
	$datf="$pdbnonhom/seq/$header.SEQ"; #print "$datf\n";exit(0);	
	if( -e $datf){ 
	    `/bin/cp $datf seq.dat`; #print "datf done!\n"; exit(0);
	}
	else{
	    #invoque mkdat.py 
	    print "ERROR: no means to obtain the sec.struc.predict file\n";
	    exit(0);
	}
    }
    else{
	print "ERROR: no means to obtain the sec.struc.predict file\n";
	exit(0);
    }
    
}
`/bin/cp $bin_dir/align align`;
`/bin/cp $bin_dir/pair65 pair`;
`/bin/cp $comm_dir/matrix1.comm .`;
`/bin/cp $comm_dir/matrix3.comm .`;
`/bin/cp $comm_dir/blosum.comm  .`;

############################################################
########### make MSA by PSIBLAST ###########################
############################################################
# make fasta sequence file ------------------>
@seqtxts=`cat seq.fasta`;
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
print "...mkpair.pl calling blastpgp...\n...running blastpgp...\n";
system("/library/yzhang/blast/blastpgp -b 5000 -j 3 -h 0.001 -d /library/yzhang/nr/nr.filter -i protein.seq -m 0 -o blast.out >/dev/null");
wait;
sleep(1);
`sync`;
print "...finished blastpgp...\n";

#generate MSA file ------------------------------->
# read PSIBALST output ==>
#printf "reading blast output file .........\n";
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
		    if($id<=0.9){ #only accept identities smaller than 90%
			$i++;
			for($j=1;$j<=$Lch;$j++){
			    $align{$i,$j}="-"; #initialize $align{$i,0..$Lch}
			}
			#read now the particular alignment
			<blast>;
			while($line=<blast>){
			    if($line=~/^Query:\s+(\d+)\s+(\S+)\s+(\d+)/){
				$i1=$1;
				$se0=$2;
				$i2=$3;
				<blast>; #don't read the line of matches and similarities(+)
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
				<blast>; #there's a blank line between two segment-alignments
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
#printf "output MSA files .....................\n";
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
#printf "removing homology templates >30% id ...........\n";
if($run eq "real"){
    `/bin/cp $libdir/list ./list`;
}
else{
    print "...mkpair.pl calling align...\n...running align...\n";
    open(listall,"$libdir/list");
    <listall>=~/(\d+)/;
    $n=$1; #total number of templates
    $k=0;
    for($i=1;$i<=$n;$i++){
	<listall>=~/(\S+)/;
	$p=$1;
	$ali=`./align protein.seq $libdir/$p\.pdb 2     >/dev/null`;
	$ali=~/Identical length\:\s*(\d+)/;
	$L_id=$1;
	$id=$L_id/$Lch;
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
    print "...finished align...\n"
}

sleep(2);

# run pair ------------------------------------------->
print "...mkpair.pl calling pair...\n...running pair...\n";
`./pair 0.3 $libdir/ >/dev/null`; #please do not change 0.3
sleep(5);
print "...finished pair...\n";

# copy back the output files. I'm unsure of pair output
`/bin/mv -f pair.3  $out_dir/pair3.dat 2>/dev/null`;
`/bin/mv -f pair.1  $out_dir/pair1.dat 2>/dev/null`;
`/bin/mv -f pair3.dat  $out_dir        2>/dev/null`;
`/bin/mv -f pair1.dat  $out_dir        2>/dev/null`;


chdir "../";
`rm -r -f $rundir`;

sleep(5);
`sync`;
`sync`;

print "...finished mkpair.pl...\n";

exit(0);

sub parse_input(){
    my $f=0;
    if(!getopts('a:b:c:d:e:i:l:r:s:o:h')){ &message(); }
    if($opt_o){ $out_dir =addAbsPath($opt_o); $f++; }
    if($opt_i){ $header  =$opt_i;                   }
    if($opt_a){ $fastaf  =addAbsPath($opt_a);       }
    if($opt_e){ $datf    =addAbsPath($opt_e);       }
    if($opt_d){ $data_dir=addAbsPath($opt_d);       }
    if($opt_b){ $bin_dir =addAbsPath($opt_b);       }
    if($opt_c){ $comm_dir=addAbsPath($opt_c);       }
    if($opt_l){ $libdir  =addAbsPath($opt_l);       }
    if($opt_r){ $rundir  =addAbsPath($opt_r);       }
    if($opt_s){ $run     =addAbsPath($opt_s);       }
    if($f<1){ &message(); }
    if($opt_h){ &message(); }
}

sub message(){
    system("clear");
    print "Usage: mkpair.pl [options]\n";
    print " Mandatory:\n";
    print " -o out_dir (where to store pair.1 and pair.3)\n";
    print " Optional:\n";
    print " -i header five letter pdb code\n";
    print " -a seq.fasta (default=$pdbnonhom/input/xxxxx)\n";
    print " -e seq.dat   (default=$pdbnonhom/input/xxxxx.SEQ)\n";
    print " -d data_dir (where to find seq.fasta and seq.dat)\n";
    print " -b bin_dir (def=$codedir/perl/tasser/tasser1.0/mkinput/bin)\n";
    print " -c comm_dir (def=$codedir/perl/tasser/tasser1.0/mkinput/data)\n";
    print " -l libdir   (default=/library/yzhang/PDB)\n";
    print " -r rundir   (temporary work dir, default=./tmp_mkpair)\n";
    print " -s benchmark|real (default=benchmark)\n";
    print "\n";
    exit(0);
}
