#!/usr/bin/perl

# usage: mkmtx.pl 1p0zA
##############################################################
# This program will generate 1p0zA.pdb, 1p0zA.mtx, 1p0zA.cnt2
# from /library/pdb/pdb1p0z.ent.
##############################################################


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

     'ASX'=>'A',
     'GLX'=>'G',
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

@AA=qw(
       GLY
       ALA
       VAL
       LEU
       ILE
       SER
       THR
       CYS
       MET
       PRO
       ASP
       ASN
       GLU
       GLN
       LYS
       ARG
       HIS
       PHE
       TYR
       TRP
       );

$s=$ARGV[0];

###############################################################
$pdblib="/library/pdb";
$blastdir="/library/yzhang/blast";
$dbname="/library/yzhang/nr/nr.filter";

#1##########################################################
############# generate '1xxxA.pdb' #########################
############################################################
$pdbname=substr($s,0,4);
$ch=substr($s,4,1);
$ch=" " if($ch eq "_");
open(pdb,"$pdblib/pdb$pdbname\.ent");
$i_atom=0;
$i_res=0;
$num_old="xxxxxxxxx";
while($line=<pdb>){
    $record=substr($line,0,6);
    $chid=substr($line,21,1);
    $alt=substr($line,16,1); #alternative, only first.
    if($record eq "ATOM  " && $chid eq $ch){
	#### get i_res:
	$num=substr($line,16,11);
	if($num ne $num_old){
	    $i_res++;
	    $num_old=$num;
	    ##### get seq:
	    $seq{$i_res}=substr($line,17,3);
	    foreach $A(@AA){
		goto pos1 if($A eq $seq{$i_res});
	    }
	    $seq{$i_res}="GLY";
	  pos1:;
	}
	#########check this atom appear before?
	if($alt ne " "){ 
	    $atom_tmp=substr($line,12,4);
	    for($i=1;$i<=$n_atom{$i_res};$i++){
		goto pos3 if($atom_tmp eq $atom{$i_res,$i});
	    }
	}
	$n_atom{$i_res}++;
	$atom{$i_res,$n_atom{$i_res}}=substr($line,12,4);
	$xyz{$i_res,$n_atom{$i_res}}=substr($line,30,50);
      pos3:;
    }
    if($i_res>1){
	if($record eq 'TER' ||$record eq 'ENDMDL'||$record eq 'MODEL'){
	    goto pos2;
	}
    }
}
 pos2:;  #read pdb finished
close(pdb);
### output ------------->
open(pdbyang,">$s\.pdb");
$k_atom=0;
$k_res=0;
for($i=1;$i<=$i_res;$i++){
    $mk=0;
    for($j=1;$j<=$n_atom{$i};$j++){
	if($atom{$i,$j}=~/CA/){
	    $k_res++;
	    $SEQ{$k_res}=$seq{$i};                 #for checking
	    $X{$k_res}=substr($xyz{$i,$j},0,8);    #for checking
	    $Y{$k_res}=substr($xyz{$i,$j},8,8);    #for checking
	    $Z{$k_res}=substr($xyz{$i,$j},16,8);   #for checking
	    $mk=1;
	}
    }
    if($mk == 1){
	for($j=1;$j<=$n_atom{$i};$j++){
	    $k_atom++;
	    printf pdbyang "%6s%5d %4s %3s  %4d    %50s\n",
	    "ATOM  ",$k_atom,$atom{$i,$j},$seq{$i},$k_res,$xyz{$i,$j};
	}
    }
}
printf pdbyang "END\n";
close(pdbyang);
####### if CA only --------------->
open(pdbdat,"$s\.pdb");
@lines=<pdbdat>;
close(pdbdat);
$n=@lines;
$n--;
if($n<=$k_res*4){
    printf "$n---------$s\n";
    `cp $s\.pdb tmp.pdb`;
    `./pulchra -vpc tmp.pdb`;
    `mv rebuilt_tmp.pdb $s\.pdb`;
}

#goto pos10;
#2##########################################################
############# generate '1xxxA.mtx' #########################
############################################################
######## make "seq.txt"------------>
open(pdb,"$s\.pdb");
$j=0;
while($line=<pdb>){
    if(substr($line,0,4) eq "ATOM" && substr($line,12,4)=~/CA/){
	$j++;
	$seq{$j}=$ts{substr($line,17,3)};
    }
}
close(pdb);
$Lch=$j;
open(seq,">$s\.seq");
for($j=1;$j<=$Lch;$j++){
    printf seq "$seq{$j}";
    if(int($j/70)*70==$j){
	printf seq "\n";
    }
}
printf seq "\n";
####### run psiblast -------------->
system("$blastdir/blastpgp -b 0 -j 3 -h 0.001 -d $dbname -i $s\.seq -C $s\.chk");
`echo $s\.chk > $s\.pn`;
`echo $s\.seq > $s\.sn`;
system("$blastdir/makemat -P $s");

pos10:;
#3##########################################################
############# generate '1xxxA.cnt2' ########################
############################################################
undef %n_atom;
open(pdb,"$s\.pdb");
################## read CA, and SG atoms ###########################
while($line=<pdb>){
    $record=substr($line,0,4);
    if($record eq "ATOM"){
	$atom=substr($line,12,4);
	$atom=~s/\s//g;
	$atom=~s/\d+//g;
	if($atom eq "CA"){    ####### for checking
	    $Lt=substr($line,22,4);
	    $Lt=~s/\s//mg;
	    $seq{$Lt}=substr($line,17,3); # template aa
	    $xc{$Lt}=substr($line,30,8);
	    $yc{$Lt}=substr($line,38,8);
	    $zc{$Lt}=substr($line,46,8);
	}
	####### record SG aotms------------->
	if($atom ne "C"){
	    if($atom ne "N"){
		if($atom ne "O"){
		    if(substr($atom,0,1) ne "H"){
			$Lt=substr($line,22,4);
			$Lt=~s/\s//mg;
			$n_atom{$Lt}++;
			$x{$Lt,$n_atom{$Lt}}=substr($line,30,8);
			$y{$Lt,$n_atom{$Lt}}=substr($line,38,8);
			$z{$Lt,$n_atom{$Lt}}=substr($line,46,8);
		    }
		}
	    }
	}
	####### SG atoms done ####################
    }
}
close(pdb);
$Lch=$Lt;

############### cc for orientation ###############
$xc{0}=$xc{1}-($xc{3}-$xc{2});
$yc{0}=$yc{1}-($yc{3}-$yc{2});
$zc{0}=$zc{1}-($zc{3}-$zc{2});
$xc{$Lch+1}=$xc{$Lch}+($xc{$Lch-1}-$xc{$Lch-2});
$yc{$Lch+1}=$yc{$Lch}+($yc{$Lch-1}-$yc{$Lch-2});
$zc{$Lch+1}=$zc{$Lch}+($zc{$Lch-1}-$zc{$Lch-2});
for($i=1;$i<=$Lch;$i++){
    $xm=$xc{$i}-$xc{$i-1};
    $ym=$yc{$i}-$yc{$i-1};
    $zm=$zc{$i}-$zc{$i-1};
    $aaa=sqrt($xm**2+$ym**2+$zm**2);
    $xm/=$aaa;
    $ym/=$aaa;
    $zm/=$aaa;
    $xp=$xc{$i+1}-$xc{$i};
    $yp=$yc{$i+1}-$yc{$i};
    $zp=$zc{$i+1}-$zc{$i};
    $aaa=sqrt($xp**2+$yp**2+$zp**2);
    $xp/=$aaa;
    $yp/=$aaa;
    $zp/=$aaa;
    $ccx{$i}=$xm-$xp;
    $ccy{$i}=$ym-$yp;
    $ccz{$i}=$zm-$zp;
    $aaa=sqrt($ccx{$i}**2+$ccy{$i}**2+$ccz{$i}**2);
    $ccx{$i}/=$aaa;
    $ccy{$i}/=$aaa;
    $ccz{$i}/=$aaa;
}

############## SG contacts ########################
undef %n_cont;
for($i=1;$i<=$Lch;$i++){
    for($j=$i+2;$j<=$Lch;$j++){
	for($m=1;$m<=$n_atom{$i};$m++){
	    for($n=1;$n<=$n_atom{$j};$n++){
		$dis=sqrt(($x{$i,$m}-$x{$j,$n})**2+
			  ($y{$i,$m}-$y{$j,$n})**2+
			  ($z{$i,$m}-$z{$j,$n})**2);
		if($dis<4.5){
		    $cc=$ccx{$i}*$ccx{$j}+$ccy{$i}*$ccy{$j}+$ccz{$i}*$ccz{$j};
		    if($cc>0.5){ 
			$orient=3; # parallel
		    }elsif($cc>-0.5){
			$orient=2;  # perpendical
		    }else{
			$orient=1; # antiparallel
		    }
		    $n_cont{$i}++;
		    $r_cont{$i,$n_cont{$i}}=$j;
		    $cc_cont{$i,$n_cont{$i}}=$orient;
		    $n_cont{$j}++;
		    $r_cont{$j,$n_cont{$j}}=$i;
		    $cc_cont{$j,$n_cont{$j}}=$orient;
		    goto pos4;
		}
	    }
	}
	pos4:;
    }
}

############# run STRIDE to get secondary structure #################
for($i=1;$i<=$Lch;$i++){
    $sec{$i}=1;
}
@lines=`/nfs/users/yzhang/bin/stride $s\.pdb 2> /dev/null`;
$k=0;
foreach $line(@lines){
    if($line=~/^ASG/){
	$k++;
	$sec_tmp=substr($line,24,1);
	$sec{$k}=2 if($sec_tmp eq "H");
	$sec{$k}=4 if($sec_tmp eq "E");
    }
    goto pos5 if($k >= $Lch);
}
 pos5:;

########### output SG contacts ####################
open(cnt,">$s\.cnt2");
printf cnt "$Lch\n";
for($i=1;$i<=$Lch;$i++){
    printf cnt "%5d %1s %1d %5d",$i,$ts{$seq{$i}},$sec{$i},$n_cont{$i};
    for($j=1;$j<=$n_cont{$i};$j++){
	printf cnt " %5d",$r_cont{$i,$j};
    }
    printf cnt "\n";
    printf cnt "%9s %5d"," ",$n_cont{$i};
    for($j=1;$j<=$n_cont{$i};$j++){
	printf cnt " %5d",$cc_cont{$i,$j};
    }
    printf cnt "\n";
}
close(cnt);

exit();
