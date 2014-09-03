#!/usr/bin/perl
######################################################################
# This program is to calculate Side-chain center of mass according to
# CA-coordinations. Comments should be addressed to zhang6@buffalo.edu
# Usage: CASG.pl CA_file_name, output will be "CA_file_name_SG"
######################################################################

use Math::Trig;
############# read side-chain parameters ###############
$ka1{GLY}=0.000;  $kb1{GLY}=0.000;  $kc1{GLY}=0.000; $ka2{GLY}=0.000;  $kb2{GLY}=0.000;  $kc2{GLY}=0.000; 
$ka1{ALA}=0.249;  $kb1{ALA}=-1.118; $kc1{ALA}=0.976; $ka2{ALA}=0.113;  $kb2{ALA}=-0.736; $kc2{ALA}=1.294; 
$ka1{SER}=0.169;  $kb1{SER}=-1.369; $kc1{SER}=1.103; $ka2{SER}=0.227;  $kb2{SER}=-0.966; $kc2{SER}=1.427; 
$ka1{CYS}=-0.073; $kb1{CYS}=-1.201; $kc1{CYS}=1.476; $ka2{CYS}=0.084;  $kb2{CYS}=-0.738; $kc2{CYS}=1.712; 
$ka1{VAL}=0.274;  $kb1{VAL}=-1.162; $kc1{VAL}=1.480; $ka2{VAL}=0.093;  $kb2{VAL}=-0.583; $kc2{VAL}=1.799; 
$ka1{THR}=0.090;  $kb1{THR}=-1.296; $kc1{THR}=1.346; $ka2{THR}=0.070;  $kb2{THR}=-0.854; $kc2{THR}=1.633; 
$ka1{ILE}=0.100;  $kb1{ILE}=-1.363; $kc1{ILE}=1.769; $ka2{ILE}=-0.105; $kb2{ILE}=-0.601; $kc2{ILE}=2.135; 
$ka1{PRO}=-0.743; $kb1{PRO}=-1.563; $kc1{PRO}=0.438; $ka2{PRO}=-0.980; $kb2{PRO}=-1.183; $kc2{PRO}=0.976; 
$ka1{MET}=-0.049; $kb1{MET}=-1.246; $kc1{MET}=2.308; $ka2{MET}=0.094;  $kb2{MET}=-0.723; $kc2{MET}=2.610; 
$ka1{ASP}=-0.221; $kb1{ASP}=-1.249; $kc1{ASP}=1.769; $ka2{ASP}=0.334;  $kb2{ASP}=-0.664; $kc2{ASP}=1.992; 
$ka1{ASN}=-0.357; $kb1{ASN}=-1.096; $kc1{ASN}=1.849; $ka2{ASN}=0.097;  $kb2{ASN}=-0.699; $kc2{ASN}=1.962; 
$ka1{LEU}=-0.057; $kb1{LEU}=-1.161; $kc1{LEU}=2.128; $ka2{LEU}=0.003;  $kb2{LEU}=-0.393; $kc2{LEU}=2.400; 
$ka1{LYS}=0.027;  $kb1{LYS}=-1.616; $kc1{LYS}=2.597; $ka2{LYS}=-0.019; $kb2{LYS}=-0.745; $kc2{LYS}=2.972; 
$ka1{GLU}=-0.013; $kb1{GLU}=-1.554; $kc1{GLU}=2.219; $ka2{GLU}=0.101;  $kb2{GLU}=-0.793; $kc2{GLU}=2.684; 
$ka1{GLN}=-0.086; $kb1{GLN}=-1.439; $kc1{GLN}=2.296; $ka2{GLN}=0.041;  $kb2{GLN}=-0.707; $kc2{GLN}=2.666; 
$ka1{ARG}=0.113;  $kb1{ARG}=-1.932; $kc1{ARG}=2.933; $ka2{ARG}=-0.020; $kb2{ARG}=-0.998; $kc2{ARG}=3.394; 
$ka1{HIS}=-0.221; $kb1{HIS}=-1.138; $kc1{HIS}=2.165; $ka2{HIS}=-0.133; $kb2{HIS}=-0.598; $kc2{HIS}=2.363; 
$ka1{PHE}=0.111;  $kb1{PHE}=-0.984; $kc1{PHE}=2.447; $ka2{PHE}=-0.363; $kb2{PHE}=-0.632; $kc2{PHE}=2.507; 
$ka1{TYR}=0.128;  $kb1{TYR}=-1.035; $kc1{TYR}=2.604; $ka2{TYR}=-0.375; $kb2{TYR}=-0.601; $kc2{TYR}=2.706; 
$ka1{TRP}=0.476;  $kb1{TRP}=-1.156; $kc1{TRP}=2.541; $ka2{TRP}=-0.058; $kb2{TRP}=-0.427; $kc2{TRP}=2.894; 
########################################################

############## read CA ################################
open(CA,$ARGV[0]); #open CA file
$Lch=0;
while($line=<CA>){
    #print $line;
    if(substr($line,0,6) eq "ATOM  " && substr($line,12,4)=~/CA/){
        #print $line;
	$Lch++;
	$seq{$Lch}=substr($line,17,3);
	$xx{$Lch}=substr($line,30,8);
	$yy{$Lch}=substr($line,38,8);
	$zz{$Lch}=substr($line,46,8);
    }
    goto pos1 if($Lch >1 && substr($line,0,3) eq "TER");
}
pos1:;
close(CA);
#print @seq;
$xx{0}=$xx{1}+($xx{2}-$xx{3});
$yy{0}=$yy{1}+($yy{2}-$yy{3});
$zz{0}=$zz{1}+($zz{2}-$zz{3});
$xx{$Lch+1}=$xx{$Lch}+($xx{$Lch-1}-$xx{$Lch-2});
$yy{$Lch+1}=$yy{$Lch}+($yy{$Lch-1}-$yy{$Lch-2});
$zz{$Lch+1}=$zz{$Lch}+($zz{$Lch-1}-$zz{$Lch-2});

######### calculate SG ########################
for($k=1;$k<=$Lch;$k++){
    ($xg{$k},$yg{$k},$zg{$k})=
	sidechain($xx{$k-1},$yy{$k-1},$zz{$k-1},
		  $xx{$k},$yy{$k},$zz{$k},
		  $xx{$k+1},$yy{$k+1},$zz{$k+1},
		  $seq{$k});
}

######### output SG ###########################
open(SG,">$ARGV[1]");
for($k=1;$k<=$Lch;$k++){
    printf SG
	"ATOM  %5s  CA  %3s  %4d    %8.3f%8.3f%8.3f\n",
	$k,$seq{$k},$k,$xx{$k},$yy{$k},$zz{$k};
    printf SG
	"ATOM  %5s  SG  %3s  %4d    %8.3f%8.3f%8.3f\n",
	$k,$seq{$k},$k,$xg{$k},$yg{$k},$zg{$k};
}
close(SG);

sub sidechain{
    my($xm,$ym,$zm,$x,$y,$z,$xp,$yp,$zp,$seq)=@_;
    ########### Angle ##########################
    my$a2=($xm-$x)**2+($ym-$y)**2+($zm-$z)**2;
    my$b2=($xp-$x)**2+($yp-$y)**2+($zp-$z)**2;
    my$c2=($xm-$xp)**2+($ym-$yp)**2+($zm-$zp)**2;
    $cosc=($a2+$b2-$c2)/(2*sqrt($a2*$b2));
    $angle=acos($cosc)/3.1415926*180;
    
    ############ vectors a,b,c#######################
    $vxm=$x-$xm;
    $vym=$y-$ym;
    $vzm=$z-$zm;
    $rm=sqrt($vxm*$vxm+$vym*$vym+$vzm*$vzm);
    $vxm/=$rm;    #vi
    $vym/=$rm;
    $vzm/=$rm;
    $vxp=$xp-$x;
    $vyp=$yp-$y;
    $vzp=$zp-$z;
    $rp=sqrt($vxp*$vxp+$vyp*$vyp+$vzp*$vzp);
    $vxp/=$rp;    #vj
    $vyp/=$rp;
    $vzp/=$rp;
    $ax=$vxm+$vxp;
    $ay=$vym+$vyp;
    $az=$vzm+$vzp;
    $aaa=sqrt($ax*$ax+$ay*$ay+$az*$az);
    $ax/=$aaa;    #a=(vi+vj)/|vi+vj|
    $ay/=$aaa;
    $az/=$aaa;
    $cx=$vxm-$vxp;
    $cy=$vym-$vyp;
    $cz=$vzm-$vzp;
    $ccc=sqrt($cx*$cx+$cy*$cy+$cz*$cz);
    $cx/=$ccc;    #c=(vi-vj)/|vi-vj|
    $cy/=$ccc;
    $cz/=$ccc;
    $bx=$cy*$az-$cz*$ay;
    $by=$cz*$ax-$cx*$az;
    $bz=$cx*$ay-$cy*$ax;
    $bbb=sqrt($bx*$bx+$by*$by+$bz*$bz);
    $bx/=$bbb;    #a(x)c
    $by/=$bbb;
    $bz/=$bbb;

    #######################################################################
    ##
    ## A=ka*a+kb*b+kc*c
    ## when a,b,c are unitary and perpenticular-->
    ## ka=A*a
    ## kb=A*b
    ## kc=A*c
    ##
    #######################################################################
    if($angle < 105){
	$xgp=$x+$ka1{$seq}*$ax+$kb1{$seq}*$bx+$kc1{$seq}*$cx;
	$ygp=$y+$ka1{$seq}*$ay+$kb1{$seq}*$by+$kc1{$seq}*$cy;
	$zgp=$z+$ka1{$seq}*$az+$kb1{$seq}*$bz+$kc1{$seq}*$cz;
    }else{
	$xgp=$x+$ka2{$seq}*$ax+$kb2{$seq}*$bx+$kc2{$seq}*$cx;
	$ygp=$y+$ka2{$seq}*$ay+$kb2{$seq}*$by+$kc2{$seq}*$cy;
	$zgp=$z+$ka2{$seq}*$az+$kb2{$seq}*$bz+$kc2{$seq}*$cz;
    }
    #printf "$xgp,$ygp,$zgp  $angle  *******\n";
    return($xgp,$ygp,$zgp);
}

exit();
