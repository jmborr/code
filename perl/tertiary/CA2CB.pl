#!/usr/bin/perl
######################################################################
# This program is to calculate C-beta coordinates according to
# CA-coordinations. Comments should be addressed to zhang6@buffalo.edu
# Usage: CASG.pl CA_file_name, output will be "CA_file_name_CB"
######################################################################

use Math::Trig;
############# read side-chain parameters ###############
$ka1{GLY}=0.000; $kb1{GLY}=0.000;  $kc1{GLY}=0.000; $ka2{GLY}=0.000; $kb2{GLY}=0.000;  $kc2{GLY}=0.000;
$ka1{ALA}=0.249; $kb1{ALA}=-1.118; $kc1{ALA}=0.976; $ka2{ALA}=0.113; $kb2{ALA}=-0.736; $kc2{ALA}=1.294;
$ka1{SER}=0.192; $kb1{SER}=-1.089; $kc1{SER}=1.016; $ka2{SER}=0.130; $kb2{SER}=-0.765; $kc2{SER}=1.276;
$ka1{CYS}=0.189; $kb1{CYS}=-1.073; $kc1{CYS}=1.025; $ka2{CYS}=0.098; $kb2{CYS}=-0.726; $kc2{CYS}=1.308;
$ka1{VAL}=0.271; $kb1{VAL}=-1.095; $kc1{VAL}=1.019; $ka2{VAL}=0.090; $kb2{VAL}=-0.659; $kc2{VAL}=1.373;
$ka1{THR}=0.218; $kb1{THR}=-1.089; $kc1{THR}=1.029; $ka2{THR}=0.145; $kb2{THR}=-0.740; $kc2{THR}=1.310;
$ka1{ILE}=0.276; $kb1{ILE}=-1.100; $kc1{ILE}=1.015; $ka2{ILE}=0.085; $kb2{ILE}=-0.654; $kc2{ILE}=1.376;
$ka1{PRO}=0.154; $kb1{PRO}=-1.091; $kc1{PRO}=0.971; $ka2{PRO}=0.096; $kb2{PRO}=-0.827; $kc2{PRO}=1.258;
$ka1{MET}=0.235; $kb1{MET}=-1.098; $kc1{MET}=1.002; $ka2{MET}=0.102; $kb2{MET}=-0.702; $kc2{MET}=1.327;
$ka1{ASP}=0.140; $kb1{ASP}=-1.049; $kc1{ASP}=1.050; $ka2{ASP}=0.093; $kb2{ASP}=-0.749; $kc2{ASP}=1.287;
$ka1{ASN}=0.047; $kb1{ASN}=-1.034; $kc1{ASN}=1.064; $ka2{ASN}=0.073; $kb2{ASN}=-0.738; $kc2{ASN}=1.292;
$ka1{LEU}=0.244; $kb1{LEU}=-1.107; $kc1{LEU}=0.996; $ka2{LEU}=0.107; $kb2{LEU}=-0.692; $kc2{LEU}=1.334;
$ka1{LYS}=0.215; $kb1{LYS}=-1.105; $kc1{LYS}=0.994; $ka2{LYS}=0.108; $kb2{LYS}=-0.718; $kc2{LYS}=1.311;
$ka1{GLU}=0.239; $kb1{GLU}=-1.110; $kc1{GLU}=0.990; $ka2{GLU}=0.109; $kb2{GLU}=-0.717; $kc2{GLU}=1.309;
$ka1{GLN}=0.216; $kb1{GLN}=-1.102; $kc1{GLN}=0.999; $ka2{GLN}=0.101; $kb2{GLN}=-0.710; $kc2{GLN}=1.316;
$ka1{ARG}=0.222; $kb1{ARG}=-1.102; $kc1{ARG}=0.996; $ka2{ARG}=0.100; $kb2{ARG}=-0.719; $kc2{ARG}=1.310;
$ka1{HIS}=0.152; $kb1{HIS}=-1.070; $kc1{HIS}=1.031; $ka2{HIS}=0.082; $kb2{HIS}=-0.729; $kc2{HIS}=1.305;
$ka1{PHE}=0.210; $kb1{PHE}=-1.085; $kc1{PHE}=1.012; $ka2{PHE}=0.097; $kb2{PHE}=-0.706; $kc2{PHE}=1.323;
$ka1{TYR}=0.207; $kb1{TYR}=-1.087; $kc1{TYR}=1.010; $ka2{TYR}=0.098; $kb2{TYR}=-0.709; $kc2{TYR}=1.322;
$ka1{TRP}=0.228; $kb1{TRP}=-1.095; $kc1{TRP}=1.006; $ka2{TRP}=0.108; $kb2{TRP}=-0.707; $kc2{TRP}=1.322;
########################################################

############## read CA ################################
open(CA,$ARGV[0]); #open CA file
$Lch=0;
while($line=<CA>){
    if(substr($line,0,6) eq "ATOM  " && substr($line,12,4)=~/CA/){
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
	"ATOM  %5s  CB  %3s  %4d    %8.3f%8.3f%8.3f\n",
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
