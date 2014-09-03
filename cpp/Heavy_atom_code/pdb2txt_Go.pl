#!/usr/bin/perl -w
use strict ;

sub test_input {
  if( @_<4 || @_>10 ){
    print "Usage:perl pdb2txt_Go.pl nat.pdb box f_f init.txt f_Go e_Go e_HB sf inf fr\n";
    print "nat.pdb : native conformation of the protein in PDB format\n";
    print "box     : side of the simulation cube\n" ;
    print "f_f     : force field file with info on atoms (mass, radius, bonded\n" ;
    print "          and nonbonded interactions)\n";
    print "init.txt: output conformation where we include hydrogen bonds in the\n" ;
    print "f_go    : (1+go_co_f)*(VdW_i+VdW_j) will be Go range for i,j. Default f_go=0.35\n" ;
    print "e_go    : interaction strenght for Go model (must be negative number!. Default e_go=-1.0)\n";
    print "e_HB    : interaction strenght for hydrogen bond (must be negative number!. Default e_HB=-4.0)\n";
    print "sf      : percentage of how much we can penetrate the hard-core, or\n" ;
    print "          go beyond the link length (default sf=0.1)\n" ;
    print "          backbone and Go interactions in the sidechains.\n";
    print "inf     : height of the extra barrier included with sf. Default=1.00002\n";
    print "fr      : energy-rotamer factor (Er=Er(frc_fld)*f_rot). Default fr=0\n";
    print "\n";
    exit(1);
  }
}


&test_input( @ARGV ) ;
#obtain input
my ($ne,$nat,$box,$f_f,$f_Go,$e_Go,$e_HB,$sf,$init,$command,$inf,$fr,$query);
$ne=@ARGV;
$nat  = $ARGV[0];
$box  = $ARGV[1];
$f_f  = $ARGV[2];
$init = $ARGV[3];
$f_Go=0.35;$e_Go=-1.0;$e_HB=-4.0;$sf=0.1; $inf=1.000002; $fr=0;
if($ne>4){$f_Go = $ARGV[4];}
if($ne>5){$e_Go = $ARGV[5];}
if($ne>6){$e_HB = $ARGV[6];}
if($ne>7){$sf   = $ARGV[7];}
if($ne>8){$inf  = $ARGV[8];}
if($ne>9){$fr   = $ARGV[9];}
print"************ P A R A M E T E R S ***********";
print"    nat=$nat\n    box=$box\n    f_f=$f_f\n    init=$init\n    f_Go=$f_Go\n    e_Go=$e_Go\n    e_HB=$e_HB\n    sf=$sf\n    inf=$inf\n    fr=$fr\n";
print"Are these the corrects parameters?(y/n): ";$query=<STDIN>;chomp($query);
if($query ne "y" && $query ne "Y"){print"TRY AGAIN!\n"; exit(0);}
$command = "./pdb2txt_relax.x $nat $box $f_f junk_pre-relaxed" ;
print "$command\n";
system($command);
$command = "./txt2txt_Go2.x junk_pre-relaxed $nat $f_f $init $f_Go $e_Go $e_HB $sf $inf $fr"; 
print "$command\n";
system($command);
$command = ' \'rm\' junk_pre-relaxed' ;
print "$command\n";
system($command);
