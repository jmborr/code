#!/usr/bin/perl

#===========================================================================
# generate a bunch of histograms. On the x-axis is the ratio of number of
# HN-HN contacts divided by the protein lenght. Each histogram is for a 
# different class of protein (alpha, beta, alpha/beta, ...)
#===========================================================================

#retrieve input, store it into memory
open(IN,"/home/jmborr/Projects/NMR_constraints/pdb1484.list");
@list=<IN>; close(IN); #list of proteins
open(IN,"/home/jmborr/Databases/scop/pdb1484.class.list");#alpha, beta, ..
while($line=<IN>){ $prot2class{substr($line,0,5)}=substr($line,7,4); }
close(IN);
open(IN,"/home/jmborr/Projects/NMR_constraints/analysis/lengths.dat");#length
while($line=<IN>){ $prot2clength{substr($line,0,5)}=substr($line,7,3); }
close(IN);

$I=30; #max number of bins;
$zoom=20; #factor to multyply the ratio #HNcontacts/prot-lentgth
$class={ #generic hash of pointers to lists
    'alph' =>[()], #[()] returns a pointer to a generic list
    'beta' =>[()],
    ' a/b' =>[()],
    ' a+b' =>[()],
    'mult' =>[()],
    'memb' =>[()],
    'smal' =>[()],
    'coil' =>[()],
    'lowr' =>[()],
    'pept' =>[()],
    'desg' =>[()]
	};

&init_hist($class,$I); #set all elements to zero. We pass class by reference

$data="/home/jmborr/Projects/NMR_constraints/data/";
foreach $d (@list){  #example : $d="0/101m_\n"
    chomp($d);
    $header=substr($d,2,5); print "$header\n";
    $c=$prot2class{$header}; #class of the protein
    $l=$prot2clength{$header}; #length of the protein
    $file=$data.$d."/".$header.".HN"; #file with HN-HN contacts
    open(IN,$file); chomp($nc=<IN>); close(IN); #number of HN-HN contacts
    $i=int($nc*$zoom/$l); #bin number (of histogram $c)
#    print "header=$header c=$c l=$l file=$file nc=$nc i=$i\n";
    $ll=$$class{$c}; #pointer to the appropriate histogram;
    $$ll[$i]++; #increment appropriate bin by one
}

$alys="/home/jmborr/Projects/NMR_constraints/analysis/";#dir to store output
&print_hists($class,$I);

#=============================================================================
sub init_hist{
    my ($histrf,$I)=@_; #print "class=$class I=$I\n"; exit;
    my ($key,$ll,$i);
    foreach $key (keys %$histrf){
#	print "init $key\n";
	$ll=$$histrf{$key}; #ref. to the list pointed to by %$histrf{$key}
        for($i=0;$i<$I;$i++){  $$ll[$i]=0;  } #initialize to zero
    }
#    print "finished init_hist\n";
}
#============================================================================
sub print_hists{
    my ($histrf,$I) = @_;
    my ($key,$ll,$out,$i);
    foreach $key (keys %$histrf){
	$ll=$$histrf{$key}; #ref. to the list pointed to by %$histrf{$key}
	$out=$alys."hnhnVSprotL".$key.".dat";   open(OUT,">$out");
	for($i=0;$i<$I;$i++){  printf(OUT "%3d%10d\n",$i,$$ll[$i]);  }
	close(OUT);
    }
}
#============================================================================
