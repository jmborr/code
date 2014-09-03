#!/usr/bin/perl
$data="/home/jmborr/Projects/NMR_constraints/data/";
$alys="/home/jmborr/Projects/NMR_constraints/analysis/";
#list of all histograms of sequence separation in an anonymous hash.
$hist={
       'HN'      =>[()],
       };

&init_hist($hist); #set all elements to zero. We pass hists by reference

while($d=<>){
    chomp($d);
    $header=substr($d,2,5); 
    $d=$data.$d."/";  #print "header=$header d=$d\n";exit(0);
    foreach $key (keys %$hist){ #for each of the three types of NOEs.
	$file=$d.$header.".".$key;
	#print "$execut $pdb $type $outf{$type}\n";
	&add_to_histo($file,$hist,$key); #update all the histograms
    }
}

&print_hists($hist);
#=============================================================================
sub init_hist{
    $histrf=shift;
    foreach $key (keys %$histrf){
	print "init $key\n";
	$l=$$histrf{$key}; #ref. to the list pointed to by %$histrf{$key}
        for($i=0;$i<400;$i++){  $$l[$i]=0;  }
    }
    print "finished init_hist\n";
}
#============================================================================
sub add_to_histo{
    $f=shift;     open(IN,$f);      $n=<IN>; #just read the first line
    ($histrf,$key)=@_; #print "f= $f histrf=$histrf key=$key\n" ; exit(0);
    $ll=$$histrf{$key};
    @l=<IN>;
    foreach $cont (@l){
	$i1=substr($cont,0,10);
	$i2=substr($cont,10,10);
	$sep=abs($i2-$i1); #print"sep=$sep\n";
	$$ll[$sep]++;
    }
    close(IN); 
}
#============================================================================
sub print_hists{
    $histrf=shift;
    foreach $key (keys %$histrf){
	$l=$$histrf{$key}; #ref. to the list pointed to by %$histrf{$key}
	$out=$alys."seqSepHist".$key;   open(OUT,">$out");
	for($i=0;$i<400;$i++){  printf(OUT "%3d%10d\n",$i,$$l[$i]);  }
	close(OUT);
    }
}
#============================================================================
