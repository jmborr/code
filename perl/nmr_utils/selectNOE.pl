#!/usr/bin/perl

#select a subset of NOE's restraints so that the histograms of the
#subset and all the restrains versus contact order are similar (except
#for a zoom factor)

use Getopt::Std;

$e=10;
&parse_input();

#read NOEs into a hash with randomized keys
open(in,$inf) or die "Can't find file $inf\n";
chomp($N=<in>); $N=0;
while($x=<in>){
    chomp($x);
    split ' ',$x;
    $i=sprintf("%d",abs($_[0]-$_[1])); #sequence separation
    unless($i<$e){ 
	$p{int(rand 1000000)}=$x; $N++; #assign a random key to this noe
	#print "$x\n";
    }
}
#print "$N\n"; exit(0);
#if($N>0){ print "$header\n"; }
close(in); #foreach $key (keys %p){ print "$key $p{$key}\n";}
if($N==0){
    open(out,">$outf"); print out "0\n"; close(out); exit(0);
}
#project the hash into a list, sorting by numeric value of the keys. This way we randomly order the noe's
@pp = sort {$a<=>$b} keys %p; #order the keys 
for($i=0;$i<@pp;$i++){ $pp[$i]=$p{$pp[$i]}; } #substitute the key by the value
#$,="\n";print @pp;exit(0);

#n is number of bins in the histogram of contacts versus contact
#order, and $M is the number of contacts to be picked. $w is bin
#width.
$M=1+sprintf("%d",$l/$m); 
$r=$M/$N; 
open(out,">$outf");
if($r>=1){#too few contacts, just dump all contacts into the output file
    print out "$N\n";
    print out @pp ;
    exit(0);
}
#print "$l $M $r\n";exit(0);
print out "$M\n";
$n=1+sprintf("%d",$M/2);
#print "n=$n\n";exit(0);
$w=($l-$e)/$n; #there are no contacts with contact order smaller than $e
#store the pais in a list of lists. Each sublist contains all the
#pairs associated to a particular bin. The index of the parent list is
#the bin index. @hn contains number of elements in the array divided by $m
#print "w=$w\n";exit(0);
foreach $x (@pp){ #x is a pair of contacting residues
    ($a,$b)=split ' ',$x;
    #sequence separation in units of the bin widht, thus $i is the bin number
    $i=sprintf("%d",abs($a-$b)/$w); #printf "%d\n",$b-$a; #find bin index for the pair
    #print "i=$i\n";
    if(defined($h[$i])){ #push returns number of elements in @{$h[$i]} 
	$hn[$i] = push @{$h[$i]}, $x ; 
    }
    else{#create a reference to a list, $x is first element 
	$hn[$i] = 1 ;
	$h[$i] = [$x] ;
    } 
}

for($i=0;$i<=$n;$i++){ 
    if(defined($hn[$i])){ $hn[$i]*=$r; }
    else{ $hn[$i]=0; }
#    print "$hn[$i]\n";
}

#$,="\n"; for($i=0;$i<=$n;$i++){ printf "%lf %lf\n",($i+1)*$w,$hn[$i]; }
#pick now $n contacts by: 
# (1) find the maximum in @hn, say it is $hn[$i]
# (2) decrease the value, $hn[$i]--;
# (3) $pair = shift @{$h[$i]}; this is the contact we pick. It is random
#     because we randomized the ordering of all the contacts
# (4) repeat (1)
while($M>0 && $N>0){ #it might be that $M > $N
    $max=-1;
    #look beginning from high bin indexes. If there's a tie, we select
    #the bin with highest index, thus we'll draw a contact with high
    #contact order.
    for($i=@h-1;$i>=0;$i--){ if($hn[$i]>$max){ $max=$hn[$i]; $im=$i; } }
    $hn[$im]--; 
    $i=shift @{$h[$im]}; #printf "%2d %7.2lf $i\n",$im,$w*$im;
    print out "$i\n";
    $M--; $N--;
}
close(out);

sub parse_input(){
    my $f=0;
    if(!getopts('i:o:m:l:e:a:h')){ &message(); }
    if($opt_i){ $inf=$opt_i; $f++; }
    if($opt_o){ $outf=$opt_o; $f++; }
    if($opt_m){ $m=$opt_m; $f++; }
    if($opt_l){ $l=$opt_l; $f++; }
    if($opt_e){ $e=$opt_e; $f++; }
    if($opt_a){ $header=$opt_a; }
    if($f<5){ &message(); } 
    if($opt_h){ &message(); }
}

sub message(){
    system("clear");
    print "Usage: selectNOE.pl -i noefile -o outfile -m m -l l -e e\n";
    print " Mandatory:\n";
    print " -i noefile: list of NOE's, and first line is N, number of NOE's\n";
    print " -o outfile: output subset of selected restraints\n";
    print " -m m: select 1+int(N/m) restraints \n";
    print " -l l: lenght of the protein\n";
    print " -e e: exclude contacts with contact order smaller than e (def=10)";
    print " Optional:\n";
    print " -a header five letter code\n";
    print "\n";
    exit(0);
}
