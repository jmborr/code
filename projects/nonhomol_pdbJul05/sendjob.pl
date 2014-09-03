#!/usr/bin/perl
use Getopt::Std;

$log='/tmp/junklog';

&parse_input();
print "l=$l\n";
print "root=$root\n";
print "log=$log\n"; 

@h=`cat $l`;
foreach $id (@h){
    print $id; chomp($id);
    if( length($id) > 5 ){ $x=$id; $id=substr($x,2,5); }
    else{ $x=substr($id,1,1)."/$id"; }
    $exec="$root/bin/mkjob.pl -i $id -c $root/common/common.tbz2 -o $root/out/$x";
#    print "$exec\n";
    `bsub -q standard \"$exec\" >> $log`; sleep 5;
}

sub parse_input(){
    my $f=0;
    if(!getopts('l:r:g:h')){ &message(); }
    if($opt_l){ $l   =addAbsPath($opt_l); $f++; }
    if($opt_r){ $root=addAbsPath($opt_r); $f++; }
    if($opt_g){ $log =addAbsPath($opt_g);       }
    if($f<2){ &message(); }
    if($opt_h){ &message(); }
}

#add absolute path to files or directories (even if inexistent)
sub addAbsPath(){
    my $c=shift;
    $c=~s/^\.\///; #remove leading './' if present
    unless($c=~/^\//){ chomp($dir=`pwd`); $c="$dir/$c"; } #add abs path
    while($c=~s/\/\w+\/\.\.//){;} #substitute '/*/..' with '/' one by one
    $c=~s/\/$//; #remove last '/', if present
    return $c;
}

sub message(){
    system("clear");
    print "Usage: sendjob.pl -l -r [-g]\n";
    print " Mandatory:\n";
    print " -l list (file with list of headers)\n";
    print " -r root directory (root/bin or root/run)\n";
    print " Optional:\n";
    print " -g log file\n";
    print "\n";
    exit(0);
}
