#!/usr/bin/perl
use strict ;

sub test_input{
  my $n_argv = @_; #print "n_argv=$n_argv\n";
  if( $n_argv<4 || $n_argv>5 ){
    system("clear");
    print "USAGE: putJobs_local.pl exec pre post list [suffix] \n";
    print "  exec: name of executable file\n";
    print "  pre : arguments before input file\n";
    print "  post: arguments after input file\n";
    print "  list: file with list of input files\n";
    print "  suffixL (optiona)\n";
    print "Jobs are sent either as \"exec pre input_file post\"\n" ;
    print "or as: \"exec pre input_file post input_file.suffix\"\n";
    print "\n";
    exit;
  }
}

&test_input( @ARGV ) ;
my ($exec,$pre,$post,$suffix,$n_argv,$command, $input,$job);
my (@input_files ) ;

$exec   = $ARGV[0] ;
$pre = $ARGV[1] ;
$post = $ARGV[2] ;
open(LIST, "$ARGV[3]"); @input_files=<LIST>; close(LIST);
$n_argv=@ARGV; if($n_argv==5){ $suffix = $ARGV[4] ; }
#print "exec=$exec\npre=$pre\npost=$post\nsuffix=$suffix\n";exit;

$job="$exec $pre input_from_list $post ";
if($n_argv==5){ $job .= "input_from_list.$suffix"; }
print "Jobs will look like:\n   $job\n";
print "Are you OK with this?(y/n): ";$job=<STDIN>;chomp($job);
unless($job eq "y" || $job eq "Y"){print"Aborted!\n";exit(1);}

foreach $input (@input_files){
  chomp($input);
  $command = "$exec $pre $input $post " ;
  if($n_argv==5){ $command .= "${input}.$suffix"; }
  print "$command\n";
  system($command);
  wait;
}
