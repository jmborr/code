#!/usr/bin/perl -w
use lib  '/project1/PROTEIN/jmborr/PROTEIN1/Code/Perl/lib/' ;
use LOAD_MACHINES;
use strict ;

my ($N,$n,$user,$exec,$pre,$post,$n_argv,$suffix,$wd,$machines,$hash_ref,$job,
    $waitime);
my @input_files ;

$waitime=1200; #20 minutes waiting time;

sub test_input{
  my $n_argv = @_; #print "n_argv=$n_argv\n";
  if( $n_argv<5 || $n_argv>6 ){
    system("clear");
    print "USAGE: putJobs_remote.pl N exec pre post list [suffix] \n";
    print "  N   : do not put more jobs until number of jobs of this user in the cluster\n";
    print "        falls below N\n";
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

#Argument is a reference to a hash. Subroutine will output name of machine
#with the least number of jobs.
sub loose_machine {
  my $hash_ref = shift ;
  my @sorted = sort {${$hash_ref}{$a}<=>${$hash_ref}{$b}} keys %{$hash_ref};
  return $sorted[0] ;
}

&test_input( @ARGV ) ;
$N=$ARGV[0]; $exec=$ARGV[1]; $pre=$ARGV[2]; $post=$ARGV[3];
open(LIST, "$ARGV[4]") or die "ERROR: I can't open file \"$ARGV[4]\"\n";
@input_files=<LIST>; close(LIST);
$n_argv=@ARGV; if($n_argv==6){ $suffix = $ARGV[5] ; }
#print "exec=$exec\npre=$pre\npost=$post\nsuffix=$suffix\n";exit;
$wd=`pwd`; chomp($wd); #current working directory
$user=`whoami`;chomp($user);
#Force user to check jobs to be remotely submitted
$job="$exec $pre input_from_list $post ";
if($n_argv==6){ $job .= "input_from_list.$suffix"; }
print "Jobs will look like:\n   $job\n";
print "Are you OK with this?(y/n): ";$job=<STDIN>;chomp($job);
unless($job eq "y" || $job eq "Y"){print"Aborted!\n";exit(1);}

#get list of running machines
$machines = &LOAD_MACHINES::machine_list( ) ;
#hash_ref is a reference to a list of machines (keys) along with the number of
#jobs for each one (values)
$hash_ref=&LOAD_MACHINES::machines_n_jobs($machines);
$n=&LOAD_MACHINES::num_jobs_per_user($user);  #print "n=$n\n";
foreach my $input (@input_files){
  my ($machine, $arguments, $command ) ;
  chomp($input) ;
  while($n>=$N){
    sleep $waitime; #print"Now go to sleep, because n=$n and N=$N\n";
    $hash_ref=&LOAD_MACHINES::machines_n_jobs($machines); 
    $n=&LOAD_MACHINES::num_jobs_per_user($user);#print "n=$n\n";
  }
  $machine=&loose_machine($hash_ref);#find machine with smallest job load
  $arguments =  "$pre $input $post " ;
  if($n_argv==6){ $arguments .= "${input}.$suffix"; }
  $command = "rsh $machine \"cd $wd; $exec $arguments;\"&";
  print "$command\n\n";
  system($command);  $n++; #another jobs is sent  print"n=$n\n";
  sleep 9; #benefit of 9 seconds for rsh to stablish connection. Important!
  $$hash_ref{$machine}++ ; #now the loose machine has one more job!
}
