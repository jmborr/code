package LOAD_MACHINES ;


$min_usage = 10 ; #a jobs is "important" if using more than 17% of the CPU
$machines="thalia water riskit jhilad urania luna melete vesta pan rushmore hypate clio erato malcolm miranda caliban ariel fuji logan rainier everest";

#below these lines is a better version of "sub machine_list". This version does
#not work when trying to "rsh" to some machines that are down.
#returns a scalar with a list of machines ordered by increasing number of jobs
#sub machine_list {
#  my ($machine, $n_jobs, $line, $cpu, $ordered_machines);
#  my (@MACHINES, @OUT, @list, @ordered_keys);
#  my %UP_MACHINES ;
#  @MACHINES=split(' ',$machines);
#now count number of jobs per machine    
#  foreach $machine (@MACHINES){
#    if( ! system("(rsh $machine uptime > /dev/null) >& /dev/null") ){
#      @OUT=`rsh $machine ps -uax`;
#      shift @OUT ; #first line of rsh is header, not data
#      $n_jobs=0;
#      foreach $line (@OUT){
#	@list=split(' ', $line);      $cpu=$list[2];
#	if( $cpu > $min_usage ){ $n_jobs++; }
#      }
#      $UP_MACHINES{$machine}=$n_jobs;
#    }
#  }
#sort machines per increasing number of jobs and return the list in a scalar
#  @ordered_keys = sort {$UP_MACHINES{$a}<=>$UP_MACHINES{$b}} keys %UP_MACHINES;
#  $ordered_machines = join(' ', @ordered_keys ) ;
#  return $ordered_machines ;
#}

#returns a scalar with a list of machines ordered by increasing number of jobs
sub machine_list {
  my ($i,$ordered_machines);
  my (@list,@list2);
  $ordered_machines="";
  @list=`ruptime -l -r`;
  for($i=0;$i<=$#list;$i++){
    $line=$list[$i];
    unless( $line =~ "down" ){#the particular machine should be "up"
      @list2=split(' ',$line); 
      $the_machine=lc($list2[0]);
      if($machines=~$the_machine){$ordered_machines.=" $the_machine";}
    }
  }
  return $ordered_machines ;
}

#returns a scalar containing a list of jobs that chew up CPU
sub machine_loads { #argument should be a SCALAR containing a list of machines
  my ($user, $machine, $load, $line, $cpu) ;
  my (@MACHINES, @OUT, @list) ;
  
  @MACHINES=  split(' ', $_[0]);
  foreach $machine (@MACHINES){
    @OUT=`rsh $machine ps -uax`;
    $load .= "\n***************************** $machine ******************************\n$OUT[0]";
    shift @OUT ;
    foreach $line (@OUT){
      @list=split(' ', $line);  $cpu=$list[2]; $owner=$list[0];
      if( $cpu > $min_usage ){ $load .= $line ; }
    }
  }
  return $load ;
}
1;# used files HAVE to return true!

#returns a reference to a hash containing the number of jobs (value) per
#machine (key). Argument should be a SCALAR containing a list machines
sub machines_n_jobs{
  my ($machine, $n_jobs, $cpu, $ref ) ;
  my (@MACHINES, @list ) ;
  my %UP_MACHINES ;
  @MACHINES=  split(' ', $_[0]);
  foreach $machine (@MACHINES){
    if( ! system("(rsh $machine uptime > /dev/null) >& /dev/null") ){
      @OUT=`rsh $machine ps -uax`;
      shift @OUT ; #first line of rsh is header, not data
      $n_jobs=0;
      foreach  (@OUT){
	@list=split ;  $cpu=$list[2];
	if( $cpu > $min_usage ){ $n_jobs++; }
      }
      $UP_MACHINES{$machine}=$n_jobs;
    }
  }
  $ref = \%UP_MACHINES ; return $ref ;
}

#returns a list with all jobs with big CPU usage for a particular user.
#Input is the username and a reference to an array where the subroutine will
#store the jobs. The last argument of each job line is the machine where the
#job resides.
sub list_jobs_per_user{ #(user,ref_to_list)
  my ($user, $list_ref,$command);
  my (@MACHINES);
  $user=shift; $list_ref=shift; #print "list_ref=$list_ref\n";
  @MACHINES=split(' ',machine_list()); #$,="-";print @MACHINES;exit(1);
  foreach $machine (@MACHINES){ #print "|$machine|\n";
    $command="rsh  $machine \"ps -fu $user\""; #print "$command\n";
    @OUT=`$command`; shift @OUT;#first line of @OUT is bogus
    foreach $line (@OUT){ #print $line;
      @list=split(' ', $line);      $cpu=$list[3];  #print "cpu=$cpu\n";
      if( $cpu > $min_usage ){ 
	chomp($line); $line .="   $machine";
	push(@{$list_ref},$line);
      }
    }
  }  #$,="\n";print "list_ref=".@{$list_ref}."\n";
}

#same as list_jobs_per_user, but there is no restriction on CPU usage
sub list_all_jobs_per_user{ #(user,ref_to_list)
  my ($user, $list_ref,$command);
  my (@MACHINES);
  $user=shift; $list_ref=shift; #print "list_ref=$list_ref\n";
  @MACHINES=split(' ',machine_list()); #$,="-";print @MACHINES;exit(1);
  foreach $machine (@MACHINES){ #print "|$machine|\n";
    $command="rsh  $machine \"ps -fu $user\""; #print "$command\n";
    @OUT=`$command`; shift @OUT;#first line of @OUT is bogus
    foreach $line (@OUT){ #print $line;
      @list=split(' ', $line);
      chomp($line); $line .="   $machine";
      push(@{$list_ref},$line);
    }
  }  #$,="\n";print "list_ref=".@{$list_ref}."\n";
}

#returns the number of jobs with big CPU usage for a particular user.
sub num_jobs_per_user{#(user)
  my ($user,$list_ref); my (@list_jobs); $list_ref=\@list_jobs;
  $user=shift; #print("user=$user\n");
  &list_jobs_per_user($user,$list_ref); #print("n=$#list_jobs\n");
  return $#list_jobs+1;

}

1; #stupid thing needed
