package small_util ;
require Exporter;

our @ISA       =qw(Exporter);
our @EXPORT    =qw(addAbsPath pastry); #exported symbols by default
our @EXPORT_OK = qw($codedir $projectsdir $db $pdbnonhom $scratchdir);
our $VERSION   = 1.00;

$pdbnonhom="/library/pdb_jan06";

#find which code mirror is accessible
@mirrors=('/gpfs1/active/jose');

$mirror='/gpfs1/active/jose';
$codedir="$mirror/code";
$projectsdir="$mirror/projects";
$db="$projectsdir/db";
$scratchdir='/gpfs1/scratch/jose';

#add absolute path to a file (allows for /.. )
sub addAbsPath {
    my $c=shift;
    if($c eq '.'){ $c='';}  #print "c=$c\n";
    $c=~s/^\.\///; #remove leading './' if present 
    unless($c=~/^\//){ 
	chomp($dir=`pwd`);
	if($c eq ''){$c=$dir;}
	else{$c="$dir/$c"; } #add abs path
    }
    while($c=~s/\/\w+\/\.\.//){;} #substitute '/*/..' with '/' one by one
    return $c;
}

sub pastry {
    my $command=shift;
    if(system($command)){ #returns non-zero if system call fails
	print STDERR "ERROR while executing: $command\n" ;
	`sync`;
	exit(1);
    }
}

1; #stupid thing needed
