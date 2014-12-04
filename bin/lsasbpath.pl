#!/usr/bin/perl -w

chomp($currdir=`pwd`);

@list=`ls -1 -R`;
foreach $line (@list){
# print "line=$line";
 if($line=~/\:\n/){
  if($line eq ".:\n"){
   $subdir=$currdir.'/';
  }
  else{
   $subdir=$currdir.'/'.substr($line,2,$#line-1).'/';
  }
 }
 elsif($line eq "\n" || $line=~/..:/){;}
 else{
  chomp($line);
  print "$subdir"."$line\n";
 }
}
