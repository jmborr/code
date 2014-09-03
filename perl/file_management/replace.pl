#!/usr/bin/perl
if(@ARGV<=2){
  system("clear");
  print("Usage: perl replaceLine.pl fileName 'oldLine' 'newLine' \n\n") ;
  print("opens file 'fileName', and replaces 'oldLine' with 'newLine' as many times as 'oldLine' appears in 'fileName'.\n\n");
  exit(1);
}

$fileName = $ARGV[ 0 ] ; open(IN, $fileName) ; @file = <IN> ; close(IN) ;

$oldLine = $ARGV[ 1 ] ;

$newLine = $ARGV[ 2 ] ;

open(OUT, ">junkReplaceLine.txt" ) ;

foreach $line (@file){

  chomp( $line ) ;

  if( $line eq $oldLine ){ print(OUT "$newLine\n") ; }
  else{ print(OUT "$line\n") ; }
  
}

system("mv -f junkReplaceLine.txt $fileName") ;
