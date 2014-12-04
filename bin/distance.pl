#!/usr/bin/perl -w

print "Enter two sets of coordenates\n" ;
chomp($_=<STDIN>) ; @r1=split(' ',$_) ;
chomp($_=<STDIN>) ; @r2=split(' ',$_) ;

printf("distance is %f\n",sqrt(($r1[0]-$r2[0])*($r1[0]-$r2[0]) + 
			       ($r1[1]-$r2[1])*($r1[1]-$r2[1]) + 
			       ($r1[2]-$r2[2])*($r1[2]-$r2[2])) ) ;
