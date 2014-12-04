#!/bin/tcsh

setenv PGI /usr/pgi
setenv fcflagp blank
setenv g98root /home/software/Gaussian
source $g98root/g98/bsd/g98.login

/home/software/Gaussian/g98_run "$*" 

