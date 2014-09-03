#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp,junkName,Bye
from spicker.spickerYangResultsManager import spickOut,TMscore

inpHand('Usage: getBestTop.py [options]\noutput combo for top cluster, and best cluster if we pass the native file',
        ' -a _RA_location directory with spicker results, or tar file containing same',
        ' -b _A_nat native file (program will then output best.pdb) (def:none)',
        ' -c _A_outd output directory for top.pdb and best.pdb (def: ".")'
        ).parse(locals(),sys.argv)
currd=os.getcwd()
if not outd: outd=currd

#parse "location" variable
ind='' #direcotory where to find spicker results
if not os.path.exists(location): sys.exit(1)               #file/directory does not exis
if location[-4:]=='.tar':                                 #tar file
    ind=junkName() ; os.system('/bin/mkdir '+ind)         #temporary directory for spicker results
    os.system('tar xf '+location+' -C '+ind)              #untar contents to temporary dir
elif location[-8:]=='.tar.bz2' or location[-5:]=='.tbz2': #compressed tar file
    ind=junkName() ; os.system('/bin/mkdir '+ind)
    os.system('/bin/cp '+location+' '+ind)                #temporary directory for spicker results
    os.chdir(ind)
    loc=os.path.basename(location) #remove path from absolute file name
    os.system('/usr/bin/bunzip2 '+loc+' && tar xf '+loc[:-4]) #uncompress and untar to temp dir
    os.chdir(currd)
else :ind=location

#get top cluster
p=spickOut(dir=ind) #instantiate a spickOut object
if p.readError: sys.exit(1)
os.system('/bin/cp '+p.combo[p.densest]+' '+outd+'/top.pdb')
if location!=ind: os.system('/bin/rm -rf '+ind) #remove temporal "ind" directory
sys.exit(0)
