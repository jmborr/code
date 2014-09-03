#!/usr/bin/python

import os,sys
from utilities.codedir import codedir,scratchdir
from utilities.small_utilities import fileDates

formatdb=codedir+'/bin/blast/formatdb'
pfilt=codedir+'/bin/blast/pfilt'
scratch_nrD=scratchdir+'/db/nr'
nrD='/library/yzhang/nr'

#***********************************************************************************
def format_archive():
    versionD=fileDates('nr.gz')['modify'] #YYYY-MM-DD-HH-MM-SS
    os.system('mkdir '+versionD+' && /bin/mv nr.gz '+versionD)
    os.chdir(versionD)
    os.system('gunzip nr.gz && sleep 1')
    os.system(pfilt+' nr > nr.filter && sleep 1')
    os.system(formatdb+' -i ./nr.filter -p T -o T')
    os.system('/bin/rm nr formatdb.log') #remove unfiltered database
    os.chdir('..')
    os.system('tar jcf '+versionD+'.tbz2 '+versionD+' && /bin/mv '+versionD+'/* .')
    os.system('bin/rmdir '+versionD)
    
#************************  MAIN  ***************************************************
os.chdir(scratch_nrD)
#os.system('wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz') #we do scp from local machine
format_archive()      #format new database



