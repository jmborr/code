#!/usr/bin/python

#EXAMPLE:
#numFolds.py -a list -b d


import sys,os,re
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.small_utilities import chomp

ih=inpHand('Usage: numFolds.py [options]\nGiven a list of headers, and a SCOP prefix, will calculate the number of distint folds in the list\n',
           ' -a _RA_list list of headers',
           ' -b _R_prefix SCOP prefix (a:helix protein, b:beta, c:alpha/beta, d:alpha+beta)',
           ' -c _A_db SCOP database file (default: projectsdir/db/dir.cla.scop.txt_1.69)'
           )
ih.parse(locals(),sys.argv)
if not db: db=addAbsPath('projectsdir/db/dir.cla.scop.txt_1.69')
folds=[]
headers=[]
ldb=open(db,'r').readlines() #load the SCOP file
pat0=re.compile('\s'+prefix+'\.')
pat1=re.compile('cf=(\d+),')
for header in chomp(open(list,'r').readlines()):
    print header
    header=header.lower() #transformation needed to conform to SCOP
    pat=re.compile(header)
    for l in ldb:          #iter over SCOP entries
         if pat.search(l):              #found entry in SCOP
             if header in headers: continue #if header is multidomain, count only one 
             if pat0.search(l): #of adequate protein class
                 x=pat1.search(l)
                 if x:
                     fold=x.group(1)
                     #print header,fold; sys.exit(1)
                     if fold in folds: continue
                     headers.append(header)
                     folds.append(fold)
print len(folds)
sys.exit(0)
