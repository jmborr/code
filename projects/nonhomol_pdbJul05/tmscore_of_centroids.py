#!/usr/bin/python

import os,sys,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,split_file_by_signal
from jobs.job import pastry
from tmScore.tm_score import gettmCA
from random import randint

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: find_template_to_centroids.py [options]',
               '  -a _AR_listH list of header files',
               '  -b _AR_root directory (include xxxxx globbing)')
    ih.parse(locals(),sys.argv)
    
    p=re.compile('/x/xxxxx') ; q=re.compile('xxxxx') ; s=re.compile('TER')
    for header in chomp(open(listH,'r').readlines()):
        sbd=header[1] ;
        branch='/'+sbd+'/'+header
        dir=p.sub(branch,root) ; dir=q.sub(header,dir) #un-globb
        if not os.path.exists(dir+'/centroids.pdb'): continue
        superp=dir+'/tmsuptoCA.pdb' #where we store the superpositions
        if os.path.exists(superp): pastry('/bin/rm '+superp)
        tmpout=dir+'/'+`randint(0,9999)`  #temporary file where we store the current superposition
        #split centroids in a list of temporary files
        tmpfiles=split_file_by_signal(dir+'/centroids.pdb',s)
        for tmpfile in tmpfiles:
            gettmCA(tmpfile,dir+'/CA',sup=tmpout)
            pastry('cat '+tmpout+' >> '+superp)
            pastry('/bin/rm '+tmpfile) #remove the temporary file
        pastry('/bin/rm '+tmpout) #remove the temporary file
