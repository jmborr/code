#!/usr/bin/python

import os,sys,re
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand
from utilities.codedir import scratchdir

#global variables
pattern=re.compile('#(\d+-\d+)')

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: filterout_by_chunk_size.py [options]',
               '  -a __minsize remove missing atom chunks below this size (def: 4)',
               )
    ih.parse(locals(),sys.argv)

    if not minsize: minsize=4
    else: minsize=int(minsize)

    file=os.path.join(scratchdir,'missing_coordinates/preparing/missingXray.080503.txt')
    for line in os.popen('grep ">" '+file).readlines():
        nsmall=0
        nchunk=0
        for match in pattern.findall(line):
            (begin,end)=match.split('-')
            size=1+int(end)-int(begin)
            if size<minsize: nsmall+=1
            nchunk+=1
        if nsmall<nchunk: print line[1:5]+line[6]

    sys.exit(0)
        
            
                      
    
