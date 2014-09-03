#!/usr/bin/python

import sys,os,re
from stats.statistics import corrcoeff
from utilities.readingWritingFiles import read_column
from inputArgs.inputArgs import inpHand


def cf(pp,comment_marker='#',xcol=1,ycol=2):
    xl=read_column(pp,xcol,comment_mark=comment_marker,isFloat=1)
    yl=read_column(pp,ycol,comment_mark=comment_marker,isFloat=1)
    return corrcoeff(xl,yl)

if __name__=='__main__':
    #Parse the command line
    inpHand('Usage: corrcoef.py',
            ' -a _RA_pp data file',
            ' -b __xcol x-column (def:1)',
            ' -c __ycol y-column (def:2)',
            ' -d __mark commentary mark (def:#)',
            ).parse(locals(),sys.argv)
    #init default values
    if not xcol: xcol=1
    else: xcol=int(xcol)
    if not ycol: ycol=2
    else: ycol=int(ycol)
    if not mark: mark='#'

    print cf(pp,comment_marker=mark,xcol=xcol,ycol=ycol)

    sys.exit(0)
    
