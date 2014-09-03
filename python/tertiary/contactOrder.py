#!/usr/bin/python

import os,sys
from inputArgs.inputArgs import inpHand

#"""
#For a given two-column file, find the contact order.
#"""

def contactOrder(pairs,L=1,cutoff=0):
    """calculate contact order of a list of pairs of integers"""
    if not pairs:
        sys.stderr.write('ERROR empty list\n')
        return None
    co=0.0  ;  npairs=0
    for (i,j) in pairs:
        d=abs(i-j)
        if d>cutoff:
            co+=d  ;  npairs+=1
    return co/(npairs*L)


if __name__=='__main__':   
    inpHand('Usage: contactOrder.py [options]',
            ' -a _AR_inmap contact map',
            ' -b __cutoff do not consider contacts i,i+k with k smaller than cutoff (def:2)',
            ' -c __L protein length',
            ).parse(locals(),sys.argv)

    if not cutoff: cutoff=2
    if not L: L=1
    
    pairs0=[ x.split() for x in open(inmap,'r').readlines() ]
    pairs=[]
    for pair in pairs0:
        i=int(pair[0]) ; j=int(pair[1])  ;  pairs.append( (i,j) )
    print contactOrder(pairs,cutoff=cutoff,L=int(L))
    sys.exit(0)
