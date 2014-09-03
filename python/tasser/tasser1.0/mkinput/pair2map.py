#!/usr/bin/python

"""
pair1.dat syntax:
            1MET
 -0.7  0.2 -0.4  3.1  2.8 -0.1 -0.1  2.8 -0.2 -0.2  3.0  2.7  3.3  4.7  3.1  1.0  2.1  1.5  2.4 -0.3 -1.0  1.9  1.2 -0.4  1.4
  2.4  2.5  3.1  4.3  0.6  1.0  0.7  3.5  0.6  0.3  1.4  0.3 -0.2  0.0  1.3  1.2  3.9  2.0  0.4  1.6  1.1  1.3  2.8  1.1  0.8
  2.7  0.7  1.7  2.1 -0.7  0.8  1.1  1.2 -0.1 -0.5  0.7  0.1 -0.7  0.2  0.3 -0.3  0.0  3.1  0.2 -0.2  1.0  0.0  1.7  1.8  2.9
  0.6  3.3  4.3  1.6  1.2  3.0  2.4  0.7  0.0  1.8  1.3  0.0  0.1  2.5  1.0 -0.2  0.5  2.2  1.4 -0.5  1.5  0.7 -0.9  0.1  3.0
 -0.3 -0.6  2.5  1.6  0.0  0.7  4.2  2.1  0.5  0.1  2.4 -0.1  1.5  2.3  2.3  0.0  1.8  1.8  2.8  1.0  0.9  0.4  0.2  2.0 -0.7
  0.1  0.7  0.4 -0.6  1.7  1.9 -0.2 -0.4  2.3  1.8  1.5  2.2  1.9  1.4  3.1  1.6  0.6  1.1  3.3  1.1  0.2  0.3  0.8  2.3  0.0
 -0.1  1.2  2.4 -0.5  1.4  1.1  1.1  3.0  3.3  0.2  1.2  0.2  0.9 -0.1  3.1  2.7
 =================================
 .
 .
 .
         166ASN
  1.4  2.9  2.8  1.7  0.1  0.1  3.2  2.4  0.0  2.6  1.6  2.2  2.7  4.0  3.0  0.4  1.3  1.5  2.4  2.4  1.3  2.6  2.3  2.9  1.4
  3.7  0.0  2.1  2.5  0.0  4.0  1.0  2.6  2.4  1.7  2.1  2.3  1.3  2.2  3.7  2.1  2.0  3.2  3.2  3.8  0.0  3.6  2.8  2.5  1.8
  2.5  4.1  1.5  0.0  0.1  1.6  3.6 -0.4  0.1 -0.1  0.9 -0.1  2.9  3.5  2.6  1.3  2.8  2.0  3.0  4.7  1.4  2.8  1.7  2.3  1.8
  2.4  1.8  0.2  3.4  2.5  2.7  0.0  2.9  1.4  2.9  4.1  1.3  0.0  0.0  3.9  0.1  3.2  1.5  4.5  2.3  0.0  0.0  1.7  0.9  2.1
  0.0  2.8  3.6  0.0  2.7  0.0  3.0  0.0  3.9  2.6  3.7  0.1  2.7  4.6  2.9  0.0  0.0  0.0  1.5  1.1  3.6  0.1  2.3  1.0  2.0
  0.0  0.0  1.0  4.1  0.0  2.9  2.1  0.1  2.2  3.0  3.5  2.5  0.4  1.7  2.9  0.0  2.9  4.6  0.0  1.5  2.4  2.1  4.4  1.2  3.6
  3.3  2.4  0.0  1.9  0.0  3.8  0.0  2.6  0.0  0.1  0.0  0.0  0.0  0.1  0.0  0.0
 =================================
"""

import sys,re,os
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp

def pair2map(pairf):
    #find sequence length first. Read all lines, reverse so that last lines becomes first, then join the lines
    lines=open(pairf,'r').readlines()  ;  lines.reverse()  ;  all=''.join(lines)
    L=int(re.compile('(\d+)[A-Z]+').search(all).group(1))
    nlines=int(L/25) #number of lines with 25 potential energy values
    nitems=[25]*nlines  ; nitems.append(L%25) #nitems[x] is number of pot.en.values in line x    
    nlines+=1 #; print nlines,nitems ; sys.exit(1)
    pin=open(pairf,'r')
    map=''
    for i in range(1,L+1):
        pin.readline() #header '            1MET'
        n=0
        for j in range(0,nlines):
            values=chomp(pin.readline()).split()            
            for k in range(0,nitems[j]):
                n+=1
                map=map+'%3d %3d %5s\n'%(i,n,values[k])
        pin.readline() #separator '================================='
    pin.close()
    return map

if __name__=='__main__':
    ih=inpHand('Usage: pair2map.py [options]',
               '  -a _AR_pairf pair file (pair1.dat or pair3.dat)',
               '  -b _AR_mapf map file')
    ih.parse(locals(),sys.argv)

    #reshape pair1.dat and output
    open(mapf,'w').write(pair2map(pairf))

    sys.exit(0)
