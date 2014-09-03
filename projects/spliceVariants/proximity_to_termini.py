#!/usr/bin/python

"""For every reference sequence, find constitutive and variant
fragments. For each fragment, calculate its proximity to the nearest
terminal."""

import os,sys,re
from utilities.small_utilities import chomp
from utilities.codedir import scratchdir
from average_disorder import read_cst_var_lines,get_bin,get_size,maxBin,minChunk

avs={'cst':{'dist':[0.0]*maxBin,'dist2':[0.0]*maxBin,'size':[0.0]*maxBin,'n':[0]*maxBin },
     'var':{'dist':[0.0]*maxBin,'dist2':[0.0]*maxBin,'size':[0.0]*maxBin,'n':[0]*maxBin }
     }
patterns={'cst':re.compile('1+'),'var':re.compile('0+')}
cmd='grep "00$" '+scratchdir+'/spliceVariants/out/success.list'
for ref in chomp(os.popen(cmd).readlines()):
    ref=ref[0:6]
    for line in read_cst_var_lines(ref):
        L=len(line)
        for type in ('cst','var'): #find first all constitutive fragments, then all variant ones
            n_link=avs[type]['n'] #a handy link
            dis_link=avs[type]['dist'] #a handy link
            dis2_link=avs[type]['dist2'] #a handy link
            size_link=avs[type]['size'] #a handy link
            for match in patterns[type].finditer(line):
                a,b=match.span() #position of the fragment, excluding "b"
                size=b-a #size of the fragment
                if size<minChunk: continue #do not process too small fragments
                bin=get_bin(size) #; print size,bin
                n_link[bin]+=1 #another fragment
                size_link[bin]+=size
                #distance from center of fragment to either termini
                d=float(min((a+b)/2,L-(a+b)/2))/L
                dis_link[bin]+=d
                #distance from edges of fragment to either termini
                d=float(min(a,L-b))/L
                dis2_link[bin]+=d

#print output
print '#1 bin-number, 2 related-fragment-size, 3 constitutive-number-of-fragments, 4 constitutive-average-fragment-size, 5 constitutive-porcentual-distance-of-center-fragment-to-termini, 6 constitutive-porcentual-distance-of-edge-fragment-to-termini, 7 variant-number-of-fragments, 8 variant-average-fragment-size, 9 variant-porcentual-distance-of-edge-fragment-to-termini, 10 variant-porcentual-distance-of-edge-fragment-to-termini'
print '#1   2    3     4    5     6     7     8    9     10'
for bin in range(2,maxBin):
    buf='%2d %3d '%(bin,get_size(bin))
    for type in ('cst','var'):
        n_link=avs[type]['n']
        size_link=avs[type]['size'] #a handy link
        buf+='%5d  '%(n_link[bin])
        if n_link[bin]>0:
            dis_link=avs[type]['dist'] #a handy link
            dis2_link=avs[type]['dist2'] #a handy link
            buf+='%3d  %4.2f  %4.2f '%(size_link[bin]/n_link[bin],
                                        dis_link[bin]/n_link[bin],dis2_link[bin]/n_link[bin])
    print buf
        
sys.exit(0)            
