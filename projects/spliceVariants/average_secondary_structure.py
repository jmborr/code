#!/usr/bin/python

"""For every reference sequence, find constitutive and variant
fragments. For each fragment, calculate alpha,beta,coil content
prediction or each and store in appropriate array. Then do the
average."""

import os,sys,re
from utilities.small_utilities import chomp
from utilities.codedir import scratchdir
from average_disorder import read_cst_var_lines,get_bin,get_size,maxBin,minChunk
from secondary.seqdat2ss import translateSec

avs={'cst':{'alpha':[0.0]*maxBin,'beta':[0.0]*maxBin,'coil':[0.0]*maxBin,
            'count':[0]*maxBin,'n':[0]*maxBin },
     'var':{'alpha':[0.0]*maxBin,'beta':[0.0]*maxBin,'coil':[0.0]*maxBin,
            'count':[0]*maxBin,'n':[0]*maxBin }
     }
s2sx={'H':'alpha','E':'beta','C':'coil'}
patterns={'cst':re.compile('1+'),'var':re.compile('0+')}
cmd='grep "00$" '+scratchdir+'/spliceVariants/out/success.list'
for ref in chomp(os.popen(cmd).readlines()):
    ref=ref[0:6]
    #obtain secondary structure and dump it in a string
    cmd = 'tar xf '+scratchdir+'/spliceVariants/out/'+ref+'/00/'+ref+'00.in.tasser.tar seq.dat'
    if not os.system(cmd):
        ss=translateSec('seq.dat')
        os.system('grep -e "[A-Z]    4 " seq.dat')
    else:
        continue
    #For each reference/isoform alignment, obtain a string on 1's and
    #0's ("1" for constitutive residues and "0" for variant residues
    for line in read_cst_var_lines(ref):
        L=len(line)
        for type in ('cst','var'): #find first all constitutive fragments, then all variant ones
            type_link=avs[type]
            count_link=avs[type]['count']
            n_link=avs[type]['n']
            #find all regions of type "type" in the alignment
            for match in patterns[type].finditer(line):
                a,b=match.span() #position of the fragment, excluding "b"
                size=b-a #size of the fragment
                if size<minChunk: continue #do not process too small fragments
                bin=get_bin(size) #; print size,bin
                count_link[bin]+=size
                n_link[bin]+=1
                for i in range(a,b):
                    s=ss[i] #find predicted secondary structure (in one-letter format)
                    type_link[ s2sx[s] ][bin]+=1

#print output
print '#1 bin-number, 2 related-fragment-size, 3 constitutive-number-of-fragments, 4 constitutive-alpha-content, 5 constitutive-beta-content, 6 constitutive-coil-content, 7 variant-number-of-fragments, 8 variant-alpha-content, 9 variant-beta-content, 10 variant-coil-content'
print '#1   2     3   4     5     6        7   8     9     10'
for bin in range(2,maxBin):
    buf='%2d %3d '%(bin,get_size(bin))
    for type in ('cst','var'):
        n_link=avs[type]['n']
        buf+='%5d  '%(n_link[bin])
        count_link=avs[type]['count']
        type_link=avs[type]
        if count_link[bin]>0:
            for sx in ('alpha', 'beta', 'coil'):
                buf+='%4.2f  '%(type_link[sx][bin]/count_link[bin])
        else:
            buf+='****  '*3
    print buf
    
sys.exit(0)
                    
        
