#!/usr/bin/python

"""For every reference sequence, find constitutive and variant
fragments. For each fragment, calculate its average disorder content"""

import os,sys,re
from utilities.small_utilities import chomp
from utilities.codedir import scratchdir
from average_disorder import read_dispro,read_vsl2,read_disembl,read_cst_var_lines

patterns={'cst':re.compile('1+'),'var':re.compile('0+')}
cmd='grep "00$" '+scratchdir+'/spliceVariants/out/success.list'
buf={'cst':'','var':''}
for ref in chomp(os.popen(cmd).readlines()):
    ref=ref[0:6]
    os.chdir(scratchdir+'/spliceVariants/out/'+ref+'/00')
    #read the three disorder predictions        
    disorders={'dispro':read_dispro(ref+'00.dispro'),
               'vsl2':read_vsl2(ref+'00.vsl2'),
               'vsl2b':read_vsl2(ref+'00.vsl2B'),
               'disembl':read_disembl(ref+'00.disembl')
               }
    #create the lines of 1's (for constitutive regions) and 0's (for variant regions)
    for line in read_cst_var_lines(ref):
        L=len(line)
        for type in ('cst','var'): #find first all constitutive fragments, then all variant ones
            for match in patterns[type].finditer(line):
                a,b=match.span() #position of the fragment, excluding "b"
                size=b-a #size of the fragment
                buf[type]+='%3d'%(size)
                for disorder in ('dispro','vsl2','vsl2b','disembl'):
                    dis=disorders[disorder] #prediction list containing disorder values
                    if dis:
                        sum=0.0
                        for i in range(a,b):sum+=dis[i]
                        buf[type]+=' %4.2lf'%(sum/(b-a))
                buf[type]+='\n'

print '# 1-size, 2-dispro, 3-vsl2, 4-vsl2b, 5-disembl'        
for type in ('cst','var'):
    print '#',type
    print '# 1  2    3    4    5'
    print buf[type]
sys.exit(0)
