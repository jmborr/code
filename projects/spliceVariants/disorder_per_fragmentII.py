#!/usr/bin/python

"""For every reference sequence, find variant fragments from
var_seq.dat. For each fragment, calculate its average disorder
content"""

import os,sys,re
from utilities.small_utilities import chomp
from utilities.codedir import scratchdir
from average_disorder import read_dispro,read_vsl2,read_disembl

varseqf=scratchdir+'/spliceVariants/preparing_database/var_seq.dat'
pin=open(varseqf,'r') 
cmd='grep "00$" '+scratchdir+'/spliceVariants/out/success.list'
bufs={'all':'','var':''}
for ref in chomp(os.popen(cmd).readlines()):
    #print ref
    os.chdir(scratchdir+'/spliceVariants/out/'+ref[0:6]+'/00')
    disorders={'dispro':read_dispro(ref+'.dispro'),
               'vsl2':read_vsl2(ref+'.vsl2'),
               'vsl2b':read_vsl2(ref+'.vsl2B'),
               'disembl':read_disembl(ref+'.disembl')
               }
    #computer average disorder of whole sequence
    L=0
    buf=''
    for disorder in ('dispro','vsl2','vsl2b','disembl'):
        dis=disorders[disorder]
        if dis:
            L=len(dis)
            sum=0.0
            for i in range(L):sum+=dis[i]
            buf+=' %4.2lf'%(sum/len(dis))
        else:
            buf+=' ****'
    bufs['all']+='%3d'%(L)+buf+'\n'
    
    #compute average disorder of its varian regions
    cmd='grep -b '+ref[0:6]+' '+varseqf+'|cut -d \':\' -f 1'
    byte_offset=int( os.popen(cmd).readline().strip() )
    pin.seek(byte_offset) ; pin.readline()
    while(True):
        l=pin.readline() #this should be a VAR_SET entry
        if l.find('VAR_SEQ') < 0: break
        items=l.split()
        a=int(items[2])-1
        b=int(items[3])
        #print a,b
        bufs['var']+='%3d'%(b-a)
        for disorder in ('dispro','vsl2','vsl2b','disembl'):
            #print disorder
            dis=disorders[disorder] #prediction list containing disorder values
            if dis:
                if b>len(dis):
                    bufs['var']+=' ****'
                    continue
                else:
                    sum=0.0
                    for i in range(a,b):sum+=dis[i]
                    bufs['var']+=' %4.2lf'%(sum/(b-a))
            else:
                bufs['var']+=' ****'
        bufs['var']+='\n'

#print output
print '# 1-size, 2-dispro, 3-vsl2, 4-vsl2b, 5-disembl'        
for type in ('all','var'):
    print '#',type
    print '# 1  2    3    4    5'
    print bufs[type]

sys.exit(0)
