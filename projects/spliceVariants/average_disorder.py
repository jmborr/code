#!/usr/bin/python

"""For every reference sequence, find constitutive and variant
fragments. For each fragment, calculate disorder prediction or each
residue and store in appropriate array. Then do the average. The
result is the average disorder prediction of a residue in a
constitutive or variant fragment of a particular size"""

import os,sys,re
from utilities.small_utilities import chomp
from utilities.codedir import scratchdir
from seq.alignYangManageResults import alignOut
from math import log
#   =============================================

def check_probs(xx):
    for x in xx:
        if x>1: print x

#   =============================================
        
def read_disembl(file):
    """ read disorder prediction from ref+'00'.disembl file"""
    #find sequence length
    try:
        pin=open(file,'r')
    except:
        return []
    while(True):
        l=pin.readline()
        if l[0]=='>':
            N=len(pin.readline().strip())
            pin.close()
            break
        
    disembl=[0.25]*N #init all to default ordered value
    pattern=re.compile('\d+-\d+')
    #keep only REM465 and HOTLOOPS, which refer to intrinsic disorder
    for l in os.popen('grep -e ">" '+file).readlines()[1:]:
        for match in pattern.finditer(l):
            a,b=[ int(x)-1 for x in match.group().split('-') ]
            disembl[a:b+1]=[0.75]*(b+1-a)
    check_probs(disembl)
    return disembl 

#   =============================================

def read_dispro(file):
    """ read disorder prediction from ref+'00'.dispro file"""
    try:
        pin=open(file,'r')
    except:
        return []
    dispro=[ float(x) for x in pin.readlines()[-1].split() ]
    check_probs(dispro)    
    return dispro

#   =============================================

def read_vsl2(file):
    """ read disorder prediction from ref+'00'.vsl2 file"""
    try:
        pin=open(file,'r')
    except:
        return []
    l=pin.readline()
    while l.find('----------------------') < 0:  l=pin.readline()
    l=pin.readline()
    vsl2=[]
    while l.find('==============') < 0:
        vsl2.append( float(l.split()[2]) )
        l=pin.readline()
    check_probs(vsl2)
    return vsl2
        
#   =============================================    
    
def read_cst_var_lines(ref):
    pidx=open(scratchdir+'/spliceVariants/preparing_database/var0_var2.alg.idx','r')
    pin=open(scratchdir+'/spliceVariants/preparing_database/var0_var2.alg','r')
    lines=[]
    """find all (reference,isoform) alignments and for each aligment,
    return a line of 1's and 0's. 1's for constitutive residues, 0's
    for variant residues """    
    pattern=re.compile('0+')
    for line in pidx.readlines():
        pos,ref2,var=line.split()
        if ref2.find(ref)==0: #found one alignment
            #print ref2,var
            pin.seek(int(pos))
            alg=alignOut(pin).alg #read one alignment
            ll=alg.isAlignedLines()[1]
            #smooth the line of 1's and 0's
            ll=ll.replace('010','000')
            ll=ll.replace('101','111')
            lines.append( ll.replace('001100','000000') )
    return lines
        
#   =============================================    

def get_bin(L):
    return int( log(L)/log(2) )

def get_size(bin):
    return int(2**(bin-1)*3) # (2**bin+2**(bin+1))/2

#   =============================================    


maxBin=9
minChunk=4

if __name__=='__main__':
    curd=os.getcwd()
    avs={'cst':{'dispro':[0.0]*maxBin,'vsl2':[0.0]*maxBin,'disembl':[0.0]*maxBin,'count':[0]*maxBin,'n':[0]*maxBin },
         'var':{'dispro':[0.0]*maxBin,'vsl2':[0.0]*maxBin,'disembl':[0.0]*maxBin,'count':[0]*maxBin,'n':[0]*maxBin }
     }
    patterns={'cst':re.compile('1+'),'var':re.compile('0+')}
    cmd='grep "00$" '+scratchdir+'/spliceVariants/out/success.list'
    for ref in chomp(os.popen(cmd).readlines()):
        ref=ref[0:6]
        os.chdir(scratchdir+'/spliceVariants/out/'+ref+'/00')
        #read the three disorder predictions
        disorders={'dispro':read_dispro(ref+'00.dispro'),
                   'vsl2':read_vsl2(ref+'00.vsl2'),
                   'disembl':read_disembl(ref+'00.disembl')
                   }
        #print read_cst_var_lines(ref)#;sys.exit(1)
        for line in read_cst_var_lines(ref): #go through all alignemnts in which "ref" participates
            for type in ('cst','var'): #find first all const. fragments, then all variant ones
                count_link=avs[type]['count'] #a handy link
                n_link=avs[type]['n'] #a handy link
                for match in patterns[type].finditer(line): #go through all fragments
                    #print type,match.group()
                    a,b=match.span() #position of the fragment, excluding "b"
                    size=b-a #size of the fragment
                    if size<minChunk: continue #do not process too small fragments
                    bin=get_bin(size) #; print size,bin
                    count_link[bin]+=size  #there are "size" residues in this fragments
                    n_link[bin]+=1 #another fragment
                    #go through every type of disorder prediction
                    for disorder in ('dispro','vsl2','disembl'):
                        list=avs[type][disorder]#a handy link
                        dis=disorders[disorder] #prediction list containing disorder values
                        if dis:
                            sum=0.0
                            for i in range(a,b):sum+=dis[i]
                            list[bin]+=sum #dump all values together


    #print output
    print '#1 bin-number, 2 L, 3 number-of-fragments, 4 constitutive-dispro, 5 const-vsl2, 6 const-disembl, 7 number-of-fragments 8 variant-dispro, 9 var-vsl2, 10 var-disembl'
    print '#1   2     3     4     5      6     7    8     9     10'
    os.chdir(curd)
    for bin in range(maxBin):
        buf='%2d %3d  '%(bin,get_size(bin))
        for type in ('cst','var'):
            count_link=avs[type]['count']
            n_link=avs[type]['n']
            buf+='%5d  '%(n_link[bin])
            if count_link[bin]>0:
                for disorder in ('dispro','vsl2','disembl'):
                    list=avs[type][disorder]
                    buf+='%4.2f  '%(list[bin]/count_link[bin])
            else:
                    buf+='      '*3
        print buf
        
    sys.exit(0)
