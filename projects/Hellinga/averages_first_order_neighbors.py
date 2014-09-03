#!/usr/bin/python

import os,sys
from utilities.small_utilities import chomp,Bye

#calculate lengths of proteins
lengths={}
pdb='/gpfs1/scratch/adrian/dat/pdb_jul07'
for file in os.popen('ls -1 '+pdb+'/CA|grep pdb').readlines():
    file=file.strip()
    n=os.popen('wc -l '+pdb+'/CA/'+file).readline().split()[0] #;Bye(n)
    lengths[file[0:5]]=int(n)

#init factors
fkeys=(0,1,2,3,4,5,6,7,8,9,1000)
factors={}
nneighbors={}
for fkey in fkeys: factors[fkey]=1+float(fkey)/10

print '#     type   factors: 1.0  1.1  1.2  1.3  1.4  1.5  1.6  1.7  1.8  1.9 101'

dd='/gpfs1/scratch/adrian/protsup'
for type in ('1_alpha', '2_beta', '3_alpha_dash_beta', '4_alpha_plus_beta'): #cycle over sec.str.
    nstructures=0
    for fkey in fkeys: nneighbors[fkey]=0
    for query in open(dd+'/scop/'+type).readlines(): #cycle over queries
        query=query.strip()
        nstructures+=1
        if query not in lengths.keys(): continue
        queryL=lengths[query]
        for line in open(dd+'/netdist0.4/'+query+'.list'): #cycle over neighbors
            target,order,tm=line.split()
            if int(order)>1: break #we finished reading first order neighbors
            if target not in lengths.keys(): continue
            for fkey in fkeys:
                if queryL*factors[fkey]>lengths[target]: nneighbors[fkey]+=1
    #calculate averages
    buf='%20s'%(type,)
    for fkey in fkeys:
        buf+=' %4d'%(nneighbors[fkey]/nstructures)
    print buf

sys.exit(0)

