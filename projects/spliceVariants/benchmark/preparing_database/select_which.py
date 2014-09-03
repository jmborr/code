#!/usr/bin/python

import sys
from utilities.small_utilities import chomp

pin=open('/gpfs1/scratch/jose/spliceVariants/benchmark/preparing_database/pdb.40.300.alns.dat','r')
pin.readline() #discard first line, since it is a comment
lines=chomp(pin.readlines())
pin.close()
homologs={}
partners={}
valid={}
for line in lines:
    list=line.split()
    h1=list[0]
    h2=list[1]
    l1=int(list[2])
    l2=int(list[3])
    lali=int(list[4]) #lenght of the alignment, including gaps
    nali=int(list[5]) #number of aligned residues
    nid=int(list[6])  #number of identical residues in the alignment
    cov=float(list[7])#coverage
    seqid=float(list[8]) #sequence identity over the aligned region
    nngap=int(list[9]) #lenght of N-gap
    ncgap=int(list[10]) #length of C-gap
    next=int(list[11]) #number of external gaps
    nint=int(list[12]) #number of internal gaps
    lbig=int(list[13]) #length of biggest internal gap segment

    #classify h1 and h2 as homologous to each other
    if not homologs.has_key(h1): homologs[h1]=[h2,]
    else: homologs[h1].append(h2)
    if not homologs.has_key(h2): homologs[h2]=[h1,]
    else: homologs[h2].append(h1)
    
    if next==0: #no external gaps
        if nint==1: #only one internal gap
            if lbig>5: #internal gap must be bigger than five residues
                if l1<200 and l2<200:
                    if seqid > 0.6:
                        if cov > 0.2:                            
                            valid[h1+' '+h2]=line
                            if not partners.has_key(h1): partners[h1]=[h2,]
                            else:  partners[h1].append(h2)
                            if not partners.has_key(h2): partners[h2]=[h1,]
                            else:  partners[h2].append(h1)


#for h1 in partners.keys():
#    for h2 in partners[h1]:
#        key=h1+' '+h2
#        if valid.has_key(key): print valid[key]
#sys.exit(1)
#print '**************************************************************'

#valid contains many redundancies of the type (h1,h2) and (h1,h3) but
#h2 and h3 are very similar. It doesn't pay to do TASSER simulations
#on h1, h2, and h3. We will impose transitivity in our set, that is,
#if the pair (h1,h2) satisfy the filtering conditions and the pair
#(h1,h3) satisfy the filtering conditions, then the pair (h2,h3)
#should also satisfy the conditions if we are to keep h1, h2, and
#h3. Hence, if (h1,h2) and (h1,h3) are in valid, then (h2,h3) should
#also be in valid

#require the transitive property
for h1 in partners.keys():
    #print h1, len(partners[h1])
    j=0
    while j<len(partners[h1])-1:
        h2=partners[h1][j]
        #print ' '+h2, j, len(partners[h1])
        k=j+1
        while k<len(partners[h1]):
            h3=partners[h1][k]
            #check for the transitive property
            if h3 in homologs[h2] and h3 not in partners[h2]:
                #print 'delete partners[',h1,'][',k,']'
                del partners[h1][k] #remove redundant entry from h1's list
                k-=1
                for l in range(0,len(partners[h3])): #same for h3's list
                    if partners[h3][l]==h1:
                        del partners[h3][l]
                        break
            k+=1
        j+=1

#output filtered list. note that valid has key h1+' '+h2 but has no
#h2+' '+h1 key. This way we avoid output the same info twice
for h1 in partners.keys():
    for h2 in partners[h1]:
        key=h1+' '+h2
        if valid.has_key(key): print valid[key]
