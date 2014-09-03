#!/usr/bin/python

import os,sys
from utilities.small_utilities import chomp,Bye

rootd='/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/tmAll2All'
#init certain variables
N=16430
factor=100
nbins=factor+1
long2short={}
headers=[]
histograms={} #number of targets with a particular TM-score to query
sectypes=('a','b','c','d')
sechist={}    #percentage of targets of a particular secondary structure type with a...
nsec={}
for sectype in sectypes:
    sechist[sectype]=[0,]*nbins #TM-score to query below a particular cutoff
    nsec[sectype]=0
    
for line in open(rootd+'/structures/structure2link.txt','r').readlines():
    long,short=line.split()
    long2short[long]=short
    headers.append(short)
    histograms[short]=[0,]*nbins

#cycle over all pairwise TM-scores
nline=0
inpf=open(rootd+'/tmAll2All.dat')
line=inpf.readline()
while line:
    nline+=1
    querylong,targetlong,tm=line.split()
    tm=float(tm)
    n=int(factor*tm)    #;print tm,n
    query=long2short[querylong]
    target=long2short[targetlong]
    histograms[query][n]+=1
    histograms[target][n]+=1    
    line=inpf.readline()

    
#output histogram and cummulative histogram for each file
os.system('/bin/mkdir -p '+rootd+'/neighbors_histograms')
for header in headers:
    sectype=header[0] #first letter of the header is the secondary structure type
    secref=sechist[sectype]
    nsec[sectype]+=1
    ref=histograms[header]
    cumhist=[0.0,]*nbins #number of targets with a TM-score to query below a particular cutoff
    cumhist[nbins-1]=ref[nbins-1]
    j=nbins-2
    for i in range(nbins-1):
        cumhist[j]=cumhist[j+1]+ref[j]
        j=j-1
    outf=open(rootd+'/neighbors_histograms/'+header+'.his','w')
    for i in range(nbins):
        secref[i]=secref[i]+cumhist[i]
        cumhist[i]=float(cumhist[i])*100/N #give as a percentage
        outf.write( '%4.2lf  %5d  %6.3lf\n'%( float(i)/factor,ref[i],cumhist[i]) )    
    outf.close()

#output cummulative histogram for each secondary structure type
for sectype in sectypes:
    print sectype,nsec[sectype]
    secref=sechist[sectype]
    buf=''
    for i in range(nbins):
        buf+=' %4.2lf %6.3lf\n'%(float(i)/factor,float(secref[i])*100/(N*nsec[sectype]),)
    open(rootd+'/neighbors_histograms/'+sectype+'.his','w').write(buf+'\n')
