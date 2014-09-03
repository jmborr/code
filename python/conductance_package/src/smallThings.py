#!/usr/bin/python
import sys,os,re,numpy,pdb
from small_utilities import chomp,Bye,junkName
from math import sqrt

def distFromStrings(v1,v2):
    [x1,y1,z1]=v1.split() ; x1=float(x1) ; y1=float(y1) ; z1=float(z1)
    [x2,y2,z2]=v2.split() ; x2=float(x2) ; y2=float(y2) ; z2=float(z2)
    return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) )

def insertBfact(pdbin,blist):
    """redefine the Bfactors according to list of numbers"""
    buf=''; bmin=min(blist);  f=100./(max(blist)-bmin)
    pin=open(pdbin); l=pin.readline().strip(); n=0
    while l:
        if l[:5]=='ATOM ':
            buf+=l[:60]+'%6.2f'%(f*(blist[n]-bmin),)
            if len(l)>=66: buf+=l[66:]
            n+=1
        else:
            buf+=l
        buf+='\n'; l=pin.readline().strip()
    return buf

def coords(pin):
    """retrieve coordinates"""
    xyz=[]; natoms=0
    l=pin.readline().strip(); n=0
    while l:
        if l[:5]=='ATOM ':
            xyz+=[float(x) for x in l[30:54].split()]; natoms+=1
        l=pin.readline()
    return numpy.array(xyz).reshape(natoms,3)

def neighborList(pin,co=4.0):
    """contact map and distances"""
    xyz=coords(pin); L=len(xyz)
    co2=co*co
    neighbors=[]; distances=[]
    for i in range(L):
        neighbors.append([]); distances.append([]);
    for i in range(L-1):
        r1=xyz[i]; neigh=[]; dist=[]; #print 'i=',i
        for j in range(i+1,L):
            rel=xyz[j]-r1; dd=(rel*rel).sum()
            if i!=j and dd<co2:
                d=sqrt(dd); neighbors[i].append(j); neighbors[j].append(i)
                distances[i].append(d); distances[j].append(d)
    #order lists by increasing distance
    for i in range(L):
        d=numpy.array(distances[i]); perm=numpy.argsort(d);
        d=d[::,perm]; n=numpy.array(neighbors[i])[::,perm]
        neighbors[i]=n; distances[i]=d
    return neighbors,distances,xyz
