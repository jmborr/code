#!/usr/bin/python
import sys,os,re,numpy,pdb
from math import sqrt

def padLines( pdbfile, format='empty' ):
    """make 'ATOM' lines 80 charadters longs"""
    buf = ''
    for line in open( pdbfile ).readlines():
        if line[ 0:4 ] == 'ATOM':
            L = len( line.strip() )
            if L < 80:
                line = line.strip()
                if format == 'empty':
                    line += ' ' * ( 80 - L )
                elif format == 'fake':
                    """add fake occupancies and bfactors"""
                    if L < 60: line += '%6.2f' % (1.00,) #occupancy
                    if L < 66: line += '%6.2f' % (1.00,) #bfactor
                    line += ' '* ( 80 - len(line) )
                line += '\n'
        buf += line
    return buf
            
def distFromStrings(v1,v2):
    [x1,y1,z1]=v1.split() ; x1=float(x1) ; y1=float(y1) ; z1=float(z1)
    [x2,y2,z2]=v2.split() ; x2=float(x2) ; y2=float(y2) ; z2=float(z2)
    return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) )

def insertBfact(pdbin,blist,byres=False, rescale=False):
    """redefine the Bfactors according to list of numbers"""
    buf=''; bmin=min(blist);  f=100./(max(blist)-bmin)
    pin=open(pdbin); l=pin.readline().strip(); n=0
    prev = None
    while l:
        if len(l)<66: l+=' '*(66-len(l)) #pad blank spaces
        if l[:5]=='ATOM ':
            bf=blist[n]
            if rescale: bf=f*(blist[n]-bmin)
            if byres:
                curr = l[22:26]
                if not prev: prev = curr
                buf+=l[:60]+'%6.2f'%(bf,)+l[66:]+'\n'
                if curr != prev:
                    prev = curr
                    n+=1
            else:
                buf+=l[:60]+'%6.2f'%(bf,)+l[66:]+'\n'
                n+=1
        else:
            buf+=l+'\n'
        l=pin.readline().strip()
    return buf

def coordsII(pin, mask=None):
    """retrieve coordinates and more!"""
    xyz=[]; natoms=0; lines=[]
    l=pin.readline().strip(); n=0
    while l:
        if l[:5]=='ATOM ':
            if mask and l[12:16] in mask:
                    xyz+=[float(x) for x in l[30:54].split()]; natoms+=1
                    lines.append(l)
        l=pin.readline()
    return {'xyz':numpy.array(xyz).reshape(natoms,3),
            'pdblines':lines
            }
    
    
def coords(pin, mask=None):
    """retrieve coordinates"""
    xyz=[]; natoms=0
    l=pin.readline().strip(); n=0
    while l:
        if l[:5]=='ATOM ':
            if mask:
                if l[12:16] in mask:
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

def CBmap(pin,co=7.5,gap=2):
    """return contact map of CB atoms. Use CA if GLY
    gap: do not count |ires-jres|<gap"""
    props=coordsII(pin,mask=" CB , CA ")
    pxyz=props['xyz']; pl=props['pdblines']
    #filter lines by CB only
    co2=co*co; xyz=[]; lines=[]; natoms=0
    for i in range(len(pl)):
        if (" CB " in pl[i]) or (" CA " and 'GLY' in pl[i]):
            xyz.append(pxyz[i]); lines.append(pl[i]); natoms+=1
    #initialice neighbors and distances
    xyz=numpy.array(xyz).reshape(natoms,3); L=len(xyz)
    neighbors=[]; distances=[]
    for i in range(L):
        neighbors.append([]); distances.append([]);
    #find neighbors
    for i in range(L-1):
        r1=xyz[i]; neigh=[]; dist=[]; #print 'i=',i
        for j in range(i+1,L):
            rel=xyz[j]-r1; dd=(rel*rel).sum()
            if j-i>gap and dd<co2:
                d=sqrt(dd); neighbors[i].append(j); neighbors[j].append(i)
                distances[i].append(d); distances[j].append(d)
    #order list or each neighbor by increasing distance
    for i in range(L):
        d=numpy.array(distances[i]); perm=numpy.argsort(d);
        d=d[::,perm]; n=numpy.array(neighbors[i])[::,perm]
        neighbors[i]=n; distances[i]=d
    return {'neighbors':neighbors,'distances':distances,
            'xyz':xyz,'pdblines':lines}
   
    
    
    
