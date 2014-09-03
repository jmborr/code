'''
Created on Mar 7, 2013

@author: jmborr
'''

import numpy

def write_from_numpy(pt,arr,fmt=' %f',ncols=7,comment=''):
    """write numpy array into file"""
    buf=''
    if comment: buf+=comment+'\n'
    b=arr.ravel(); L=len(b) #reshape to one dimension
    index=0
    while True:
        #loop to write a whole line
        icol=0
        while icol<ncols:
            if index==L: break
            buf += fmt%(b[index])
            index += 1
            icol += 1
        buf+='\n'
        if index==L: break
    if pt:
        pt.write(buf)
        pt.flush()
    return buf

def histogram(a,bins=10,rg=None,normed=False):
    a = numpy.asarray(a).ravel()
    if not numpy.iterable(bins):
        if rg is None:
            rg = (a.min(), a.max())
        mn, mx = [mi+0.0 for mi in rg]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = numpy.linspace(mn, mx, bins, endpoint=False)

    # best block size probably depends on processor cache size
    block = 65536
    n = numpy.sort(a[:block]).searchsorted(bins)
    for i in range(block,len(a),block):
        n += numpy.sort(a[i:i+block]).searchsorted(bins)
    n = numpy.concatenate([n, [len(a)]])
    n = n[1:]-n[:-1]
    
    if normed:
        db = bins[1] - bins[0]
        return 1.0/(a.size*db) * n, bins
    else:
        return n, bins

def insertBfact(pdbin,blist,byres=False):
    """redefine the Bfactors according to list of numbers"""
    buf=''; bmin=min(blist);  f=100./(max(blist)-bmin)
    pin=open(pdbin); l=pin.readline().strip(); n=0
    prev = None
    while l:
        if len(l)<66: l+=' '*(66-len(l)) #pad blank spaces
        if l[:5]=='ATOM ':
            if byres:
                curr = l[22:26]
                if not prev: prev = curr
                buf+=l[:60]+'%6.2f'%(f*(blist[n]-bmin),)+l[66:]+'\n'
                if curr != prev:
                    prev = curr
                    n+=1
            else:
                buf+=l[:60]+'%6.2f'%(f*(blist[n]-bmin),)+l[66:]+'\n'
                n+=1
        else:
            buf+=l+'\n'
        l=pin.readline().strip()
    return buf
