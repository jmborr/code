#!/usr/bin/python

import sys,os,re,pdb

def resolveType(val):
    """read a string and return int, float, or string"""
    try:
        val=int(val); return val
    except:
        try:
            val=float(val); return val
        except:
            return val

def importPropFile(file):
    """read file containing lines like this:
    
    zwit..name ZGLY
    zwit..charge 0
    zwit..multiplicity 0
    id GLY

    will return a dictionary like this:
    { 'id':GLY, 'zwit':{ 'name':'ZGLY', 'charge':0, 'multiplicity':0}
    """
    props={}
    for l in open(file).readlines():
        key,val=l.split()
        val=resolveType(val)
        h=key.split('..');

        L=len(h)-1; d={h[L]:val}
        while L:
            L-=1
            d={h[L]:d}
            
        L=0; dd=props
        while L<len(h):
            if h[L] in dd.keys():
                d=d[ h[L] ]   #go down one level
                dd=dd[ h[L] ] #go down one level
                L+1
            else:
                dd.update(d)
                break
    return props

#coln: col number. The first column would have coln=1
def read_column(filename,coln,comment_mark='#',isInt=0,isFloat=0):
    col=[]
    n=coln-1
    pattern=re.compile(comment_mark)
    p=open(filename,'r')
    l=p.readline()    
    while l:
        if not pattern.match(l): #the line is not a comment
            v=l.split() #split by space(s)
            if isInt: v[n]=int(v[n])
            elif isFloat: v[n]=float(v[n])
            col.append(v[n])
        l=p.readline()
    p.close()
    return col

def read_to_cols(filename,comment='#',separator=' ',multiple=True,
                 outFmt='dictionary',outComm=False,xtype='str',
                 nskip=0):
    """read filename of fixed number of columns into a data structure
       outFmt: output Format (dict, list, numpy.array
       outComm: output comments found in filename)
       xtype: format to type (str, int, float)
       nskip: skip reading the first nskip lines (comment lines count too)"""
    from numpy import array as narray
    cols={}
    if outFmt in ('list', 'numpy.array'): cols = []
    buf = ''
    iskip=0
    for line in open(filename).readlines():
        if iskip<nskip or line[0]==comment:
            buf += line
            iskip +=1
            continue
        line=line.strip();
        items=[]
        if multiple:
            for x in line.split(separator):
                if x: items.append(x)
        else:
            items=line.split(separator);

        if xtype == 'int':
            X = [ int(x) for x in items ]
            items = X
        elif xtype == 'float':
            X = [ float(x) for x in items ]
            items = X
        L=len(items)
        if not cols:
            for i in range(L):
                if outFmt == 'dict': cols[i+1]=[]
                else: cols.append( [] )
        for i in range(L):
            if outFmt == 'dict': cols[i+1].append(items[i])
            else: cols[i].append(items[i])
    if outFmt == 'numpy.array':
        cols = narray( cols )
        cols = cols.transpose() #row index first, column index second
    if outComm:
        return cols,buf
    else:
        return cols
    
def write_lists(outf,lists,comment_line=''):
    "write to file a set of lists as columns"""
    p=open(outf,'w')
    p.write(comment_line+'\n')
    N=len(lists)
    M=len(lists[0])
    for i in range(M):
        l=''
        for j in range(N): l=l+' '+str(lists[j][i])
        p.write(l+'\n')
    p.close()

def arethere(files):
    """check if all files in the passed list exists"""
    for file in files:
        if not os.path.exists(file): return False
    return True

def read_to_numpy(pt,n,shape=False,comment='#',cformat=''):
    """read n numbers into a numpy array of shape=shape
    cformat: field for every number in printf format"""
    from numpy import zeros,array
    if cformat:
        L=int( re.compile('(\d+)').search(cformat).group(1) )
    arr2=[]
    line=pt.readline()
    while line and line[0]==comment: line=pt.readline()
    curr=pt.tell()
    m=0 #number of numbers read
    while line and m<n:
        curr=pt.tell()
        items=[]
        if cformat: #there may be two numbers glued together
            while len(line)>=L: items.append(line[0:L]); line=line[L:]
        else: items=line.split()
        for x in items:
            arr2.append( float(x) )
            m+=1
        line=pt.readline()
        while line and line[0]==comment: line=pt.readline()
    pt.seek(curr)
    arr = array( arr2 )
    if arr.any() and shape: arr=arr.reshape( shape )
    return arr

def write_from_numpy(pt,arr,format=' %f',ncols=7,comment=''):
    """write numpy array into file"""
    buf=''
    if comment: buf+=comment+'\n'
    b=arr.ravel(); L=len(b) #reshape to one dimension
    index=0
    while True:
        #loop to write a whole line
        for jndex in range(ncols):
            if index==L: break
            buf += format%(b[index])
            index += 1
        buf+='\n'
        if index==L: break
    if pt:
        pt.write(buf)
        pt.flush()
    return buf
