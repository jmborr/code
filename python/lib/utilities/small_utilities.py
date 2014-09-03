import re,os,sys
from random import randint

three2one={'ALA':'A','VAL':'V','PHE':'F','ILE':'I','LEU':'L','PRO':'P',
           'MET':'M','ASP':'D','GLU':'E','LYS':'K','ARG':'R','SER':'S',
           'THR':'T','TYR':'Y','HIS':'H','CYS':'C','ASN':'N','GLN':'Q',
           'TRP':'W','GLY':'G',}
one2three={'A':'ALA','V':'VAL','F':'PHE','I':'ILE','L':'LEU','P':'PRO',
           'M':'MET','D':'ASP','E':'GLU','K':'LYS','R':'ARG','S':'SER',
           'T':'THR','Y':'TYR','H':'HIS','C':'CYS','N':'ASN','Q':'GLN',
           'W':'TRP','G':'GLY',}

res={'access':'', 'modify':'', 'change':''}
   
def listToDict(list): #return a dictionary out of a list
    x=len(list)-1
    i=0
    dic={}
    while i<x:
        dic[list[i]]=list[i+1]
        i=i+2
    return dic

def chomp(object,replace=''):
    t=`type(object)`  #;print t # string representation of type
    if re.compile('str').search(t):
        if object:
            if(object[-1]=='\n'):  return object[:-1]+replace
        return object
    elif re.compile('list').search(t):
        for i in range(len(object)):
            if object[i]:
                if(object[i][-1]=='\n'):  object[i]=object[i][:-1]+replace
        #for line in object: no good, because line isn't a reference to an element of object
            #if(line[-1]=='\n'):
                #line=line[:-1]
        return object
    elif re.compile('tuple').search(t):
        stderr.write("ERROR: chomp can't change in place a tuple !")

def mysplit(s,append='',delimiter='',chunksize=None,count=None):
    l =[]
    if not delimiter and not chunksize:  #split into as many characters in s
        for char in s: l.append(char+append)    #run standard split function
    elif delimiter: 
        if count: segments=s.split(delimiter,count)        
        else    : segments=s.split(delimiter)
        if append:
            for segment in segments:  l.append=segment+append    
    else:                               #split into chunks of chunksize
        i=0 ; n=len(s) 
        while i<n:
            l.append(s[i:i+chunksize]+append)
            i+=chunksize
    return l

#'signal' is a regular expression to match, 'include=1' includes the
#line matching the signal in the previous file, 'skip' means do not
#check for signal in the first skip lines

def split_file_by_signal(filename,signal,include=1,skip=0,tmproot=`randint(0,99999)`):
    tmpfiles=[]
    if not os.path.exists(filename):
        sys.stderr.write('ERROR: file '+filename+' does not exists')
        return []
    id=0
    tmproot=os.path.dirname(filename)+'/'+tmproot #put temp files in same dir as filename
    nlines=0
    tmpfile=tmproot+'.%04d'%(id) ; outp=open(tmpfile,'w')
    inp=open(filename,'r')
    line=inp.readline()
    #do not check the first 'skip' lines
    for times in range(0,skip):
        outp.write(line) ; nlines=nlines+1
        line=inp.readline()
    #check
    while line:
        if signal.match(line):
            if include: outp.write(line) ; nlines=nlines+1 #include in previous file
            outp.close() ; tmpfiles.append(tmpfile) ; nlines=0
            id=id+1
            tmpfile=tmproot+'.%04d'%(id) ; outp=open(tmpfile,'w')
            if not include: outp.write(line) ; nlines=nlines+1 #include in next file
        else:
            outp.write(line) ; nlines=nlines+1
        line=inp.readline()        
    inp.close()
    outp.close()
    if nlines>0: tmpfiles.append(tmpfile)
    else: os.system('/bin/rm '+tmpfile)
    return tmpfiles

def my_dict_sort(dict,kv='key',mode='ascending'):
    kvs=('key','value')
    modes=('ascending','descending')
    if kv not in kvs: kv='keys'
    if mode not in modes: mode='ascending'
    sorted=[]
    if kv=='value':
        values=dict.values()
        if mode=='ascending':
            values.sort()
            for value in values:
                for key in dict.keys():
                    if dict[key]==value:
                        sorted.append((key,value))
    return sorted

def replace_at(string,position,replacement):
    l=len(string)-1 #last index
    if position==0:
        if l>0: return replacement+string[1:]
        else  : return replacement
    if position==l: return string[0:l]+replacement
    return string[0:l-1]+replacement+string[l:]

def junkName():
    junk='junk'+`randint(0,999)`
    while os.path.exists(junk): junk='junk'+`randint(0,999)`
    return junk

def todayDate():
    return chomp(os.popen('date +%F').readline())

pat=':\s(\d+-\d+-\d+\s\d+:\d+:\d)'
res={'access':re.compile('Access'+pat), 'modify':re.compile('Modify'+pat), 'change':re.compile('Change'+pat)}
def fileDates(file):
    times={}
    if not os.path.exists(file): return times
    all=''.join( os.popen('stat '+file).readlines() ) #dump everything into a string
    for key in res.keys():
        date=res[key].search(all).group(1)
        times[key]=date.replace(':','-').replace(' ','-')
    return times

#we pass string to "abort", but any object to "Bye" function (see below "abort")
def abort(message):
    sys.stderr.write(message)
    sys.stderr.write('\n')
    sys.exit(1)

def Bye(object):
    print object
    sys.exit(1)

def unTARme(f,wd='',fileL=[],taropts=''):

    """untar some file (it may be also compressed) to some directory

    unTARme(wd='')
    wd: optional directory where to untar. Otherwise create temp. dir.
    function returns the directory where untarred results are
    fileL: list of files to untar, in case we do not want to untar all files"""

    #create directory
    if not wd: wd=os.getcwd()+'/'+junkName()
    os.system('/bin/mkdir -p '+wd)
    #check tar method used
    opts='xf'
    if '.bz2' in f or '.tbz2' in f or '.bz' in f: opts='jxf'
    if '.gz'  in f or '.tgz' in f: opts='zxf'
    if taropts: opts=taropts
    #untar
    files=''
    if fileL: files=' '.join(fileL) #we pass a list of files to untar
    cmd='/bin/tar '+opts+' '+f+' -C '+wd+' '+files
    if os.system(cmd): #error when trying to untar to directory
        os.system('/bin/rm -r '+wd)
        wd=''
    return wd
