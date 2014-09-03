#!/usr/bin/python
import os,sys
import re,string
from small_utilities import listToDict,Bye
from codedir import prefixd

#pick arguments, associate pairs of flags and values. Flags with no values are
#associated a null ('') value
def flagsToDict(list,dic=None):
    if not dic: dic={}
    if not list: return dic #we passed an empty list
    fre=re.compile('^-([a-z,A-Z])$') #detect the '-' sign and store matched letter
    while fre.search(list[0]): #list[0] has the '-' sign
        l=list[0]
        list=list[1:]
        l=fre.sub(r"\1",l) #remove the '-' sign by keeping only matched letter
        dic[l]=''
        if not list: return dic #return if we "ate" the whole list
        if not fre.search(list[0]): #the flag contained an argument
            dic[l]=list[0]
            list=list[1:]
            if not list: break
    return dic

def dictToFlags(dic):
    flags=''
    for (flag,val) in dic.items(): flags+=' -'+flag+' '+val+' '#lots spaces,just in case
    return flags
    
#update options by overwriting default options
def update_opts(opt,opt0):
    optL=opt.split()  ;  opt0L=opt0.split()
    dic=flagsToDict(opt0L) #initialize dictionary of passed options
    dic=flagsToDict(optL,dic=dic) #overwrite the dictionary
    return dictToFlags(dic)

def does_exists(file):
    if not os.path.exists(file):
        sys.stderr.write('ERROR: '+file+' does not exists\n')
        sys.exit(1)

def addAbsPath(local): #add absolute path name
    if local=='.': local=''
    p=re.compile('^\.\/') #remove leading './', if present
    local=p.sub('',local,count=1)
    p=re.compile('\/$') #remove last '/', if present
    local=p.sub('',local,count=1)
    #substitute prefix keywrd for appropriate path
    for prefix in prefixd.keys():
        p=re.compile(prefix)
        if p.match(local):#only match at the beginning of the 'local' string
            local=p.sub(prefixd[prefix],local,count=1)
            break
    p=re.compile('^\/')#begins with '/', ie, already global
    if not p.match(local): #no match, thus no global
        if local: #subdirectory of current directory
            local=os.getcwd()+'/'+local
        else: #local is just the current directory
            local=os.getcwd()
    #translate all the  trailing ts (xxx/..)
    p=re.compile('.*\/\w+\/\.\.')
    while p.match(local):
        q=re.compile('\/\w+\/\.\.')
        local=q.sub('',local,count=1)
    return local

    
#deglobb substitute yyyyy for xxxxx, y/yyyyy for x/xxxxx
#deglobb substiture xxxxx for the header, and x/xxxxx for header[1]/header
pre=re.compile('/x/xxxxx')
pre2=re.compile('xxxxx')
pre3=re.compile('/y/yyyyy')
pre4=re.compile('yyyyy')
def deglobb(globbed,header):
    subd='/'+header[1]+'/'+header
    #if passed xxxxx, then switch to header
    degl=pre.sub(subd,globbed)
    degl=pre2.sub(header,degl)
    #if passed yyyyy, then switch to xxxxx
    degl=pre3.sub('/x/xxxxx',degl)
    degl=pre4.sub('xxxxx',degl)
    return degl

def isglobbed(globbed):
    if  pre2.search(globbed) or pre4.search(globbed): return 1
    return 0

def extractHeader(globbed,instance):    
    """
    given a globbed name and one instance (the globbed is replaced
    with a header), extract the header
    """
    subdirs_g=globbed.split('/')
    nitems=len(subdirs_g)
    subdirs_i=instance.split('/')
    if nitems!=len(subdirs_i):
        Bye('ERROR: '+globbed+' and '+instance+' do not have same number of subdirectories')
    position=0
    for i in range(nitems):
        if subdirs_g[i]=='xxxxx':
            position=i
            break
    return subdirs_i[position]

class inpHand: 
    def __init__(self,*descriptor):
        self.dic={}     #flags-name pairs
        self.optdic={}  #flags-special pairs (A,R)
        self.f=0        #number of required flags
        fre=re.compile('\s*-([a-z])\s+')        #find flag
        wre=re.compile('\s+_([A-Z]*)_(\w+)\s+') #find leading _??_name
        self.message=descriptor[0]+'\n'
        reqbuf=" Required arguments\n"
        optbuf=" Optional arguments\n"
        fo=0 #number of optional flags
        for line in descriptor[1:]:
            flag=fre.search(line,1).group(1)
            pattern2=wre.search(line,1)
            self.dic[flag]=pattern2.group(2)
            self.optdic[flag]=pattern2.group(1)
            if re.compile('R').search(self.optdic[flag]): #required options
                self.f=self.f+1
                line3=wre.sub(' ',line,1)[0:60]+'\n'
                line2=wre.sub(' ',line,1)[60:]
                while(len(line2)>0):
                    line3=line3+'    '+line2[0:55]+'\n'
                    line2=line2[55:]
                reqbuf+=line3
                #self.message=self.message+line3
                    #self.message=self.message+wre.sub(' ',line,1)+"\n"
            else:
                fo=fo+1
                line3=wre.sub(' ',line,1)[0:60]+'\n'
                line2=wre.sub(' ',line,1)[60:]
                while(len(line2)>0):
                    line3=line3+'    '+line2[0:55]+'\n'
                    line2=line2[55:]
                optbuf+=line3
        if self.f>0: self.message+=reqbuf
        if fo>0: self.message+=optbuf

        
    def abort(self,error='',clear=1):
        if clear==1: os.system("clear")
        if error: print error
        print self.message
        sys.exit(1)

    def parse(self,initvar,args):
        #introduce __doc__ of initvar
        lines=self.message.split('\n')
        if initvar['__doc__']:
            self.message=lines[0]+'\n'+initvar['__doc__']+'\n'+'\n'.join( lines[1:] )
        if re.compile('--help|-h').search(''.join(args[1:])): self.abort()
        #first put args[1:] to a dictionary, removing the '-' from the flags        
        idic=flagsToDict(args[1:])
        #idic=listToDict(map((lambda l: re.compile('^-([a-z])$').sub(r"\1",l)), args[1:]))
        #did we input all required flags?
        n=0 #number of required flags passed
        for flag in idic.keys():
            if not flag in self.optdic.keys():
                self.abort()                
            if re.compile('R').search(self.optdic[flag]):
                n=n+1
            if re.compile('A').search(self.optdic[flag]):
                idic[flag]=addAbsPath(idic[flag])
        if n != self.f: #the number of required flags passed is insufficient
            self.abort()
        #abort if -h flag present in the passed flags
        if 'h' in idic.keys() :
            self.abort()
        #instantiate variable names in caller namespace
        for flag in self.dic.keys():
            if flag in idic.keys():
                initvar[self.dic[flag]]=idic[flag]
            else:
                initvar[self.dic[flag]]=''


def resolveArgument(arg,outfmt='',out=None):

    """resolve the type of the argument

    resolveArgument(arg,outfmt='',out=None) returns one of these keywords:
    _FILEPOINTER_ : arg is a pointer to file
    _FILENAME_    : arg is the name of a file
    _FILESTRING_  : arg is a string
    _LIST_        : arg is a list

    In addition, we can specify some action to be done in arg if we
    want to change the argument type. For example, if arg is the name
    of a file and we want the file to be stored in a list, then we can
    set outfmt=_LIST_ and out=mylist. The function will open file arg
    and store its contents in list mylist."""

    infmt=None
    if isinstance(arg,'list'): infmt='_LIST_'
    elif isinstance(arg,'file'): infmt='_FILEPOINTER_'
    elif isinstance(arg,'sting'):
        if '\n' in arg: infmt='_FILESTRING_'
        else: infmt='_FILENAME_'

    if outfmt: #we want some output
        if infmt==outfmt:
            out=arg
        elif infmt=='_FILENAME_':
            if outfmt=='_FILEPOINTER_' : out=open(arg,'r')
            elif outfmt=='_FILESTRING_': out=''.join(open(arg,'r').readlines())
            elif outfmt=='_LIST_' :      out=open(arg,'r').readlines()
        elif infmt=='_FILEPOINTER_':
            if outfmt=='_FILESTRING_': out=''.join(arg.readlines())
            elif outfmt=='_LIST_' :    out=arg.readlines()
        elif infmt=='_LIST_':
            if outfmt=='_FILESTRING_': out=''.join(arg)
        elif infmt=='_FILESTRING_':
            if outfmt=='_LIST_': out=arg.split('\n')
        else:
            sys.stderr.write('ERROR(resolveArgument) cannot transform' +infmt+' to '+outfmt+'\n')
    return infmt


#self-tester
if __name__ == '__main__':
    inpHand('Usage: myfirst.py -i -o [-p -x]',
            '  -i _AR_input      listfile'    ,
            '  -o _AR_outdir   output dir'   ,
            '  -p __optdir  optdir ja ja ja',
            '  -x __optflag  optional flag'   ).parse(locals(),sys.argv)

    print (input,outdir,optdir,optflag)
    sys.exit(0)
