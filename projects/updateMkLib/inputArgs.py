#!/usr/bin/python
import os,sys
import re,string
from small_utilities import listToDict
from codedir import codedir,libseq,prefixd

#pick arguments, associate pairs of flags and values. Flags with no values are
#associated a null ('') value
def flagsToDict(line):
    dic={}
    if not line: return dic #we passed an empty line
    fre=re.compile('^-([a-z])$') #detect the '-' sign and store matched letter
    while fre.search(line[0]): #line[0] has the '-' sign
        l=line[0]
        line=line[1:]
        l=fre.sub(r"\1",l) #remove the '-' sign by keeping only matched letter
        dic[l]=''
        if not line: return dic #return if we "ate" the whole line
        if not fre.search(line[0]): #the flag contained an argument
            dic[l]=line[0]
            line=line[1:]
            if not line: return dic
    return dic

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
    #translate all the  trailing dots (xxx/..)
    p=re.compile('.*\/\w+\/\.\.')
    while p.match(local):
        q=re.compile('\/\w+\/\.\.')
        local=q.sub('',local,count=1)
    return local
        
class inpHand: 
    def __init__(self,*descriptor):
        self.dic={}     #flags-name pairs
        self.optdic={}  #flags-special pairs (A,R)
        self.f=0        #number of required flags
        fre=re.compile('\s*-([a-z])\s+')        #find flag
        wre=re.compile('\s+_([A-Z]*)_(\w+)\s+') #find leading _??_name
        self.message=descriptor[0]+"\n Required arguments\n"
        opt=" Optional arguments\n"
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
                self.message=self.message+line3
                    #self.message=self.message+wre.sub(' ',line,1)+"\n"
            else:
                fo=fo+1
                line3=wre.sub(' ',line,1)[0:60]+'\n'
                line2=wre.sub(' ',line,1)[60:]
                while(len(line2)>0):
                    line3=line3+'    '+line2[0:55]+'\n'
                    line2=line2[55:]
                opt=opt+line3
        if self.f==0: self.message=self.message+'    none\n'
        if fo>0:
            self.message=self.message+opt

        
    def abort(self,error=''):
        os.system("clear")
        if error: print error
        print self.message
        sys.exit()

    def parse(self,initvar,args):
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
        
#self-tester
if __name__ == '__main__':
    inpHand('Usage: myfirst.py -i -o [-p -x]',
            '  -i _AR_input      listfile'    ,
            '  -o _AR_outdir   output dir'   ,
            '  -p __optdir  optdir ja ja ja',
            '  -x __optflag  optional flag'   ).parse(locals(),sys.argv)

    print (input,outdir,optdir,optflag)
    sys.exit(0)
