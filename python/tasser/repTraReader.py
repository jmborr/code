#!/usr/bin/python

import os,sys,re
from utilities.small_utilities import chomp,junkName,Bye
from inputArgs.inputArgs import addAbsPath,inpHand
from random import randint
from tmScore.tm_score import gettmCA

res={}
res['L=']=re.compile('L=(\d+)')
res['old L']=re.compile('\s+(\d+)\s')

###############   auxiliary functions ########################

def readL(ptr):
    pos=ptr.tell()
    header=ptr.readline()
    ptr.seek(pos)
    L=0
    m=res['L='].search(header)
    if not m:
        m=res['old L'].match(header)
        if not m: return 0
    return int(m.group(1))
##############################################################


    
##############################################################
class repTraSnapshot:

    def __init__(self,ptr,L=''):
        '''
        ptr: pointer to trajectory file
        '''
        self.error=''
        self.header=''
        self.xyzs=[]
        self.L=L
        
        #read header, determine L if not passed
        if not L:
            self.L=readL(ptr)
            if not self.L:
                self.error+='repTraSnapshot::__init__: header not contain number of residues\n'
                return
        
        #read header and line coordinates
        self.header=chomp(ptr.readline())
        if not self.header:
            self.error+='repTraSnapshot::__init__: reached end of file'
            return
        for i in range(self.L):
            self.xyzs.append( chomp(ptr.readline()) )

    ##########################################################
            
    def writepdb(self,outf,pdbf,extraHeader='',extraFooter=''):
        '''
        write snapshot to a file, passing a pdb file as template
        '''
            
        CAs=os.popen('grep "ATOM" '+pdbf+' | grep " CA "').readlines() #CA atom lines
        if len(CAs)!=self.L:
            self.error+='repTraSnapshot::write: '+pdbf+' and snapshot contain different number of residues\n'
            return
            
        buf=extraHeader
        for i in range(self.L):
            CA=CAs[i]
            x,y,z=self.xyzs[i].split()
            x=float(x) ; y=float(y) ; z=float(z)
            pdbxyz='%8.3lf%8.3lf%8.3lf'%(x,y,z)
            buf+=CA[0:30]+pdbxyz+CA[54:]
        buf+=extraFooter
        open(outf,'w').write(buf)
            
    ##########################################################
            
    def write2string(self):
        return self.header+'\n'+'\n'.join(self.xyzs)+'\n'
            
    ##########################################################

    def addTMscore(self,strf):
        junk=junkName()
        self.writepdb(junk,strf) #write coords to temporary file
        tm=gettmCA(junk,strf)    #calculate tm
        self.header+=' TM=%5.3lf'%(tm)
        os.system('/bin/rm '+junk)

##################################################################

    
class repTraReader:
    """
    This class has utilities to read a rep*tra files
    """

    def __init__(self,reptraf='./rep1.tra'):
        self.error=''
        self.reptraf=addAbsPath(reptraf)
        self.natf=''
        self.reptrap=open(reptraf,'r') #pointer to file
        self.L=readL(self.reptrap)
        if not self.L: self.error='repTraReader::__init__: L not successfully read'
        self.lastshot=None             #last snapshot read
            
            
    ##################################################################        
        
    def getNatFileName(self,natf):
        natf=addAbsPath(natf)
        if not os.path.exists(natf):
            self.error+='ERROR from repTraReader::readNatFile: '+natf+' does not exists\n'
        self.natf=natf

    ##################################################################        

    def readNextSnapShot(self):
        self.lastshot=repTraSnapshot(self.reptrap)
        if self.lastshot.error: return False
        return True
        
    ##################################################################        
        
    def addTMscore(self,strf='',outf='_SAME_'):
        '''
        addTMscore(self,strf=self.natf,outf=\'_SAME_\')
        For every snapshot, calculate TMscore to strf and store in trajectory file outf
        '''

        #prepare pointer to output file
        ptw=None
        rewrite=False
        if outf=='_SAME_': #overwrite the rep*tra file
            rewrite=True
            outf=junkName()
        ptw=open(outf,'w')

        #check existence of strf
        if not strf:
            strf=self.natf
            if not strf:
                self.error+='ERROR from repTraReader::addTMscore: No structure passed, and native file has not been loaded\n'
                return
        
        self.reptrap.seek(0) #rewind pointer to beginning of file        
        self.lastshot=repTraSnapshot(self.reptrap,L=self.L) #read next snapshot
        n=1
        while not self.lastshot.error:
            self.lastshot.addTMscore(strf)             #add TMscore
            ptw.write( self.lastshot.write2string() )      #write output
            self.lastshot=repTraSnapshot(self.reptrap,L=self.L)
            n+=1

        if rewrite:
            self.reptrap.close()
            os.system('/bin/mv '+outf+' '+self.reptraf)
            self.reptrap=open(self.reptraf,'r')
        
    ##################################################################        

    def write2pdb(self,dir='./tmp',prefix='a',strf=''):
        
        if not strf:
            strf=self.natf
            if not strf:
                self.error+='ERROR from repTraReader::write2pdb: No structure passed, and native file has not been loaded\n'
                return
            
        if not os.path.exists(dir): os.system('/bin/mkdir -p '+dir)
        
        index=0
        self.reptrap.seek(0)
        self.lastshot=repTraSnapshot(self.reptrap,L=self.L)
        while not self.lastshot.error:
            outf=os.path.join(dir,prefix+'%05d'%(index))
            self.lastshot.writepdb(outf,strf)
            self.lastshot=repTraSnapshot(self.reptrap,L=self.L)
            index+=1
