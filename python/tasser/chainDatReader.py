#!/usr/bin/python

import os,sys,re
from utilities.small_utilities import chomp,junkName,Bye
from inputArgs.inputArgs import addAbsPath,inpHand
from random import randint

res={}
res['TMalign']=re.compile('TM-score=(\d\.\d+)')
res['TMscore']=re.compile('TM-score\s+=\s(\d\.\d+)')


def buildScaffold(natf,fastaf=''):
    "return all ATOM ... CA lines in natf, but without coordinates"
    #first line empty because we want indexes in native file to coincide with indexes of list
    scaf=['',]
    if fastaf:
        from seq.fastaManager import importFastaEntry
        from seq.letters import one2three
        header,seq=importFastaEntry( open(fastaf,'r') )
        n=1
        for x in seq:
            scaf.append('ATOM    %3d  CA  %s _ %3d    '%(n,one2three[x],n))
            n+=1
        return scaf
    if natf:
        flag=1
        pt=open(natf,'r')
        for l in pt.readlines():
            if l[0:4]=='ATOM' and l[12:16]==' CA ':
                if flag: #check integrity of native pdb file
                    if l[10]!='1':
                        sys.stderr('ERROR from buildScaffold in chainReader\nWe need a native file with onlyCA atoms,\nand the atom number of first CA must be "1"')
                        sys.exit(1)
                    flag=0
                scaf.append(l[0:30])
        pt.close()
        return scaf


class chainDatReader:
    """
    This class has utilities to read a chain.dat
    Attributes:
     chainf: absolute path to chain.dat filename (default=./chain.dat)
     natf:   absolute path to native pdb file (default=./CA)
     L:      sequence length
     nt:     number of templates inside chain.dat
     type:   type of target ('easy','medm',hard')
     scaf:   list with all ATOM ... CA lines in natf, but without coordinates
     templ:  list template coordinates in  pdb stile, helped by scaf
     cov  :  coverage of current templ
    """
        
    def __init__(self,chainf='chain.dat',natf='',fastaf=''):
        self.chainf=addAbsPath(chainf)
        self.pt=open(chainf,'r')
        [self.nt,self.type]=self.pt.readline().split()  ;  self.nt=int(self.nt)
        self.natf=''
        if natf: self.natf=addAbsPath(natf)
        #Bye('%s\n%s\n%s\n%s'%(self.chainf,self.nt,self.type,self.natf))
        self.scaf=buildScaffold(self.natf,fastaf=fastaf)
        self.L=len(self.scaf)-1
        self.template_header=''
        self.templ=self.scaf[:] #shallow copy
        self.chaindat=[]
        self.cov=0


    def readTempl(self):
        """read another template"""
        items=self.pt.readline().split()
        if not items: return 0
        nitems=len(items)
        n=int(items[0])
        if nitems>2: self.template_header=items[2]
        self.cov=float(n)/self.L
        self.templ=self.scaf[:] #initialize the template
        rec=re.compile('\s*(\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)')
        self.chaindat=[]
        for i in range(n):
            line=chomp(self.pt.readline()) #;print line
            match=rec.match(line)
            self.chaindat.append(line.strip())
            if match:
                [m,x,y,z]=match.groups() #; print [m,x,y,z] ;sys.exit(1)
                m=int(m) ; x=float(x)  ; y=float(y); z=float(z)
            self.templ[m]=self.templ[m][0:30]+'%8.3f%8.3f%8.3f'%(x,y,z)
        return 1


    def templToChaindat(self,fname='chain.dat'):
        """ write self.chain to file"""
        buffer='1 '+self.type+'\n'+`len(self.chaindat)`+'\n'+'\n'.join(self.chaindat)
        open(fname,'w').write(buffer)
        

    def TMscore(self,struct=''):
        """calculate TM-score of current template"""
        if not struct: struct=self.natf
        jptf=self.templToFile()
        cmd='TMscore '+jptf+' '+struct+'|grep -e "TM-score    ="'
        match=res['TMscore'].search(os.popen(cmd).readline())
        tm=0.0
        if match: tm=float( match.group(1) )
        os.system('/bin/rm '+jptf)
        return tm
        
    def templToFile(self,fname=''):
        buffer=''
        for line in self.templ:
            if len(line) > 30: buffer+=line+'\n'
        if fname: jptf=fname
        else:     jptf=junkName()
        open(jptf,'w').write(buffer)
        return jptf

        
    def TMalign(self,struct=''):
        """calculate TM-score of current template"""
        if not struct: struct=self.natf
        jptf=self.templToFile()
        cmd='TMalign '+jptf+' '+struct+'|grep -e "TM-score="'
        match=res['TMalign'].search(os.popen(cmd).readline())
        if match: tm=float( match.group(1) )
        os.system('/bin/rm '+jptf)
        return tm

        
if __name__=='__main__':
    
    ih=inpHand('Usage: chainDatReader.py [options]',
               ' -a __natf native file name (default="./CA")',
               ' -b __chaf chain.dat type file (default="./chain.dat")')
    ih.parse(locals(),sys.argv)
   
    if not natf: natf='./CA'
    if not chaf: chaf='./chain.dat'

    ch=chainDatReader(chainf=chaf,natf=natf)
    while(ch.readTempl()):
        ch.TMalign()
