#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp,Bye
from seq.letters import valid_A
from copy import deepcopy

rs={}
rs['size1']=re.compile('Length of sequence 1:\s+(\d+)')
rs['size2']=re.compile('Length of sequence 2:\s+(\d+)')
rs['algnL']=re.compile('Aligned length:\s+(\d+)')
rs['idL'  ]=re.compile('Identical length:\s+(\d+)')
rs['seqid']=re.compile('Sequence identity:\s+(\d\.\d+)')
rs['seq1' ]=re.compile('Sequence identity:[^\n]+\n+([^\n]+)\n')
rs['seq2' ]=re.compile('Sequence identity:[^\n]+\n+[^\n]+\n[^\n]+\n([^\n]+)')

def extract(key,all):
    match=rs[key].search(all)
    if match: return match.group(1)
    else: return ''

#sequence containing gaps
class gappedSeq:

    def __init__(self,gs,gapsymbol='-'):
        self.gs=gs
        self.gsL=len(gs)
        self.s=gs.replace(gapsymbol,'')
        self.sL=len(self.s) #sequence without gaps
        self.gapsymbol=gapsymbol
        #correspondence between gapped and ungapped indexes
        self.gsi2si=[-1]*self.gsL #initialize to non-correspondence
        self.si2gsi=[-1]*self.sL  #initialize to non-correspondence
        n=0
        for m in range(0,self.gsL):
            if self.gs[m]!=self.gapsymbol: #position occupied by one amino acid
                self.gsi2si[m]=n #from gapped index to ungapped index
                self.si2gsi[n]=m #from ungapped index to gapped index
                n+=1

    def __iadd__(self,extrags):

        """append a gapped sequence"""
        
        gs2=extrags.gs
        oldsymbol=extrags.gapsymbol
        newsymbol=self.gapsymbol
        if oldsymbol != newsymbol: gs2=gs2.replace(oldsymbol,newsymbol)
        self=gappedSeq(self.gs+gs2)
        return self


    def __add__(self,extrags,newsymbol=''):

        """new gapped sequence from the merging of two"""

        if not newsymbol: newsymbol=self.gapsymbol

        gs1=self.gs
        if self.gapsymbol != newsymbol:
            gs1=gs1.replace(self.gapsymbol,newsymbol)
        gs2=extrags.gs
        if extrags.gapsymbol != newsymbol:
            gs2=gs2.replace(extrags.gapsymbol,newsymbol)
        return gappedSeq(gs1+gs2,gapsymbol=newsymbol)


    def __getslice__(self,b,e):

        """new gapped sequence extracted from chunk"""
        
        gs=self.gs[b:e]
        return gappedSeq(gs)


    def __getitem__(self,i):
        return self.gs[i]

    def __repr__(self):
        return self.gs
    
    def extractGaps(self):
        """make a list of all the gaps segments"""
        return re.compile(self.gapsymbol+'+').findall(self.gs)

    def gapIndexesGappedSeq(self):
        """return a list of begin and end positions for each gap. Example
        [[4,9],[56,89]..]. The reference is the gapped sequence"""
        list=[]
        for match in re.compile(self.gapsymbol+'+').finditer(self.gs): #iterate over all matches
            bounds=match.span() ; bounds[1]=bounds[1]-1
            list.append(bounds)
        return list

    #return a list of positions where gap begins. The reference is the
    #ungapped sequence. Thus "n" means a gap begins right after "n" residues    
    def gapPositionsUngappedSeq(self):
        pos=[]
        for match in re.compile(self.gapsymbol+'+').finditer(self.gs):
            x=match.start() #;print 'x=',x
            if x==0: pos.append(-1)
            else:    pos.append(self.gsi2si[x-1])
        return pos 


    def changeGapSymbol(self,newsymbol):

        """change the gap symbol"""

        oldsymbol=self.gapsymbol
        self.gapsymbol=newsymbol
        if oldsymbol != newsymbol: self.gs=self.gs.replace(oldsymbol,newsymbol)
        
        
#alignment of to gapped sequences. It can be a chunk of a parent alignment
class alignment:

    
    def __init__(self,gs1,gs2,gapsymbol=''):
        """
        gs1:first gapped sequence, an object of class gappedSeq
        """

        if not gapsymbol: gapsymbol=gs1.gapsymbol
        self.gapsymbol=gapsymbol        
        self.gs={}
        self.gs[1]=deepcopy(gs1)
        self.gs[1].changeGapSymbol(gapsymbol)
        self.gs[2]=deepcopy(gs2)
        self.gs[2].changeGapSymbol(gapsymbol)
        self.l=0         #length of the chunk
        self.parent=None #the aligned chunk does not come from other aligned chunk
        self.b=0         #b is the index in the parent chunk where the current chunk begins
        if gs1==None or gs2==None: return None
        if gs1.gsL != gs2.gsL:
            sys.stderr.write('ERROR from alignment.__init__: gapped sequences ')
            sys.stderr.write('have different lengths\n')
            return None
        #derived data
        self.l=gs1.gsL #;print self.l

    def changeGapSymbol(self,newsymbol):

        """change the gapsymbol in the alignment"""
        
        self.gapsymbol=newsymbol
        self.gs[1].changeGapSymbol(newsymbol)
        self.gs[2].changeGapSymbol(newsymbol)

    
    def extract(self,b,e):

        """extract a portion and return as new alignment"""

        alg=alignment(self.gs[1][b:e+1],self.gs[2][b:e+1])
        alg.parent=self
        return alg


    def __getslice__(self,b,e):
        return self.extract(b,e-1) #notice the shift in the last index


    def extractII(self,b,e,i=1,appendGaps=False,prependGaps=False):

        """return a chunk of the alignment that begins in ungapped index 'b'
        and finish in ungapped index 'e' for sequence i

        extractII(self,b,e,i=1,appendGaps=False,prependGaps=False)
        appendGaps example:
           LALALA----LILILI
           LALALALALALILILI
        we want to extract from 0 to 6 for sequence 1. But then, which is correct?
           LALALA  or  LALALA----
           LALALA      LALALALALA
        both alignments have first and last ungapped indexes equal to 0 and 6, respectively.
        The default behaviour is the first alignment, but we can include the gaps if we set
        prepenGaps to True
        Similar thing for prepend gaps        
        """

        bg=self.gs[i].si2gsi[b] #switch to gapped indexes
        if prependGaps:
            if b>0: bg=self.gs[i].si2gsi[b-1]+1
            else: bg=0
        eg=self.gs[i].si2gsi[e]
        if appendGaps:
            lastUngappedIndex=self.gs[i].sL-1
            if e<lastUngappedIndex: eg=self.gs[i].si2gsi[e+1]-1
            else: eg=self.gs[i].gsL-1
        return self.extract(bg,eg)

        
    def __iadd__(self,alg):

        """append an alignment"""

        gs1=alg.gs[1].gs     
        gs2=alg.gs[2].gs
        if alg.gapsymbol != self.gapsymbol: #armonize gap symbols
            gs1=gs1.replace(alg.gapsymbol,self.gapsymbol)
            gs2=gs2.replace(alg.gapsymbol,self.gapsymbol)
        gs1=self.gs[1].gs+alg.gs[1].gs     
        gs2=self.gs[2].gs+alg.gs[2].gs
        self=alignment( gappedSeq(gs1),gappedSeq(gs2))
        return self
        

    def merge(self,alg,newsymbol=''):

        """new alignment by merging two

        we allow for changes in the gapsymbol"""

        if not newsymbol: newsymbol=self.gapsymbol

        gs11=self.gs[1].gs     
        gs12=self.gs[2].gs
        oldsymbol=self.gapsymbol
        if oldsymbol != newsymbol: #armonize gap symbols
            gs11=gs11.replace(oldsymbol,newsymbol)
            gs12=gs12.replace(oldsymbol,newsymbol)

        gs21=alg.gs[1].gs     
        gs22=alg.gs[2].gs
        oldsymbol=alg.gapsymbol
        if oldsymbol != newsymbol: #armonize gap symbols
            gs21=gs21.replace(oldsymbol,newsymbol)
            gs22=gs22.replace(oldsymbol,newsymbol)
        
        gs1=gs11+gs21
        gs2=gs12+gs22
        self=alignment( gappedSeq(gs1),gappedSeq(gs2),gapsymbol=newsymbol)
        
        
    def __add__(self,alg):

        """new alignment by merging two"""
        
        return self.merge(alg)
        
        
    def addSplicingEvent(self,b,e,variation='',i=1):

        """insert an splicing event

        we chop the alignment in three chunks, modify the internal
        alignment, then merge the chunks"""

        #extract alignment prior to the splicing event
        algprev=None
        if b>0:
            algprev=self.extractII(0,b-1,i,prependGaps=True,appendGaps=True)
        else:
            algprev=alignment(gappedSeq(''),gappedSeq(''))
        L=self.gs[i].sL-1 #last index in the alignment
        if e<L:
            algnext=self.extractII(e+1,L,i,prependGaps=True,appendGaps=True)
        else:
            algnext=alignment(gappedSeq(''),gappedSeq(''))
        #splice the remaining chunk of the alignment
        alg=self.extractII(b,e,i)        
        gap=alg.gapsymbol
        seq1=alg.gs[1].s
        seq2=alg.gs[2].s
        if i==1:
            gs1=seq1+gap*len(variation) #seq1---------
            gs2=gap*len(seq1)+variation #----variation (substitue seq2 with variation)
        else:
            gs1=variation+gap*len(seq2) #variation---- (substitute seq1 with variation)
            gs2=gap*len(variation)+seq2 #---------seq2
        alg=alignment( gappedSeq(gs1), gappedSeq(gs2) )
        #fuse the three alignments into one
        return mergeAlignments( [algprev,alg,algnext] )
        
    
    #return indexes where gapped sequence "i" begins and ends
    def boundingIndexes(self,i=1):
        gs=self.gs[i]  ;  gaps=self.gapsymbol  #;  print 'gs=\n',gs
        n=0
        while(1):
            if gs[n]!=gaps: break
            n+=1
        begin=n
        n=gs.gsL-1
        while(1):
            if gs[n]!=gaps: break
            n-=1
        end=n
        return begin,end

       
    def identitiesLine(self,idsymbol=':',valid_aa=valid_A):
        gs1=self.gs[1].gs ; gs2=self.gs[2].gs
        idline=''
        for i in range(0,self.l):
            if gs1[i]==gs2[i]:
                if gs1[i] in valid_aa and gs2[i] in valid_aa: #identities must be known
                    idline+=idsymbol
            else: idline+=' '
        return idline

    #return number of identical residues
    def identicalL(self,idsymbol=':'):
        full=self.identitiesLine(idsymbol=idsymbol)       
        return len(full)-len(full.replace(idsymbol,''))


    def alignedLine(self,algsymbol=':',valid_aa=valid_A):
        """gapped line marking aligned fragments"""
        gs1=self.gs[1].gs ; gs2=self.gs[2].gs
        idline=''
        for i in range(0,self.l):
            if gs1[i] in valid_aa and gs2[i] in valid_aa:
                idline+=algsymbol
            else: idline+=' '
        return idline

    def nAligned(self,valid_aa=valid_A):
        """number of aligned residues"""
        return len( self.alignedLine(valid_aa=valid_aa).replace(' ','') )
        
    def presentAlignedFragments(self,algsymbol=':',valid_aa=valid_A):
        l=self.alignedLine(algsymbol=algsymbol,valid_aa=valid_aa)
        return re.compile(algsymbol+'+').findall( l )

    
    def isAlignedLines(self,idsymbol=':',valid_aa=valid_A, algtype='identical'):
        """return two strings, one for each sequence, with 1 for
        aligned residue at that position in the ungapped sequence, 0
        otherwise"""
        #init ll[1] as a list of 0's with length equal to length of
        #ungapped first sequence
        ll={1:['0']*self.gs[1].sL, 2:['0']*self.gs[2].sL}
        n=0 #index of gapped sequence
        if algtype=='identical':
            line=self.identitiesLine(idsymbol=idsymbol,valid_aa=valid_aa)
        elif algtype=='aligned':
            line=self.alignedLine(algsymbol=idsymbol,valid_aa=valid_aa)
        for x in line:
            if x==idsymbol:
                #if position corresponds to two aligned residues, find
                #in the gapped sequence "id" the corresponding index in
                #its ungapped sequence, and mark it as aligned                
                for id in ll.keys(): ll[id][ self.gs[id].gsi2si[n] ]='1'
            n+=1
        for id in ll.keys(): ll[id]=''.join(ll[id]) #convert list to string
        return ll

    #make a list of all the gaps segments
    def extractGaps(self,gapsymbol='-'):
        gapList=re.compile(gapsymbol+'+').findall(self.gs[1].gs)
        gapList+=re.compile(gapsymbol+'+').findall(self.gs[2].gs)
        return gapList


    def extractInternalGaps(self,gapsymbol='-'):
        pattern='^'+gapsymbol+'*([A-Z\-]+[A-Z])'+gapsymbol+'*$' #exclude beginning and ending gaps
        m=re.compile(pattern).search(self.gs[1].gs)
        gapList=re.compile(gapsymbol+'+').findall(m.group(1))
        m=re.compile(pattern).match(self.gs[2].gs)
        gapList+=re.compile(gapsymbol+'+').findall(m.group(1)) #; print gapList
        return gapList
        
    #extract gap segment corresponding to begginiing of unaligned region
    def Ngap(self,gapsymbol='-'):
        gap=re.compile('('+gapsymbol+'*)').match(self.gs[1].gs).group(1)
        if not gap: gap=re.compile('('+gapsymbol+'*)').match(self.gs[2].gs).group(1)
        return gap

    #extract gap segment corresponding to ending unaligned region
    def Cgap(self,gapsymbol='-'):
        gap=re.compile('('+gapsymbol+'*)$').search(self.gs[1].gs).group(1)
        if not gap: gap=re.compile('('+gapsymbol+'*)$').search(self.gs[2].gs).group(1)
        return gap
        
    #return biggest gap
    def extractBiggestGap(self,gapsymbol='-'):
        n=0 ; biggestGap=''
        for gap in self.extractGaps(gapsymbol):
            if len(gap)>n:
                n=len(gap)
                biggestGap=gap
        return biggestGap

    #return biggest internal gap
    def extractBiggestInternalGap(self,gapsymbol='-'):
        n=0 ; biggestInternalGap=''
        for gap in self.extractInternalGaps(gapsymbol):
            if len(gap)>n:
                n=len(gap)
                biggestInternalGap=gap
        return biggestInternalGap

    #return number of gaps among the two gapped sequences
    def ngaps(self,gapsymbol='-'):
        return len( ''.join(self.extractGaps(gapsymbol)) )
        
    #return number of aligned residues
    def alignedL(self):
        return self.l-self.ngaps()

    #return number of aligned residues if their identities are known
    def alignedL20(self,valid_aa=valid_A):
        nalg=0
        gs1=self.gs[1].gs ; gs2=self.gs[1].gs
        for i in range(0,self.l):
            if gs1[i] in valid_aa and gs2[i] in valid_aa: nalg+=1
        return nalg
    
    #return coverage
    def coverage(self):
        return (1.0*self.alignedL())/self.l
    
    #return sequence identity over aligned region
    def seqIdentity(self):
        return (1.0*self.identicalL())/self.alignedL()

    #return the projection of one sequence onto another. For instance, we have
    #seq1: 1234---89
    #seq2: 12-45678-
    #the projection of seq2 onto seq1 is 12-48-, and seq1 onto seq2 is 124---8
    #project gs[index1] onto gs[index2]
    def project(self,index1,index2):
        gs=''  ;  gs1=self.gs[index1]  ;  gs2=self.gs[index2]
        for i in range(0,self.l):
            if gs2[i]==self.gapsymbol: continue
            gs+=gs1[i]
        return gappedSeq(gs)
        
    def info(self):
        algL=self.alignedL20()  #;print 'algL=',algL#;sys.exit(1)
        idl=self.identicalL() #;print 'idl=',idl
        seqid=0
        if algL: seqid=(1.*idl)/algL
        all= 'Length of sequence 1:%5d -> no-name\n'%(self.gs[1].sL)
        all+='Length of sequence 2:%5d -> no-name\n'%(self.gs[2].sL)
        all+='Aligned length:%5d\nIdentical length:%5d\n'%(algL,idl)
        all+='Sequence identity: %5.3lf (=%5d/%5d)\n'%(seqid,idl,algL)
        all+=self.gs[1].gs+'\n'+self.identitiesLine()+'\n'+self.gs[2].gs
        return all
        
    def __repr__(self):
        return self.info()
    
    def extractIdenticalChunks(self,minSize=1,idsymbol=':'):
        idChunks=[]
        il=self.identitiesLine()
        for match in re.compile(idsymbol+'+').finditer(il): #iterate over all matches
            bounds=match.span() #beginning and ending+1 indexes of the match
            if bounds[1]-bounds[0] >= minSize:
                idChunks.append( self.extract(bounds[0],bounds[1]-1) )
        return idChunks

    #if a parent is defined, return for each of the two
    #self.gapped-sequences-objects, the begining and ending index of the
    #parent.ungapped-sequences
    def parentUngappedSpans(self):
        parent=self.parent
        gp=parent.gapsymbol #symbol used for gaps in the gapped sequence
        if not parent: return None
        spans={}
        for i in (1,2): #we have to gapped sequences in the
            #prev is the chunk of the parent.gs[i] gapped sequence object coming before self.gs[i]
            prev=parent.gs[i][0:self.b] 
            begin=prev.sL #length of the ungapped sequence
            end=begin+self.gs[i].sL-1
            spans[i]=[begin,end]
        return spans

    #   *************************************
    def listIndexPairsIdentical(self):
        """list of pair of ungapped-indexes for identical residues

        [ (i1,j1), (i2,j2), ...,(in,jn)] where residue in i1 of seq1 identical
        to residue j1 of seq2. Indexes correspond to ungapped sequences."""
        pairs=[]
        idl=self.identitiesLine() ;
        gs1=self.gs[1] ; gs2=self.gs[2]
        for i in range(self.l):
            if idl[i]==':': pairs.append( (gs1.gsi2si[i],gs2.gsi2si[i]) )
        return pairs


def mergeAlignments( list, newsymbol='' ):

    """merge a list of alignments into a new one"""
    
    if not newsymbol: newsymbol=list[0].gapsymbol
    algFinal=alignment( gappedSeq(''),gappedSeq(''),gapsymbol=newsymbol)
    for alg in list:
        algFinal+=alg
    return algFinal


#####################################################################

def loadAlg(inpf):
    alg=[]
    if isinstance(inpf,file):
        #search for the first line of the entry
        while(1):
            line=inpf.readline()
            if not line: return [] #reached end of file
            if line.find('Length of sequence 1')>=0:
                alg.append(line) ; break
        #keep on reading until we find other entry or reach end of file
        while(1):
            pos=inpf.tell()
            line=inpf.readline()
            if not line: break #end-of-file reached
            if line.find('Length of sequence 1')>=0:
                inpf.seek(pos) #rewind to beginning of the newly found entry
                break
            alg.append(line)
    return alg


class alignOut:
    def __init__(self,inpf):
        '''inpf can be the following:
        (1) handle to input file containing results of some Yang\'s align program
        (2) string containing name of previous file
        (3) list containing the lines of the previous file
        '''
        self.inpf=None
        self.all=None
        self.size={}
        self.readError=''
        self.alg=None
        
        #type checking
        if isinstance(inpf,str): #name of the alignment file
            self.inpf=inpf
            lines=open(inpf,'r').readlines()
            self.all=''.join(lines) #join all lines of the input file
        elif isinstance(inpf,list): #alignment contained in a list
            for i in range(0,len(inpf)):
                if inpf[i][-1] != '\n': inpf[i]+='\n' #we need these end-of-line characters
            self.all=''.join(inpf)
        elif isinstance(inpf,file): #file handle to the alignment file
            lines=loadAlg(inpf)
            self.all=''.join(lines)
        else: sys.stderr.write('ERROR: wrong input to instantiate an alignOut object\n')

        
        #gather a few facts of the aligment
        self.size[1]=extract('size1',self.all) #size of first sequence
        if self.size[1]=='':
            self.readError+='ERROR (alignOut:__init__:  reading size of first sequence\n' 
            
        self.size[2]=extract('size2',self.all) #size of second sequenc
        if self.size[2]=='':
            self.readError+='ERROR (alignOut:__init__:  reading size of second sequence\n' 

        self.alignedL=extract('algnL',self.all)#length of aligned region
        if self.alignedL=='':
            self.readError+='ERROR (alignOut:__init__:  reading size length of aligned region\n' 
       
        self.idL=extract('idL',self.all)       #length of identical region
        if self.idL=='':
            self.readError+='ERROR (alignOut:__init__:  reading size length of identical region\n' 

        self.seqid=extract('seqid',self.all)   #sequence identity of the aligned region
        if self.seqid=='':
            self.readError+='ERROR (alignOut:__init__:  reading sequence identity of the aligned region\n' 

        #gather the two sequences including gaps, dump them into an alignedChunks object
        gs1=extract('seq1',self.all)
        if gs1=='':
            self.readError+='ERROR (alignOut:__init__:  reading gs1\n' 
        
        gs2=extract('seq2',self.all)
        if gs2=='':
            self.readError+='ERROR (alignOut:__init__:  reading gs2\n' 

        if self.readError!='': return None

        gs1=gappedSeq( gs1 ) #sequence with gaps object
        gs2=gappedSeq( gs2 )
        
        self.alg=alignment(gs1,gs2)


        
