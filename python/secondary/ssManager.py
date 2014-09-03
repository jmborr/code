#!/usr/bin/python
"""secondary structure object"""

import os,sys,re
from copy import deepcopy
from secondary.seqdat2ss import translateSec,nn2ss

minChunkL_0={'coil':1, 'helix':4, 'strand':2} #mimimal chunk size

class ssec:
    ######################################################################
    def __init__(self,seqdat='',sseq='',dict={},ss2=None):
        self.dict={'coil':'C','helix':'H','strand':'E'}
        self.sseq=''
        self.error=''
        self.chunks={}
        self.minChunkL=minChunkL_0
        self.content={}
        self.fracontent={}
        if seqdat:
            self.sseq=translateSec(seqdat,translator=nn2ss)
        elif sseq and dict:
            self.dict=dict
            self.sseq=sseq
        elif ss2:
            self=deepcopy(ss2)
        self.L=len(self.sseq)

    #######################################################################
    def initChunks(self,minChunkL=''):
        """for each secondary type, store a list of chunks of valid size

        Optional arguments:
        minlengths='' dictionary with minimal chunk sizes. Will use default
                      self.minChunkL if not passed
        """
        if not minChunkL: minChunkL=self.minChunkL
        for type in self.dict.keys():
            minL=minChunkL[type]
            self.chunks[type]=[]
            chunks=re.findall(self.dict[type]+'+',self.sseq)
            for chunk in chunks:
                if len(chunk) > minL: self.chunks[type].append(chunk)

    ######################################################################
    def getFreshChunks(self,type,minL):
        """return number of chunks of a particular secondary type
        without taking into account self.chunks

        Required arguments:
        type, a valid secondary type
        minL, a minimum chunk size
        """
        if type not in self.dict.keys():
            self.error='ERROR in ssManager::ss::getFreshChunks. '+\
                        type+' not a valid secondary type'
            return False
        list=re.findall(self.dict[type]+'+',self.sseq)
        chunks=[]
        for chunk in list:
            if len(chunk)>=minL: chunks.append(chunk)
        return chunks
    
    ######################################################################
    def numberOfChunks(self,type,minL=0):
        """ return number of chunks of particular type

        Required arguments:
        type: a valid secondary type
        
        Optional arguments:
        minL=0 minimal required chunk size. If not used, will use self.minChunkL[type]
        """
        if type not in self.dict.keys():
            self.error='ERROR in ssManager::ss::number. '+type+' not a valid secondary type'
            return -1
        if minL==0: minL=self.minChunkL[type]
        if minL<self.minChunkL[type]:
            return len( self.getFreshChunks(type,minL) )
        if not self.chunks: self.initChunks()
        chunks=self.chunks[type]
        nchunks=0
        for chunk in chunks:
            if len(chunk)>=minL: nchunks+=1
        return nchunks

    ######################################################################
    def initContent(self):
        """initalize the fractional content of every secondary type"""
        for type in self.dict.keys():
            self.content[type]=len( re.findall(self.dict[type],self.sseq) )
        return True


    ######################################################################
    def initFractionalContent(self):
        """initialize content and fractional content of every secondary type"""
        self.initContent()
        for type in self.dict.keys():
            self.fracontent[type]=(1.0*self.content[type])/self.L
        return True
    
