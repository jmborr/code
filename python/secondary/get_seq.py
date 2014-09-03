#!/usr/bin/python

import sys,re
from utilities.small_utilities import Bye
from letters import failsFormat,symbolstr0,valid_ss0,getValidSS,symbolstrTOnn2ss
from seqdat2ss import translateSec
from seq.fastaManager import importFastaEntry

standardInputTypes=['fasta','seq.dat']
res={}
res['header']=re.compile('>\s*(\S+)') #extract the header out of the fasta header line

#   *****************************************

def importFastaEntry(pin,firstLetterHeader=True,valid_ss=valid_ss0,symbolstr=''):
    seq='' ;  header='' ; pos=0
    #override valid_ss with symbolstr
    if symbolstr: valid_ss=getValidSS(symbolstr)
    #find beginning of fasta entry
    while not header:
        header=pin.readline()
        if not header: return ['',''] #end-of-file reached
        match=res['header'].match(header)
        if match:
            header=match.group(1)
            if firstLetterHeader: header=header.split()[0]
            break
    #read the sequence
    line=pin.readline().strip()
    while line[0] in valid_ss: #line is a chunk of the sequence
        seq+=line
        pos=pin.tell()
        line=pin.readline().strip()
    if line: pin.seek(pos) #go back one line, since we over-read in the previous loop
    
    return [header,seq]
        
#   *****************************************

def get_seq(input,symbolstr=symbolstr0,inputType=''):
    if failsFormat(symbolstr): return seq
    #Determine inputType
    if not inputType and isinstance(input,str):
        for query_type in standardInputTypes: #contains 'fasta' or 'seq.dat' strings
            if query_type in input:
                inputType=query_type
                break
    if isinstance(input,file): inputType='file'
    if not inputType: inputType='seq'
    #Obtain sequence
    seq=''
    if inputType=='seq':
        validSS=getValidSS(symbolstr) #secondary structure symbols
        for X in input:
            if X not in validSS:
                sys.stderr.write('ERROR get_seq could not determine type of '+input+'\n')
                return ''
        seq=input
    elif inputType=='fasta':
        header,seq=importFastaEntry(open(input,'r'))
    elif inputType=='seq.dat':
        seq=translateSec(open(input,'r'),translator=symbolstrTOnn2ss(symbolstr))
    elif inputType=='file':
        seq=open(input,'r').readline().strip()

    return seq    
        
