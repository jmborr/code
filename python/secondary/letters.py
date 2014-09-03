#!/usr/bin/python

import re

symbolstr0='alpha:H beta:E coil:C gap:-'
ssdict0={'alpha':'H', 'beta':'E', 'coil':'C', 'gap':'-'}
sstypes=ssdict0.keys(); sstypes.sort(); #list of alphabetically sorted sec. str. types
l=ssdict0.values(); l.sort(); valid_ss0=''.join(l) #string  of alphabetically sorted symbols

sblseqdat='alpha:2 beta:4 coil:1 gap:-'

#diverse regular expressions
res={}
res['format']=re.compile('[a-z]+\:[0-9,A-Z,-]') #required format for a pair like 'alpha:H'

#   ***********************************************

def missingSymbol(symbolstr):
    """check that the symbol string contains all secondary structure types"""
    for type in sstypes:
        if type not in symbolstr: return True
    return False

#   ***********************************************

def failsFormat(symbolstr):
    """check that the symbol string complies with required format

    Format like this: 'alpha:H beta:E coil:C gap:-'"""
    #check all needed secondary structure types are present
    if missingSymbol(symbolstr): return True
    #check every pair complies to required format
    pairs=symbolstr.split()
    rex=res['format']
    for pair in pairs:
        if not rex.match(pair): return True        
    return False
    
#   ***********************************************

def returnDict(symbolstr):
    if failsFormat(symbolstr): return {}
    ssdict={}
    """return a dictionary of symbols

    Required format like this: 'alpha:H beta:E coil:C gap:-'"""
    pairs=symbolstr.split()
    for pair in pairs:
        for type in sstypes:
            if type in pair:
                key,val=pair.split(':')
                ssdict[key]=val #last character of pair is the symbol

    return ssdict

#   ***********************************************

def getValidSS(symbolstr):
    """return the one-letter symbols in a single string, and sorted"""
    v=returnDict(symbolstr).values()
    v.sort()
    return v

#   ***********************************************

def transDict(symbolstr1,symbolstr2):
    """provide a dictionary that link corresponding symbols"""
    td={}
    d1=returnDict(symbolstr1)
    if symbolstr2!=symbolstr1:
        d2=returnDict(symbolstr2)
    else:
        d2=d1
    for sstype in sstypes: td[ d1[sstype] ] = d2[sstype]
    return td

#   ***********************************************

def translate(seq,symbolstr1,symbolstr2):
    """translate one sequence from one set of symbols to the other"""
    m=transDict(symbolstr1,symbolstr2)
    seq2=''
    for x in seq: seq2+=m[x]
    return seq2

#   ***********************************************

def symbolstrTOnn2ss(symbolstr):
    """ translate symbolstr to nn2ss as used in seqdat.py"""
    td=transDict(sblseqdat,symbolstr)
    del td['-'] #remove gap symbol in dictionary
    return td 
