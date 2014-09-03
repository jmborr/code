import re,os,sys

three2one={'ALA':'A','VAL':'V','PHE':'F','ILE':'I','LEU':'L','PRO':'P',
           'MET':'M','ASP':'D','GLU':'E','LYS':'K','ARG':'R','SER':'S',
           'THR':'T','TYR':'Y','HIS':'H','CYS':'C','ASN':'N','GLN':'Q',
           'TRP':'W','GLY':'G',}
one2three={'A':'ALA','V':'VAL','F':'PHE','I':'ILE','L':'LEU','P':'PRO',
           'M':'MET','D':'ASP','E':'GLU','K':'LYS','R':'ARG','S':'SER',
           'T':'THR','Y':'TYR','H':'HIS','C':'CYS','N':'ASN','Q':'GLN',
           'W':'TRP','G':'GLY',}
def listToDict(list): #return a dictionary out of a list
    x=len(list)-1
    i=0
    dic={}
    while i<x:
        dic[list[i]]=list[i+1]
        i=i+2
    return dic

def chomp(object):
    t=`type(object)` ; #print t # string representation of type
    if re.compile('str').search(t):
        if object:
            if(object[-1]=='\n'):  return object[:-1]
        return object
    elif re.compile('list').search(t):
        for i in range(len(object)):
            if object[i]:
                if(object[i][-1]=='\n'):  object[i]=object[i][:-1]
        #for line in object: no good, because line isn't a reference to an element of object
            #if(line[-1]=='\n'):
                #line=line[:-1]
        return None
    elif re.compile('tuple').search(t):
        stderr.write("ERROR: chomp can't change in place a tuple !")

def replace_at(string,position,replacement):
    l=len(string)-1 #last index
    if position==0:
        if l>0: return replacement+string[1:]
        else  : return replacement
    if position==l: return string[0:l]+replacement
    return string[0:l-1]+replacement+string[l:]

