#!/usr/bin/python

"""create pickle dumps containing basic info. Create also MAC list"""

import os,sys,re
from mac import scroot,scout,filteredlist,maclist,MAC,updColletion
from utilities.codedir import scratchdir
from utilities.small_utilities import chomp,Bye

#global variables
fheaders=os.popen('grep ">" '+scroot+'/preparing/missingXray.080503.txt').readlines()#fasta headers
pat=re.compile('(\d+-\d+)')

def parse_chunks(line):
    """split the chunks in the fasta header into a nice list of lists
    instance of line: '>1A06:1 #1-9, #54-63, #164-181, #317-332'
    """
    chunks=[]
    for chunk in pat.findall(line):
        chunks.append( [int(x) for x in chunk.split('-')] )
    return chunks

#process each sequence from the filtered set
macl=[]
headers= chomp(open(filteredlist).readlines())
remaining=len(headers)
for header in headers:
    print remaining
    header2=header[0:4]+':'+header[4:]+' '
    for fheader in fheaders: #search header in the fasta headers
        if header2 in fheader:
            dir=os.path.join(scout,header)
            for chunk in parse_chunks(fheader):
                x=MAC(pdbheader=header,b_e=chunk)
                x.update('L',x.b_e[1]-x.b_e[0]+1)
                x.pickleDump(dir) #also generates macheader attribute
                macl.append(x.macheader)
            break
    remaining-=1
updColletion() #gather new pickle files into the collection file
print macl
open(maclist,'w').write( '\n'.join(macl)+'\n') #dump MAC headers into file

sys.exit(0)
