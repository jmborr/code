"""output the size of every internal gap, tail gap, and identical chunk
   found in the gaps_and_ident.dat file
"""
gid='/gpfs1/scratch/jose/spliceVariants/refsAboveC1/gaps_and_ident.dat'
#gid='/gpfs1/scratch/jose/spliceVariants/refsAboveC1/junk.dat'

import os,sys,re
from utilities.small_utilities import chomp

ig='' #list of internal gap sizes
tg='' #list of tail gap sizes
id='' #list of identical chunk sizes

re1=re.compile(':(\s+)(?=:)') #find all internal gaps
re2=re.compile('^(\s+):')     #find N-terminal gap
re3=re.compile(':(\s+)$')     #find C-terminal gap
re4=re.compile('(:+)')        #find all identical chunks
for line in chomp(open(gid,'r').readlines()):
    list=re1.findall(line)
    for match in list: ig+=`len(match)`+'\n'
    match=re2.search(line)
    if match: tg+=`len(match.group())`+'\n'
    match=re3.search(line)
    if match: tg+=`len(match.group())`+'\n'
    list=re4.findall(line)
    for match in list: id+=`len(match)`+'\n'
print '#internal gaps sizes\n'+ig+'#tail gap sizes\n'+tg+'#identical chunk sizes\n'+id
