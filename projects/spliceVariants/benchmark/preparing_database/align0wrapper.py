#!/usr/bin/python

# wrapper for fasta "align" program

import os,sys,string

prog='~/bin/align0'
#b62='/net/dell/02/users/adrian/bin/fasta34/blosum62.mat'

def align(seq1,seq2,gap_open=-10,gap_ext=-2):
    fw=open('seq1.tmp','w')
    fw.write('>seq1\n')
    for i in range(len(seq1)/60+1):
	fw.write(seq1[i*60:i*60+60]+'\n')
    fw.close()
    fw=open('seq2.tmp','w')
    fw.write('>seq2\n')
    for i in range(len(seq2)/60+1):
	fw.write(seq2[i*60:i*60+60]+'\n')
    fw.close()
    f=os.popen('%s -f %d -g %d seq1.tmp seq2.tmp 2> /dev/null' % 
    (prog,gap_open,gap_ext))
    
    while 1:
	line=f.readline()
	if string.find(line,'Global alignment score')>0: break
	
    als1=''
    als2=''
    while 1:
	line=f.readline() # empty line before next block
	if not line: break
	line=f.readline() # residue numbers of seq1
	line=f.readline() # seq1 aligment
	als1=als1+string.strip(line[7:67])
	line=f.readline() # line showing matches
	line=f.readline() # seq2 alignment
	als2=als2+string.strip(line[7:67])
	line=f.readline() # residue numbers of seq2
    f.close()
    
    map12={}
    map21={}
    rn1=0
    rn2=0
    for i in range(len(als1)):
	if als1[i]!='-':
	    rn1=rn1+1
	if als2[i]!='-':
	    rn2=rn2+1
	if als1[i]!='-':
	    if als2[i]!='-':
		map12[rn1]=rn2
	    else:
		map12[rn1]=-1
	if als2[i]!='-':
	    if als1[i]!='-':
		map21[rn2]=rn1
	    else:
		map21[rn2]=-1
    
    return (map12,map21)
    
if __name__=='__main__':
    f=open(sys.argv[1])
    f.readline()
    seq1=''
    while 1:
	line=f.readline()
	if not line:
	    break
	seq1=seq1+string.strip(line)
    f.close()
    f=open(sys.argv[2])
    f.readline()
    seq2=''
    while 1:
	line=f.readline()
	if not line:
	    break
	seq2=seq2+string.strip(line)
    f.close()
    (map12,map21)=align(seq1,seq2)
    leni=len(seq1)
    lenj=len(seq2)
    iden=0
    nali=0
    for ri in map12.keys():
	rj=map12[ri]
	if rj==-1: continue
	nali=nali+1
	if seq1[ri-1]==seq2[rj-1]:
	    iden=iden+1

    print '%5d %5d %5d %5d' % (iden,nali,leni,lenj)
    
    
    
