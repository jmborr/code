#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,Bye,junkName
from utilities.codedir import scratchdir
from secondary.seqdat2ss import translateSec
from seq.fastaManager import importFastaEntry
from seq.alignYangManageResults import alignOut
from math import sqrt

#initialization of several variables
infinity=99
minChi2=infinity
mindifsum=infinity
average=0
std=0
Zscore=infinity
begin=-1
end=-1
pdbBest='xxxxx'
nhits=0
bestChunks=''
bestSeq=''
percentDev=1.0

ih=inpHand('Usage: mkfasta.py [options]',
           ' -a _AR_seqdatf seq.dat file',
           ' -b _A_sstxt ss.txt file (def:scratchdir/foldsOrigin/ss.txt)'
           )
ih.parse(locals(),sys.argv)

if not sstxt: sstxt=scratchdir+'/foldsOrigin/ss.txt'

patt=re.compile( '[^H,E]') #anything but H and E
patu=re.compile( 'E+|H+')  #only chunks of H's and E's
paa=re.compile('[H,E,\s,S,T,G,B,I,-]')
pr={'H':re.compile('H+'),'E':re.compile('E+'),'C':re.compile('C+')}

def printComment():
    print '# nhits(1) Nchunks(2) minDifSum(3) minChi2(4) Zscore(5) average(6) std(7) pdbBest(8) begin(9) end(10)'
    print '# (1) (2)  (3)   (4)      (5)    (6)    (7)   (8)  (9) (10)'

def smoothSeq(seq):
    #remove singlets
    seq=seq.replace('CHC','CCC')
    seq=seq.replace('CEC','CCC')
    seq=seq.replace('HCH','HHH')    
    seq=seq.replace('ECE','EEE')
    return seq
    
def reduceSeq(seq):
    reduced=seq
    for symbol in pr.keys():
        reduced=pr[symbol].sub(symbol,reduced)
    return pr['C'].sub('',reduced)
    
#check alignSS.x is present in the working directory
if not os.path.exists('./alignSS.x'): Bye('ERROR: copy alignSS.x to working directory, please')

#write target secondary structure to temporary file in fasta-like format
target_seq=translateSec(seqdatf) #;Bye(target_seq)
target_seq=smoothSeq(target_seq)
target_seqR=reduceSeq(target_seq) #;Bye(target_seqR)
target_seqRL=[]  #store the position and extent of all chunks of H's and E's
iterator=patu.finditer(target_seq)
for match in iterator:
    (b,e)=match.span()
    target_seqRL.append( (b,1+e-b) )
#print target_seq
#print target_seqR
#print target_seqRL

Ntarget=len(target_seq) #;Bye(`Ntarget`)
Nchunks=len(target_seqR)

#read all appropriate sequences of ss.txt file to memory, and simplify symbols
sss=open(sstxt,'r') #file containing secondary structure assignments
(header,seq)=importFastaEntry(sss,paa=paa)
while header:
    seq=patt.sub('C',seq) #substitute anything but H and E with C
    seq=smoothSeq(seq)
    seqR=reduceSeq(seq) #reduce sequence to series of helixes and strands
    if seqR.find(target_seqR)>=0: #at least one occurrence of target_seqR in seqR
        #print seq
        #print seqR
        seqRL=[]
        it=patu.finditer(seq) #store the position and extent of all chunks of H's and E's
        for match in it:
            (b,e)=match.span()
            seqRL.append( (b,1+e-b) )
        #find all occurrences of target_seqR in seqR
        d=seqR.find(target_seqR)
        n=d
        while d>=0:
            #print header,'n=',n
            nhits+=1
            Chi2=0.0
            difsum=0
            for i in range(Nchunks):
                #print 'i=',i
                #print 'target_seqRL[i]=',target_seqRL[i]
                #print 'n+i=',n+i,'len(seqR)=',len(seqR),'len(seqRL)=',len(seqRL),'seqRL[n+i]=',seqRL[n+i]
                dif=float( target_seqRL[i][1] - seqRL[n+i][1] )
                Chi2+=dif*dif/target_seqRL[i][1]
                difsum+=abs(dif)
            average+=Chi2
            std+=Chi2*Chi2
            if Chi2<minChi2:
                bestChunks=seqR[n:n+Nchunks]
                pdbBest=header
                minChi2=Chi2
                mindifsum=difsum

                begin=seqRL[n][0] #beginning position of the best match
                end=seqRL[n+Nchunks-1][0]+seqRL[n+Nchunks-1][1] #ending position of the best match
                bestSeq=seq[begin:end+1]
            d=seqR[n+1:].find(target_seqR) #next position for beginning of match
            n=n+1+d
            
            #print '%5d %6.3lf %6.3lf %7.3lf %6.3lf %6.3lf %s %3d %3d'%(nhits,Chi2,minChi2,Zscore,average,std,header,begin,end)

    (header,seq)=importFastaEntry(sss,paa=paa)

if nhits: #assuming we found a match in sstxt
    percentDev=float(mindifsum)/len(target_seq.replace('C',''))
    mindifsum=float(mindifsum)/Nchunks
    minChi2=minChi2/Nchunks #avoid scaling of Chi2 with number of secondary elements
    average=float(average)/(nhits*Nchunks)
    std=sqrt(float(std)/nhits-average*average)/Nchunks
    if std: #if there's only one hit, then std==0
        Zscore=(minChi2-average)/std

printComment()
print '%5d %2d %4.1lf %4.2f %6.3lf %7.3lf %6.3lf %6.3lf %s %3d %3d'%(nhits,Nchunks,mindifsum,percentDev,minChi2,Zscore,average,std,pdbBest,1+begin,1+end)
#print target_seq
#print bestSeq

sys.exit(0)
