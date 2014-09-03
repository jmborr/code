#!/usr/bin/python
'''
Select those alignments with sequence identity in between 0.5 and 0.99 after global pairwise
alignment
'''
import os,sys,re
from utilities.small_utilities import chomp,junkName
from seq.alignYangManageResults import alignOut

#global variables
rootc='/gpfs1/active/jose/code/projects/spliceVariants/benchmark/preparing_database'
root='/gpfs1/scratch/jose/spliceVariants/benchmark/preparing_database'
ph=re.compile('>\s*(.*)\s*\n') #extract the header out of the fasta header line
paa=re.compile('[G,A,V,L,I,S,T,C,M,P,D,N,E,Q,K,R,H,F,Y,W,B,Z,X]') #check for amino acid letter
listAlgs=''
minSeqIdent=0.5
minmaxgap=5 #minimum gap

#return sequence and header of current fasta entry where pin is pointing to
def importFastaEntry(pin):
    seq='' ;  header=''

    #check end-of-file nad beginning of fasta entry
    while not header:
        header=pin.readline()
        if not header: return '','' #check end-of-file
        header=ph.match(header).group(1) #check beginning of fasta entry

    #read fasta sequence
    line=chomp( pin.readline() )
    while paa.match(line): #correspond to an amino acid letter
        seq+=line
        n=pin.tell()
        line=chomp( pin.readline() )

    if line: pin.seek(n) #go back one line, since we over-read in the previous loop
    return header,seq

#write fasta file to current working directory
def writeFasta(header,seq):
    tmp=seq
    pout=open(header,'w')
    pout.write('>'+header+'\n')
    while tmp:
        pout.write(tmp[0:60]+'\n')
        tmp=tmp[60:]        
    pout.close()

#run fasta on the header and sequence, generating several files
def runfasta(h,s):
    writeFasta(h,s) #dump sequence to fasta file
    cmd=fasta+' '+h+' '+db+' -b 500 -E 0.001 -Q -m 10 -p -w 50 -O '+h+'.out >/dev/null'
    os.system(cmd)    #run fasta
    cmd=fastaout+' '+h+' '+h+'.out 0.99 0.35 0.35 0.99'
    os.system(cmd)    #parse fasta output
    return h+'.log'

#filter hits by sequence identity (remember fasta does a local aligment)
def filterHits(header,log):
    headers2=[]
    pin=open(log,'r')
    for i in range(0,4): pin.readline() #first 4 lines uninteresting
    for line in pin.readlines():
        list=line.split()
        header2=list[1]
        seqIdent=float(list[2])
        if seqIdent>minSeqIdent and seqIdent<0.999: #filter by sequence identity
            headers2.append(header2)
    pin.close()
    return headers2

#impose several filters to accept the aligment as useful
def passFilter(h1,h2,alg):
    seqident=alg.seqIdentity() ; #print ' seqident',seqident
    if seqident>minSeqIdent and seqident<0.999:
        if listAlgs.find(h2+' '+h1) < 0:
            return 1
    return 0

def someInfo(alg):
    out={}
    l1=alg.gs[1].sL                   # L1
    l2=alg.gs[2].sL                   # L2
    l=alg.l                           # length_alignment
    naligned=alg.alignedL()           # number-aligned-residues
    nidentical=alg.identicalL()       # number-identical-residues
    coverage=(1.*naligned)/l           # coverage
    seqident=(1.*nidentical)/naligned # sequence-identity((7)/(6))
    ngap=alg.Ngap()                   # N-terminal gap
    nngap=len(ngap)                   # length of N-terminal gap
    cgap=alg.Cgap()                   # C-terminal gap
    ncgap=len(cgap)                   # length of C-terminal gap
    next=0                            # number of external gaps segments [0,2]
    if nngap>0:next+=1
    if ncgap>0:next+=1
    intgaps=alg.extractInternalGaps()
    nintgaps=len(intgaps)             # number of internal gaps
    intbiggap=0                       # size of biggest interior gap
    for gap in intgaps:
        if len(gap) > intbiggap: intbiggap=len(gap)
    
    return '%3d %3d %3d %3d %3d %5.3lf %5.3lf %3d %3d %3d %3d %3d'%(l1,l2,l,naligned,nidentical,coverage,seqident,nngap,ncgap,next,nintgaps,intbiggap)

#********** EXECUTION BEGINS HERE *********

#load fasta files into memory
seqs={}
pin=open(root+'/pdb.40.300','r')
while 1:
    header,seq=importFastaEntry(pin) #; print header,seq
    if not header: break #reached end-of-file
    seqs[header]=seq
pin.close()

#run fasta for each sequence to output the top 500 hits. Then, for
#every hit do global alignment
pout=open('pdb.40.300.alns','w')      #store alignments
pout2=open('pdb.40.300.alns.idx','w') #store binary indexes
pout3=open('pdb.40.300.alns.dat','w') #store some info about the alignments
tmpd=junkName() ; os.system('/bin/mkdir '+tmpd) ; os.chdir(tmpd) #work in temporary directory
pout3.write('(1)header1 (2)header2 (3)L1 (4)L2 (5)length_alignment (6)number-aligned-residues (7)number-identical-residues (8)coverage (9)sequence-identity (10)length-Ngap (11)length-Cgap (12)number-external-gaps (13)number-internal-gaps (14)length-biggest-internal-gap\n')
fasta=rootc+'/fasta34/fasta34'
fastaout=rootc+'/fastaout'
db=root+'/pdb.40.300'
headers=seqs.keys() ; N=len(headers) #; print 'N=',N
n=0
for header in headers:
    seq=seqs[header]  ;  n+=1 ;  print header,n
    log=runfasta(header,seq) #log file containing hits
    headers2=filterHits(header,log) #; print headers2 #list of filtered headers
    for header2 in headers2:
        #print header,header2 #; print headers2
        seq2=seqs[header2]
        cmd='$HOME/bin/align '+seq+' '+seq2+' 3'
        list=os.popen(cmd).readlines() #; print type(list) ; print list
        alg=alignOut( list ).alg #alignment object from align output
        if passFilter(header,header2,alg)>0: #several filters to consider alignment
            #print alg.info()
            listAlgs+=header+' '+header2+'\n'
            index=pout.tell()            
            pout.write( alg.info()+'\n' ) #write alignment
            #print '%7d %s %s'%(index,header,header2)
            pout2.write( '%9d %s %s\n'%(index,header,header2) ) #write index location
            #print header+' '+header2+' '+someInfo(alg)
            pout3.write(header+' '+header2+' '+someInfo(alg)+'\n')

pout.close()
pout2.close()
pout3.close()
os.system('/bin/rm -rf '+tmpd) #a little bit of cleaning
sys.exit(0)
