#!/usr/bin/python

import os,sys,re
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,junkName
from utilities.codedir import codedir
from seq.blastManager import blastRun
from seq.fastaManager import importFastaEntry
from seq.msa import gen_msa
from seq.letters import one2three

#create frequency file
def createFreqFile(seq,msa):
    '''
    seq: query sequence

    msa: multiple sequence alignment in a list, in the style of 3590
    and e10. msa does not contain the query sequence, and each
    template sequence is a list item    
    '''
    
    aas='ARNDCQEGHILKMFPSTWYV' #amino acid order to output, assumed by jnet
    freq=[]                   #freq[i][j] 0<i<len(seq), 0<j<19
    N=len(seq)
    for i in range(N): freq.append([0]*20) #initialize
    for i in range(N): #add amino acids of query sequence to the frequency table
        m=aas.find(seq[i]) #return position in aas corresponding to amino acid seq[i]
        if m<0:continue    #seq[i] character not found in aas (like "X" or "-")
        freq[i][m]+=1        
    for seq2 in msa:  #do same for sequences in the alignment
        for i in range(N):
            m=aas.find(seq2[i]) #return position in aas corresponding to amino acid seq[i]
            if m<0:continue    #seq[i] character not found in aas
            freq[i][m]+=1
    blastffreq=junkName()+'.jnetfreq'
    out=open(blastffreq,'w')
    for i in range(N):
        line=''  ; m=0
        for j in range(20):m+=freq[i][j] #total number of counts for amino acid at position "i"
        for j in range(20):
            freq[i][j]=round( (10.0*freq[i][j])/m ) #rounded to nearest integer
            line+='%3d'%(freq[i][j])
        out.write(line+'\n') #write amino acid frequencies for amino acid at position "i"
    out.close()
    return blastffreq
       
#########################################################
#       SCRIPT  EXECUTION STARTS HERE
#########################################################
inpHand('Usage: inprosp.py [options]',
        '  -a _RA_fastaf sequence file in fasta format',
        '  -b _A_blasttarf tarred blast output (contains xxxxx.blast,xxxxx.pssm). If not provided, jnet.py will do a psiblast run',
        '  -c _A_outf output file name (def: ./seq.dat)',
        ).parse(locals(),sys.argv)

currd=os.getcwd() #current directory

if not outf: outf=currd+'/seq.dat' #output file

workd=currd+'/'+junkName() ;os.system('/bin/mkdir -p '+workd) ; os.chdir(workd) #temp directory

header,seq=importFastaEntry(open(fastaf,'r'))
pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)

#Retrieve/create  psiblast outputs
if blasttarf:#we passed a *.tar file containing psiblast report and pssm file
    blastOuts={'outf':'', 'blast':'', 'chk':'', 'fasta':fastaf, 'pssm':''}
    os.system('tar xf '+blasttarf)
    blastOuts['blast']=os.popen('ls -1 *.blast').readline().strip() #get name of blast report
    blastOuts['pssm']=os.popen('ls -1 *pssm.').readline().strip()   #get name of pssm file
else: blastOuts=blastRun(fastaf) #run blast with default options

#create multiple sequence alignment projected to the query sequence (3590 or e10 style)
msa=gen_msa(seq,header,Eco=0.0001,maxId=0.85,minId=0.10,red=0.75,blastOuts=blastOuts)['alignments']

#find frequency table from msa, and output to "freq" file. I did this
#function because PSIBLAST has evolved and the perl script "getfreq"
#provided in the jnet distro does not work. I could use newer script
#"parse_psi -freq" but it requires a psiblast report obtained with
#blastpgp -m 6, instead of blastpgp -m 0. I don't want to run PSIBLAST
#twice and running with -m 6 gives some humongously-sized reports
blastffreq=createFreqFile(seq,msa)

#remove too short sequences, thenk keep only first M sequences
msa2=[]  ;  N=len(seq)  ;  Lhalf=N/2
for seqgapped in msa:
    sequngapped=seqgapped.replace('-','')
    if len(sequngapped) > Lhalf: msa2.append(sequngapped)
M=int(1000*200/N) #maximum number of sequences to use
msa=msa2[0:M] #reclaim memory space by liberating the gapped sequences list

#output alignment as suitable use for clustalw
rootname=junkName()
fastasf=rootname+'.aln'  ;  fpt=open(fastasf,'w')
fpt.write('>00000\n'+seq+'\n')
for i in range(len(msa)): fpt.write('>'+'%05d\n'%(i+1)+msa[i]+'\n')
fpt.close()

#run clustalw
os.system('clustalw -OUTORDER=INPUT -INFILE='+fastasf+' -OUTPUT=GCG >/dev/null')
msf=rootname+'.msf'
if not os.path.exists(msf): Bye('ERROR: not msf file generated in jnet.py')

#run perl scripts to create the various inputs required by jnet
pastry('/bin/cp -r '+codedir+'/bin/jnet/perl .')
pastry('/bin/cp -r '+codedir+'/bin/jnet/bin .')
os.system('./perl/msf2jnet '+msf) ; msffa=msf+'.fa'                  #multiple sequence alignment
os.system('./perl/gethmm '+msf+' >/dev/null') ; msfhm=msf+'.hmmprof' #hidden-Markov model
pssmf=blastOuts['pssm']
os.system('./perl/getpssm '+pssmf+' > '+pssmf+'.jnetpssm') ; pssmf=pssmf+'.jnetpssm'

#run jnet and parse to generate seq.dat
jnetout=junkName()
os.system('./bin/jnet -p '+msffa+' '+msfhm+' '+pssmf+' '+blastffreq+' > '+jnetout)
pattern=re.compile(':\s(\S+)\n')  ;  final=''  ;  conf=''
ss2nn={'-':1, 'H':2, 'E':4}
for line in os.popen('grep -P "\sFINAL\t" '+jnetout).readlines():
    final+=pattern.search(line).group(1)
for line in os.popen('grep -P "\sCONF\t" '+jnetout).readlines():
    conf+=pattern.search(line).group(1)
out=open(outf,'w')
for i in range(N):
    out.write(  '%5d%6s%5d%5d\n'%( i+1, one2three[seq[i]], ss2nn[final[i]], int(conf[i]) )  )    
out.close()

#clean-up working directory
os.chdir(currd)
#os.system('/bin/rm -rf '+workd)
sys.exit(0)
