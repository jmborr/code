#!/usr/bin/python

'''
Generate the fasta files for each of the homologous sequences
'''
import os,sys
from utilities.small_utilities import junkName
from seq.fastaManager import importFastaEntry,writeFasta
from seq.alignYangManageResults import alignOut

#batch_entrez contains all hits by psiblast to g8000
pin=open('/gpfs1/scratch/jose/G8/input/batch_entrez.dat','r')
pout=open('/gpfs1/scratch/jose/G8/input/genbank2header.dat','w')

g8000='/gpfs1/scratch/jose/G8/input/g8000'
junk=junkName()
n=1
header,seq=importFastaEntry(pin) #first sequence actually contains g8000, thus avoid
while header:
    header,seq=importFastaEntry(pin) #; print header 
    writeFasta(junk,seq) #write fasta entry to temporary file
    alignResults=os.popen('$HOME/bin/align '+g8000+' '+junk+' 0').readlines()#global alignment
    alignment=alignOut(alignResults).alg #alignment object
    si=alignment.seqIdentity()
    begin,end=alignment.boundingIndexes() #bounding indexes of first sequence
    seq=alignment.extract(begin,end).gs[2].s #projection of junk onto g8000, remove gaps
    h='g8%03d'%(n) #create a header name
    writeFasta(h,seq) #write projection to fasta file
    print h+' %4.2lf '%(si)+header
    pout.write(h+' %4.2lf '%(si)+header+'\n') #relate header and genbank accesion number
    n+=1
os.system('/bin/rm '+junk) #remove temporary file
pout.close()
pin.close()
sys.exit(0)
