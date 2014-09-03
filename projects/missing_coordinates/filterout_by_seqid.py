#!/usr/bin/python

import os,sys,re
from utilities.small_utilities import chomp,junkName,Bye
from inputArgs.inputArgs import inpHand
from utilities.codedir import scratchdir
from seq.align import do_alignment,seqId
from seq.fastaManager import importFastaEntry

#global variables
pattern=re.compile('#(\d+-\d+)')
alignx=os.path.join('/tmp/jose',junkName())
os.system('/bin/cp ~/bin/align '+alignx) #bring align to local node

def filter_homolog_chains(homologs,maxseqid):
    """filter out homolog sequences"""
    nseqs=len(homologs) #;Bye(nseqs)
    headers=homologs.keys() #;Bye(headers)
    if nseqs>1:
        i=0
        while i<len(headers)-1:
            header_i=headers[i]
            seq_i=homologs[header_i]
            j=i+1
            while j<len(headers):
                header_j=headers[j]
                seq_j=homologs[header_j]
                seqid=float(os.popen(alignx+' '+seq_i+' '+seq_j+' 3|grep "Sequence identity"|cut -d ":" -f 2|cut -d "(" -f 1').readline())
                if seqid > maxseqid: #we find a homolog
                    headers.remove(header_j)
                else:
                    j+=1
            i+=1
    non_homologs={}
    for header in headers: non_homologs[header]=homologs[header]
    return non_homologs
    
#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: filterout_by_seqid.py [options]',
               ' -a __maxseqid remove homologs below this sequence identity (def: 0.9)',
               ' -b _A_outf output file (def: STDOUT)',
               )
    ih.parse(locals(),sys.argv)

    #init default variables
    if not maxseqid: maxseqid=0.9
    else: maxseqid=float(maxseqid)
    if not outf: outp=sys.stdout
    else: outp=open(outf,'w')

    #set has following structure: set={ pdb1=={header1:seq1,header2:seq2}, pdb2.. }
    #set--pdb1---header1:seq1
    #  |     |___header2:seq2
    #  |__pdb2---header1:seq1
    #        |___header2:seq2
    set={}
    file=os.path.join(scratchdir,'missing_coordinates/preparing/missingXray.080503.txt')
    pin=open(file,'r')
    remaining=int( os.popen('grep ">" '+file+'|wc -l').readline().split()[0] )
    while (True):
        header,seq=importFastaEntry(pin) #read one fasta entry
        remaining-=1
        if header:
            pdb=header[0:4] #pdb id, which is composed of four letters
            if pdb not in set.keys():
                set[pdb]={header:seq}
            else:
                set[pdb][header]=seq
        else:
            break #end-of-file reached
        print remaining

    #within each "pdb", remove those chains that are homologs.
    remaining=len(set)
    for pdb in set.keys():
        set[pdb]=filter_homolog_chains(set[pdb],maxseqid)
        remaining-=1
        print remaining
        
    #print filtered set
    buf=''
    for pdb in set.keys():
        non_homologs=set[pdb]
        for header in non_homologs.keys(): buf+=header.replace(':','')+'\n'
    outp.write(buf)
    
    sys.exit(0)
        
            
                      
    
