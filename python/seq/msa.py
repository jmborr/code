#!/usr/bin/python

import sys,re,os
from seq.letters import three2one
from seq.get_seq import get_seq
from utilities.codedir import pdbdir,nrfiltdb,ncbidir
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,mysplit,abort,Bye
from jobs.job import pastry
from random import randint
from seq.blastManager import blastRun,parse_blast
from seq.fastaManager import importFastaEntry 

#returns a dictionary with a multiple sequence alignment
#msa={'header':???,     #some header for this alignment
#     'target':???,     #string containing target sequence
#     'alignments':[???,???,...,???]  #list of strings, each one alignment
#    }
def gen_msa(sequence,header='protein',Eco=1.0,maxId=0.9,minId=0.35,red=0.9,blastOuts={}):
    '''Eco: E=value cut-off
     maxId: max seq ident to sequence
     minId: min seq ident to sequence
     red: no two seq can have more than red redundancy after projected onto sequence
     '''
    msa={}
    msa['header']=header
    msa['target']=sequence
    if not blastOuts:
        blastOuts=blastRun([header,sequence]) #run blast
    parsed_algs=parse_blast(blastOuts['blast'],Eco,maxId,minId,red,querySeq=sequence)
    if not parsed_algs:
        sys.stderr.write('ERROR seq.msa.gen_msa parse_blast returned empty alignments\n')
    msa['alignments']=parsed_algs
    return msa

#output file format as that of mkpair.pl
def output_msa(msa,outf):
    algnList=msa['alignments']  #list with the alignments
    n_msa=len(algnList)         #number of sequences aligned
    lch=len(msa['target'])      #length of the target sequence
    h=msa['header']             #a header for the alignment
    fpt=open(outf,'w')
    fpt.write('%5d %5d %s\n'%(n_msa+1,lch,h))
    #break sequence into 50 character-chunks. Append '*\n' at end of each chunk, except last chunk
    lines=mysplit(msa['target'],append='*\n',chunksize=50)
    lines[-1]=lines[-1].replace('*','')
    fpt.writelines(lines)
    for algn in algnList: fpt.writelines(mysplit(algn,append='\n',chunksize=50))
    return n_msa

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: msa.py [options]',
               '  -a _A_fastaf sequence file in fasta format',
               '  -b _A_seqf raw sequence file',
               '  -c _A_pdbf pdb file',
               '  -d __pdbh pdb header (one of either -a -b -c -d)',
               '  -n __outf output file name(default:msa.aln)',
               '  -o _A_outd output directory(default:current dir)',
               '  -h __help returns a file with a multiple sequence alignment againt PSI-BLAST search')
    ih.parse(locals(),sys.argv)

    #compose absolute path to output file
    if not outf: outf='msa.aln'
    else: outf=os.path.basename(outf) #remove absolute path
    if not outd: outd=os.getcwd()     #default is current directory
    outf=outd+'/'+outf                #compose absolute output filename

    #check we passed sufficient input
    if not fastaf and not seqf and not pdbf and not pdbh:
        abort('ERROR: get_seq.py needs either one of -a -b -c -d options')

    #produce the multiple sequence alignment, by passing the target sequence and a header
    #msa={'header':???,'target':???,'alignments':[???,???,...,???]}
    if fastaf: header,sequence=importFastaEntry( open(fastaf,'r') )
    elif seqf: sequence=get_seq(seqf,inputType='seq')
    elif pdbf:
        sequence=get_seq(pdbf,inputType='pdb')  ;  header=pdbf
    elif pdbh:
        sequence=get_seq(pdbh,inputType='header')  ;  header=pdbh

    if not header: header='protein'
    msa=gen_msa(sequence,header=header)

    #output alignment to file
    output_msa(msa,outf)
    
    sys.exit(0)
