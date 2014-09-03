#!/usr/bin/python

import sys,re,os
from seq.letters import three2one,is_prot_seq
from utilities.codedir import pdbdir
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,abort
from seq.fastaManager import importFastaEntry

def parse_the_pdbfile(fpt,chainID=''):
    sequence=''
    if chainID=='_': chainID=''
    ding=re.compile('ATOM\s+\d+\s+CA\s+(\S+)\s+'+chainID)
    line=fpt.readline()    #;print 'line='+line
    while line and not ding.match(line):
        line=fpt.readline()        
    end=re.compile('[END|TER]')
    while line and not end.match(line):
        p=ding.match(line)
        if p: sequence=sequence+three2one[p.group(1)]
        line=fpt.readline()
    return sequence

#return a string containing whole sequence. Give one of the four different inputs
def get_seq(geninput,inputType=''):
    sequence=''
    if not inputType: #generic input.
        if re.compile('\.fasta$').search(geninput):
            sequence=get_seq(geninput,inputType='fasta')
        elif re.compile('\.seq$').search(geninput):
            sequence=get_seq(geninput,inputType='seq')
        elif re.compile('\.pdb$').search(geninput):
            sequence=get_seq(geninput,inputType='pdb')
        else: sequence=get_seq(geninput,inputType='header')

    if inputType=='fasta':        
        header,sequence=importFastaEntry( open(geninput,'r') )
    elif inputType=='seq':
        sequence= ''.join( chomp(open(geninput,'r').readlines()) )
    #for the pdb file, will read only the first chain, or the first model
    elif inputType=='pdb':
        sequence=parse_the_pdbfile( open(geninput,'r') )
    #assume we'll find corresponding pdbfile
    elif inputType=='header':
        putative=pdbdir+'/pdb'+geninput[0:4]+'.ent'
        if not os.path.exists(putative):abort('ERROR: no '+putative+' pdb file found in '+pdbdir)
        sequence=parse_the_pdbfile(open(putative,'r'),chainID=geninput[4])

    if not sequence:abort('ERROR: get_seq could not produce a single-line sequence\n')
    return sequence

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: get_seq.py [options] returns a file with a single-line sequence',
               '  -a _A_fastaf sequence file in fasta format',
               '  -b _A_seqf raw sequence file',
               '  -c _A_pdbf pdb file',
               '  -d __pdbh pdb header (one of either -a -b -c -d)',
               '  -e _A_geninput (generalized input, either one of the previous)',
               '  -n __outf output file name(default:seq.raw)',
               '  -o _A_outd output directory(default:current dir)')
    ih.parse(locals(),sys.argv)

    #compose absolute path to output file
    if not outf: outf='seq.raw'
    else: outf=os.path.basename(outf) #remove absolute path
    if not outd: outd=os.getcwd()     #default is current directory
    outf=outd+'/'+outf                #compose absolute output filename

    #check we passed sufficient input
    if not fastaf and not seqf and not pdbf and not pdbh and not geninput:
        sys.stderr.write('ERROR: get_seq.py needs either one of -a -b -c -d options')
        sys.exit(1)

    #produce the one-line single-sequence
    sequence=''  ;  inputType=''
    if fastaf:
        geninput=fastaf ;  inputType='fasta'
    elif seqf:
        geninput=seqf   ;  inputType='seq'
    elif pdbf:
        geninput=pdbf   ;  inputType='pdb'
    elif pdbh:
        geninput=pdbh   ;  inputType='header'
    sequence=get_seq(geninput,inputType=inputType)
    open(outf,'w').write(sequence)
    
    sys.exit(0)
