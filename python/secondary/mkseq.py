#!/usr/bin/python

import re,sys,os
from secondary.psipred.runpsipred import run_psipred2
from seq.letters import one2three
from seq.get_seq import get_seq
from inputArgs.inputArgs import inpHand
from utilities.codedir import libseq
from utilities.small_utilities import Bye
from jobs.job import pastry

#obtain secondary structure prediction with format suitable for
#TASSER in a list from a single-line containing the sequence
#Threee consecutive lines with info look like this:
#Conf: 985147885798999999999884598999983838999999608998999852378988
#Pred: CCCEEEEEECCHHHHHHHHHHHHHCCCEEEEECCHHHHHHHHHHCCCCEEEEECCCCCCC
#  AA: MTKKILIIEDEKNLARFVSLELQHEGYEVIVEVNGREGLETALEKEFDLILLDLMLPEMD
def make_sec(sequence,prediction=None):
    l=[] ; index=1
    #run psipred and store output in a list
    if not prediction:prediction=run_psipred2(sequence)
    for i in range(len(prediction)):
        line=prediction[i]
        #found line beginning with "Conf" ?
        p=re.compile('Conf:\s+(\d+)').match(line)
        if p:
            conf=p.group(1) #obtain confidence numbers
            i=i+1
            line=prediction[i]
            #obtain secondary-structure single letter code
            pred=re.compile('Pred:\s+(\S+)').match(line).group(1)
            i+=1
            #obtain single-letter sequence
            line=prediction[i]
            aa=re.compile('\s+AA:\s+(\S+)').match(line).group(1)
            for j in range(len(aa)):
                n=1
                m=int(conf[j])
                if m>=1:
                    if pred[j] == 'H':   n=2
                    elif pred[j] == 'E': n=4
                x=one2three[aa[j]]
                l.append('%5d   %3s%5d%5d\n'%(index,x,n,m))
                index+=1
        i+=1
    return l

def output_sec(sequence,outf,prediction=None):
    '''prediction is the output from secondary/psipred/runpsipred.py/run_psipred'''
    open(outf,'w').writelines(make_sec(sequence,prediction=prediction))

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: mkseq.py [options]',
               '  -a _A_fastaf sequence file in fasta format',
               '  -b __rawseq one-single string sequence ',
               '  -c _A_seqfile raw sequence file ',
               '  -d _A_pdbf pdbfile ',
               '  -e __pdbh pdb header(one of -a -b -c -d -e) ',
               '  -f __outf output file name(default:seq.dat)',
               '  -o _A_outd output directory(default:current dir)')
    ih.parse(locals(),sys.argv)    
    #absolute path of output file
    if not outf: outf='seq.dat'
    if not outd: outd=os.getcwd()
    outf=outd+'/'+outf
    #check if we passed some input
    if not fastaf and not rawseq and not seqfile and not pdbf and not pdbh:
        sys.stderr.write('ERROR: required one option of a/b/c/d/e')
    #if passed a pdbheader, then see if already made
    if pdbh:
        putative=libseq+'/'+pdbh+'.SEQ'
        if os.path.exists(putative):
            pastry('/bin/mv '+putative+' '+outf)
            sys.exit(0) #finished!
    #obtain a single-line sequence from the input
    if   fastaf  : sequence=get_seq(fastaf,inputType='fasta')
    elif rawseq  : sequence=rawseq
    elif seqfile : sequence=get_seq(seqfile,inputType='seq')
    elif pdbf    : sequence=get_seq(pdbf,inputType='pdb')
    elif pdbh    : sequence=get_seq(pdbh,inputType='header')
    #run the prediction
    if not sequence:
        sys.stderr('ERROR:could not produce a single-sequence line')
        sys.exit(1)
    else:
        output_sec(sequence,outf)
    sys.exit(0)
