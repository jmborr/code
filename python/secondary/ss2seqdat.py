#!/usr/bin/python
"""

Translate secondary structure to a seq.dat file

"""
import os,sys
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand
from seq.fastaManager import importFastaEntry
from seq.letters import one2three

lett2numb0='C 1 H 2 E 4'

def unfold(lett2numb):
    l2ndict={}
    items=lett2numb.split()
    for i in range(0,len(items),2): l2ndict[ items[i] ]= int(items[i+1])
    return l2ndict

def genSeqdat(sequence,secondary, lett2numb=lett2numb0):
    '''
    given a one-letter sequence of amino acids and a one-letter
    secondary structure assignment, produce a seq.dat file. We need
    some translator between the secondary structure assignments and
    the numbers 1,2,4 denoting secondary structure in the seq.dat
    files
    '''
    buf=''
    seqs={'aa':'','ss':''}
    correspondende={'aa':sequence,'ss':secondary}
    #we treat in the same way the amino acid sequence and the secondary sequence
    for key in seqs.keys():
        value=correspondende[key]
        if isinstance(value,str): #string, either a file or a one-letter sequence
            s1=''
            if os.path.exists(value): #it's a file
                line=open(value,'r').readline()
                if line[0]=='>': #it's a fasta format
                    header,s1=importFastaEntry( open(value,'r') )
                else: s1=line.strip() #it's a single-line file
            else: s1=value.strip()   #it's a string, assume one-letter code of amino acids
            seqs[key]=s1
    if len(seqs['aa'])!=len(seqs['ss']):
        sys.stderr.write('Error from ss2seqdat.genSeqdat: different lengths\n')
        return ''
    l2ndict=unfold(lett2numb) #secondary letters to seq.dat 1,2,4 types  #print dict;sys.exit(1)
    for i in range( len(seqs['aa']) ):
        a=seqs['aa'][i] ; b=seqs['ss'][i] #;print a,one2three[a],b,l2ndict[b]#;sys.exit(1)
        buf+='%4d  %s  %d  10\n'%(i+1,one2three[a],l2ndict[b])
    return buf

###################################################################################

if __name__=='__main__':
    inpHand('Usage: ss2seqdat.py [options]',
            ' -a _RA_ssf file with secondary structure, either fasta or single-line file',
            ' -b _RA_seqf file with amino acid sequence, either fasta or single-line file',
            ' -c _A_outf output file (def:./seq.dat)',
            ' -d __l2n letter to number (def:"C 1 H 2 E 4")'
            ).parse(locals(),sys.argv)

    if not outf:outf='seq.dat'
    if not l2n:l2n=lett2numb0
    open(outf,'w').write( genSeqdat(seqf,ssf,lett2numb=l2n) )

    sys.exit(0)
