#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq

ih=inpHand('Usage: mkrmsinp.py -io [-fs]',
           ' -i _R_header (five-letter code)',
           ' -f _A_fastaf fasta file',
           ' -s _A_seqf secondary structure propensity file',
           ' -o _AR_outd output directory'
           )
ih.parse(locals(),sys.argv)
print '...running mkrmsinp.py...'

if not fastaf and not seqf:
    #print codedir+'/python/seq/mkfasta.py -i '+header+' -d .';sys.exit(0)
    print '...mkrmsing.py calling mkfasta.py...'
    os.system(codedir+'/python/seq/mkfasta.py -i '+header+' -d .')
    fastaf=os.getcwd()+'/seq.fasta'
n=0

if fastaf:
    handle=open(fastaf,'r')
    handle.readline()
    for line in handle.readlines():
        n=n+len(line)
        if line[-1]=='\n': n=n-1
    handle.close()
else: n=len( open(seqf,'r').readlines() )

if not os.path.exists(outd): os.makedirs(outd);
outf=open(outd+'/rmsinp','w',)
outf.write('1 '+`n`+'\n'+`n`+'\n'+header+'\n')
outf.close()

print '...finished mkrmsinp.py...'
sys.exit(0) 
