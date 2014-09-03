#!/usr/bin/python

# This is a simple script which will carry out all of the basic steps
# required to make a PSIPRED V2 prediction. Note that it assumes that the
# following programs are in the appropriate directories:
# blastpgp - PSIBLAST executable (from NCBI toolkit)
# makemat - IMPALA utility (from NCBI toolkit)
# psipred - PSIPRED V2 program
# psipass2 - PSIPRED V2 program

import sys,os,re
from random import randint
from inputArgs.inputArgs import inpHand
from utilities.codedir  import nrfiltdb,ncbidir,psipredd
from utilities.small_utilities import Bye,junkName
from jobs.job import pastry
from seq.blastManager import blastRun
from seq.fastaManager import importFastaEntry

#run psipred on a fasta formatted file, and output the prediction to a list
def run_psipred(fastaf,mtxf=''):
    '''mtxf: file with profile done by makemat
    '''
    header,seq=importFastaEntry( open(fastaf,'r') )
    #create temporary directory and copy fasta file there
    tmpdir=junkName()
    pastry('mkdir -p '+tmpdir)
    pastry('/bin/cp -f '+fastaf+' '+tmpdir+'/'+header+'.fasta')
    if mtxf: pastry('/bin/cp '+mtxf+' '+tmpdir+'/'+header+'.mtx')
        
    os.chdir(tmpdir)

    if not mtxf:
        blastOuts=blastRun(header+'.fasta')
        pastry('echo '+blastOuts['chk']  +' > '+header+'.pn')
        pastry('echo '+blastOuts['fasta']+' > '+header+'.sn')
        pastry(ncbidir+'/makemat -P '+header)
        
    mtxf=header+'.mtx'

    #run first pass
    print 'Pass1 ...'
    pastry('/bin/cp '+psipredd+'/data/data.tbz2 .')    #copy data dir to tmpdir
    pastry('tar jxf data.tbz2')
    pastry(psipredd+'/bin/psipred '+mtxf+' weights.dat weights.dat2 weights.dat3 weights.dat4 > psitmp.ss')

    #run second pass
    print 'Pass2 ...'
    cmd=psipredd+'/bin/psipass2 weights_p2.dat 1 1.0 1.3 psitmp.ss2 psitmp.ss'
    prediction=os.popen(cmd).readlines()

    #remove temporary directory
    os.chdir('..')
    pastry('/bin/rm -rf '+tmpdir)
    
    return prediction

#run prediction with a raw sequence, and maybe a header
def run_psipred2(sequence,header='>xxxxx',n=60):
    #create temporary fasta file
    fastaf='protein'+`randint(0,1000000)`+'.fasta'
    pastry('echo "'+header+'" > '+fastaf)
    while len(sequence):
        pastry('echo '+sequence[0:n]+' >> '+fastaf)
        sequence=sequence[n:]
    prediction= run_psipred(fastaf)
    pastry('/bin/rm '+fastaf)
    return prediction

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: runpsipred.py [options]',
               '  -a _A_fastaf sequence file in fasta format',
               '  -b __rawseq one-single string sequence (needed -a or -b)',
               '  -c __outf output file name(default:protein.horiz)',
               '  -o _A_outd output directory(default:current dir)')
    ih.parse(locals(),sys.argv)

    #compose absolute path to output file
    if not outf: outf='protein.horiz'
    else: outf=os.path.basename(outf) #remove absolute path
    if not outd: outd=os.getcwd()     #default is current directory
    outf=outd+'/'+outf        #compose absolute output filename

    if not fastaf:
        if not rawseq:
            sys.stderr.write('ERROR: either one of -a or -b options')
            sys.exit(1)
        prediction=run_psipred2(rawseq)
    else:  prediction=run_psipred(fastaf)
    open(outf,'w').writelines(prediction)
    os.system('sync')
    sys.exit(0)
