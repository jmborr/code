#!/usr/bin/python

import sys,re,os
from seq.msa import gen_msa,output_msa
from seq.get_seq import get_seq
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.small_utilities import junkName,Bye
from utilities.codedir import codedir
from jobs.job import pastry
from random import randint

root=codedir+'/perl/tasser/tasser1.0/mkinput'

def gen_pair(msa,seqdatf,outd,co=0.3,libd='/library/yzhang/PDB/'):
    seqdatf=addAbsPath(seqdatf)  ;  outd=addAbsPath(outd)
    if libd[-1]!='/':libd=libd+'/' #dummy neccessity
    #create temporary workding dir
    tmpd=junkName()  ;  pastry('/bin/mkdir -p '+tmpd)  ;  os.chdir(tmpd)
    
    output_msa(msa,'msa.aln')
    pastry('/bin/cp '+libd+'/list .')
    pastry('/bin/cp '+root+'/bin/pair65 pair')
    pastry('/bin/cp '+root+'/data/*comm .')
    os.system('/bin/cp '+seqdatf+' seq.dat')
    pastry('./pair '+`co`+' '+libd)

    #remove temporary directory
    if not os.path.exists(outd): pastry('/bin/mkdir -p '+outd)
    pastry('/bin/mv pair.1 '+outd+'/pair1.dat')
    pastry('/bin/mv pair.3 '+outd+'/pair3.dat')
    os.chdir('..')  ;  pastry('/bin/rm -rf '+tmpd)
    
#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: mkpair.py [options]\nreturns pair1.dat and pair3.dat files',
               '  -a _A_fastaf sequence file in fasta format',
               '  -b _A_seqf raw sequence file',
               '  -c _A_pdbf pdb file',
               '  -d __header pdb header (one of either -a -b -c -d)',
               '  -e _A_seqdatf seq.dat-like  file containing formatted secondary structure',
               '  -f __co sequence cutoff (default=0.3)',
               '  -g _A_libd library of pdb,mtx,cnt2 files (def:/library/yzhang/PDB)',
               '  -o _A_outd output directory(default:current dir)')
    ih.parse(locals(),sys.argv)

    if not outd: outd=os.getcwd()     #default is current directory
    #produce the multiple sequence alignment, by passing the target sequence and a header
    #msa={'header':???,'target':???,'alignments':[???,???,...,???]}
    if fastaf:
        sequence=get_seq(fastafile=fastaf)
        if not header:
            header=re.compile('>\s*(\S*)\s*').match(open(fastaf,'r').readline()).group(1)
    elif seqf:
        sequence=get_seq(seqfile=seqf)
    elif pdbf:
        sequence=get_seq(pdbfile=pdbf)
        header=pdbf
    elif header:
        sequence=get_seq(pdbheader=pdbh)
    else:
        sys.stderr.write('ERROR: supply either one of -a -b -c -d')
        sys.exit(1)        
    if not header: header='protein'
    if not co: co=0.3
    else: co=float(co)
    
    #generate or find secondary structure prediction if not supplied
    if not seqdatf:
        if pdbh: pastry('mkseq.py -e '+pdbh) #maybe seq.dat already done
        else: #generate a seq.dat file
            seqdatf=os.getcwd()+'/seq.dat'
            from secondary.mkseq import output_sec
            output_sec(sequence,'seq.dat')

    if not libd: libd='/library/yzhang/PDB/'
    
    msa=gen_msa(sequence,header=header)
    print 'finished generating alignment'
    gen_pair(msa,seqdatf,outd,co,libd=libd)

    sys.exit(0)
