#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.codedir import codedir
from utilities.small_utilities import junkName,Bye,chomp
from jobs.job import pastry
from seq.fastaManager import importFastaEntry,writeFasta
from seq.blastManager import blastRun,blastAlignment

wgtf=codedir+'/perl/tasser/tasser1.0/mkinput/data/wgt.tar.bz2'
#solve=codedir+'/perl/tasser/tasser1.0/mkinput/bin/solve' #deprecated
solve=codedir+'/c/pdb_utils/solve'

ww=['0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70',\
    '0.75','0.80','0.85']

def gen_solve(fastaf,outd='.'):
    outd=addAbsPath(outd)
    header,seq=importFastaEntry(open(fastaf,'r'))
    tmpd=junkName()  ;  pastry('/bin/mkdir '+tmpd)  ;  os.chdir(tmpd) ; pastry('/bin/cp '+solve+' .')
    writeFasta('protein',seq,filename='protein.seq')     #appropriate fasta file
    pastry('/bin/cp '+wgtf+' . && tar jxf wgt.tar.bz2') #bring weight files
    
    blastOuts=blastRun('protein.seq',opt='-Q protein.mat3 -o protein.out -b 999 -e 1e-3 -m 6') #blast run with appropriate options

    res={}
    res['Lch']=re.compile('(\d+)')
    res['exp']=re.compile('(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)')

    exp={}
    for w in ww:
        cmd='./solve weight.'+w+' protein'  ;  os.system(cmd)
        neu=open('protein.neu','r') ; Lch=int( chomp(neu.readline()) )
        exp[w]=[]
        for i in range(0,Lch): exp[w].append( int(chomp(neu.readline())[-1]) )
        neu.close()

    expp=open('exp.dat','w') #note we do one directory up
    expp.write(`Lch`+'  ')
    for w in ww: expp.write(' %5s'%(w))
    expp.write('\n')
    for i in range(0,Lch):
        expp.write('%5d'%(i+1))
        for w in ww: expp.write(' %5d'%(exp[w][i]))
        expp.write('\n')
    expp.close()

    pastry('/bin/mv -f exp.dat '+outd)
    os.chdir('..')  ;  pastry('/bin/rm -rf '+tmpd)  #clean-up


    
if __name__=='__main__':
    
    #parse input options
    inpHand('Usage: mksolve.py [options]',
            '  -a _RA_fastaf sequence file in fasta format',
            '  -b _A_outd output directory (def: current dir)'
            ).parse(locals(),sys.argv)

    if not outd: outd=os.getcwd()
    
    gen_solve(fastaf,outd=outd)

