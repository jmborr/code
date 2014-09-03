#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand
from utilities.codedir  import codedir,pdbnonhom
from random import randint
from jobs.job import pastry

ih=inpHand('Usage: mkseq.py -w -i -s -r -s',
           '  -w _A_workd temporaty working directory (def=./tmp_mkseq$RANDOM)',
           '  -p _A_packd directory where mkseq.pl and related are located (def=codedir/perl/tasser/tasser1.0/mkinput)',
           '  -o _RA_outd output directory',
           '  -i __pdbh file letter pdb-header',
           '  -s _A_rawf raw sequence file',
           '  -f _A_fastaf sequence file in fasta format',
           '  -l __log  log file')
ih.parse(locals(),sys.argv)
if not workd: workd=os.getcwd()+'/tmp_mkseq'+`randint(0,1000000)`
if not packd: packd=codedir+'/perl/tasser/tasser1.0/mkinput'

print '...running mkseq.py...'

#print 'fastaf='+fastaf
if (not pdbh) and (not rawf) and (not fastaf):
    ih.abort('  !!Needed at least one of flags i,s !!')
cmd='/bin/mkdir -p '+workd; pastry(cmd)
cmd='/bin/mkdir -p '+outd;  pastry(cmd);
currd=os.getcwd() #print workd,outd,currd
#obtain raw sequence file from different sources
if pdbh:
    #see if secondary structure file is already made
    cmd='find '+pdbnonhom+' -name "'+pdbh+'.SEQ"' # ;print cmd;sys.exit(1)
    seqf=os.popen(cmd).readline()[:-1]
    if seqf:
        #print 'copy the found file!';sys.exit(1)
        cmd='cp '+seqf+' '+outd+'/seq.dat';
        pastry(cmd)
        pastry('/bin/rm -rf '+workd)
        sys.exit(0)  #will this stop also the parent processes?
    if (not rawf) and (not fastaf):
        cmd=codedir+'/python/seq/mkraw.py -i '+pdbh+' -o seq.raw -d '+workd
if fastaf:
    cmd=codedir+'/python/seq/mkraw.py -f '+fastaf+' -o seq.raw -d '+workd
if rawf:
    cmd='/bin/cp '+rawf+' '+workd+'/seq.raw';
pastry(cmd)
#move all executables and related to the working directory
os.chdir(packd)          
cmd='/bin/cp -rf mkseq.pl bin/ data/ '+workd;
print '...mkseq.py calling mkseq.pl...'
pastry(cmd)
os.chdir(workd)
#create the secondary structure file
cmd='./mkseq.pl  -r ./ -s seq.raw'  ;pastry(cmd);

#move output and remove temporary working directory, unless workd and currd
#are the same
cmd='/bin/mv seq.dat '+outd      ;pastry(cmd)
os.chdir(currd)
if workd != currd: pastry('/bin/rm -rf '+workd)
print '...finished mkseq.py...'
sys.exit(0)
