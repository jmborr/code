#!/usr/bin/python

import os,sys,re
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp
from jobs.job import job

from seq.fastaManager import importFastaEntry

inpHand('Usage: dukka_pulchra.py [options]',
        '  -a _A_fastaf sequence file in fasta format',
        '  -b _RA_spktar output from spicker',
        '  -c _RA_outdir output directory',
        '  -d __header header (def:none)',
        '  -e __infix infix, append after header and and extra dot (def:none)',
        '  -f __tarname give a *.tar.bz2 file name if you desire to tar & compress results and then send them to outdir (def:none)'
        ).parse(locals(),sys.argv)
header2=''
if fastaf:
    pastry('/bin/cp '+fastaf+' .') ;fastaf=os.path.basename(fastaf)
    header2,seq=importFastaEntry(open(fastaf,'r'))
if not header and fastaf: header=header2
prefix=''
if header and infix: prefix=header+'.'+infix+'.'

taropts='xf'
if 'bz2' in spktar: taropts='jxf'
pastry('/bin/tar '+taropts+' '+spktar)

#find list of cluster numbers
clusterlist=chomp(os.popen('ls -1 combo*.pdb').readlines())
if clusterlist==None:sys.exit(1)
recombo=re.compile('combo(\d+)\.pdb')
for i in range(len(clusterlist)):
    clusterlist[i]=int( recombo.search(clusterlist[i]).group(1) )

#now refine on every combo
pastry('/bin/cp /gpfs1/u/jose/bin/dukka.x /gpfs1/u/jose/bin/pulchra .')
for i in clusterlist:
    closcf='closc'+`i`+'.pdb'
    pastry('pulchra -b -s '+closcf) #only keep the CA's
    combof='combo'+`i`+'.pdb'
    pastry('./dukka.x '+closcf+'.rebuilt '+combof)
    pastry('./pulchra best.pdb')
    pastry('/bin/mv best.pdb.rebuilt '+prefix+'combo'+'%02d'%(i)+'.pdb.rebuilt')
    pastry('pulchra '+closcf)
    closcfr0=closcf+'.rebuilt'
    closcfr =prefix+'closc'+'%02d'%(i)+'.pdb.rebuilt'
    if closcfr0!=closcfr: pastry('/bin/mv '+closcfr0+' '+closcfr)

if tarname:
    pastry('/bin/tar jcf '+tarname+' '+prefix+'combo*.pdb.rebuilt '+prefix+'closc*.pdb.rebuilt')
    pastry('/bin/mv '+tarname+' '+outdir)
else:
    pastry('/bin/mv '+prefix+'combo*.pdb.rebuilt '+prefix+'closc*.pdb.rebuilt '+outdir)

sys.exit(0)
