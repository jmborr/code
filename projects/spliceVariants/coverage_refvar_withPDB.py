#!/usr/bin/python

import os,sys
from string import lowercase
from utilities.codedir import scratchdir
from utilities.small_utilities import Bye,junkName,chomp
from seq.align import getYangManageResultsAlignmentFormat
from seq.fastaManager import importFastaEntry

if os.popen('hostname').readline().strip() != 'cng0002':
    sys.stderr.write('You must be on cng0002 to run this script\n')
    sys.exit(1)
    
pin=open( os.path.join(scratchdir,'spliceVariants/preparing_database/refs300_nonduplicated_withPDB.dat'),'r')
l=pin.readline()

junkf=junkName()

while l and l[0:2]=='AC':
    ref=l[3:11]
    root=ref[0:6]
    drs=[]
    pdbheaders=[]
    l=pin.readline()
    while l and l[0:2]=='DR': #gather all annotated PDB's for this sequence
        drs.append(l)
        pdbheaders.append(l[10:14].lower())
        l=pin.readline()
    cmd='ls -1 '+scratchdir+'/spliceVariants/input/'+root+'*'
    for ref_fasta in chomp(os.popen(cmd).readlines()):
        ref=ref_fasta[-8:]
        print 'AC',ref
        ref_header,ref_seq=importFastaEntry(open(ref_fasta,'r'))
        for i in range(len(pdbheaders)):
            pdbheader=pdbheaders[i]
            dr=drs[i]
            pdbf='/local/images/pdb-2007062900/pdb/pdb'+pdbheader+'.ent'
            os.system('pdbseq -o '+junkf+' '+pdbf+' &>/dev/null') #extract all chains in the PD
            pin2=open(junkf,'r') ; l2=pin2.readline()
            id=0.0 ; max_id_cov=0.0
            while l2:
                if l2.find('sequence')>=0:
                    chain_seq=pin2.readline().strip()
                    try:
                        x=getYangManageResultsAlignmentFormat(ref_seq,chain_seq)
                        id_cov=float(x.identicalL())/len(ref_seq)
                    except:
                        id_cov=0.0
                    if id_cov>max_id_cov:
                        max_id_cov=id_cov
                        id=x.seqIdentity()
                l2=pin2.readline()
            print dr.strip(),'id=%4.2lf, id_cov=%4.2lf'%(id,max_id_cov)

os.system('/bin/rm '+junkf)

sys.exit(0)
