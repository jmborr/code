#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.small_utilities import three2one #three2one dict for aa codes

ih=inpHand('Usage: mkraw.py  -i|-p|-f -c -o -d',
        '  -i __pdbh five-letter pdb header',
        '  -p _A_pdbf pdb file',
        '  -f _A_fastaf fasta file',
        '  -c __chId chain id',
        '  -o __outf output file with raw sequence (def=\'seq.raw\')',
        '  -d __outd output dir (def=curr dir)',
        '  -r _A_rootd root directory, if outd=\'RHEADER\', rootd/2/1234',
        '  -h __help Examples:\n-i 101m_, fetch /library/pdb/pdb101m.ent, output ./seq.raw\n-i 101m_ -o 101m_raw, fetch /library/pdb/pdb101m.ent, output ./101m_.raw\n-f fasta, fetch ./fasta, output ./seq.raw\n-i 101m_ -f fasta, same as before\n-i 101m_ -f fasta -o 101m_raw , fetch ./fasta, output ./101m_.raw\n-i 101m_ -f fasta -d HEADER, fetch ./fasta, output ./0/101m/seq.raw\n-i 101m_ -f fasta -o 101m_raw -d HEADER, fetch ./fasta, output ./0/101m/101m_.raw\n-i 101m_ -f fasta -o 101m_raw -d R_HEADER -r /tmp/jose/../, fetch ./fasta, output /tmp/0/101m/101m_.raw\n-i 101m_ -o 101m_raw -d R_HEADER -r /tmp/jose/../, fetch /library/pdb/pdb101m.ent, output /tmp/0/101m/101m_.raw\n-i 101m_ -p this.pdb, fetch ./this.pdb, output ./seq.raw\n-i 102dB -p this.pdb, fetch ./this.pdb, look for chain B, output ./seq.raw\n-i 102dB -c A -p this.pdb, fetch ./this.pdb, look for chain A, output ./seq.raw\n-c A -p this.pdb, same as before'
        )        

#special lines, since not a single flag is allways required
ih.parse(locals(),sys.argv)
if (not pdbh) and (not pdbf) and (not fastaf):
    ih.abort('At least one of (-i|-p|-f) flags needed!')

#parse input flags
if pdbh: #we passed a pdb header
    if not chId:  chId=pdbh[-1]
    if (not pdbf) and (not fastaf):
        pdbf='/library/pdb/pdb'+pdbh[0:4]+'.ent'
if (not chId) | (chId=='_'): chId=' '

#parse output flags
if not outf:
    outf='seq.raw' #default output file name
if outd=='HEADER':
    if not pdbh: ih.abort('!!! pdb header needed !!!')
    outd=addAbsPath('./')+'/'+pdbh[1:2]+'/'+pdbh
elif outd=='R_HEADER':
    if not pdbh: ih.abort('!!! pdb header needed !!!')
    if not rootd: ih.abort('Need flag -r rootd if -d R_HEADER')
    outd=addAbsPath(rootd)+'/'+pdbh[1:2]+'/'+pdbh
elif outd:
    outd=addAbsPath(outd);
if not outd:
    outd=addAbsPath('./')

outf=outd+'/'+outf
#print "pdbh=",pdbh,"\npdbf=",pdbf,"\nfastaf=",fastaf,"\noutf=",outf,"\noutd=",outd,"\nchId=",chId
if not os.path.exists(outd): #create output directory, if non-existent
    os.makedirs(outd)

outl=''
if pdbf: #get raw sequence from PDB file
    if not os.path.exists(pdbf): ih.abort(' !!!Non existent PDB file '+pdbf)
    #pick first chId if chId is blank space
    if chId==' ':
        for line in open(pdbf,'r').readlines():
            if line[0:5]=='ATOM ':
                chId=line[21]
                break
    #go to appropriate line in pdb file
    pdbhandle=open(pdbf,'r')
    line=pdbhandle.readline()
    while line[0:5]!='ATOM ' and line[21]!=chId:
        line=pdbhandle.readline()
    #create the single-symbol sequence
    while line[0:4]!='TER ' and line[0:4]!='END ':
        if line[0:5]=='ATOM ' and line[12:16]==' CA ':
            outl=outl+three2one[line[17:20]]
        line=pdbhandle.readline()
    pdbhandle.close()
elif fastaf: #get raw sequence from FASTA file
    try: fastah=open(fastaf,'r')
    except IOError: ih.abort(' !!! fasta file '+fastaf+' non-existent')
    fastah.readline()
    for line in fastah.readlines():
        if line[len(line)-1] == "\n": line=line[:-1]
        outl=outl+line
open(outf,'w').write(outl) #handle to file        

