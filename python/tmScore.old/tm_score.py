#!/usr/bin/python

import sys,os
from inputArgs.inputArgs import inpHand
from tmScore.tm_score_cpp import allocI1,allocD1,allocD2,lowerNotAligned,getCAcoords2,getrmsd,gettm,coverage,rmsd_cov,supCAs
from seq.get_seq import parse_the_pdbfile
from random import randint
from utilities.small_utilities import junkName,Bye
from utilities.codedir import codedir
from seq.alignYangManageResults import alignOut
from seq.letters import valid_A

N=2000 #maximum number of amino acids

def TMscore(f1,f2):
    #print f1,f2
    cmd=codedir+'/f77/TM-score/TMscore '+f1+' '+f2+' | grep "TM-score    ="'
    return float(os.popen(cmd).readline().split()[2])

def processSeq(pdb1,alg):
    line='REMARK        10|       20|       30|       40|       50|       60|       70|\n'
    pdbf=open(pdb1,'r')
    seq2=parse_the_pdbfile(pdbf)    #;print seq2 #get sequence in one line*/
    seq=lowerNotAligned(seq2,alg)   #lower case those residues above d0

    while seq:
        line=line+'REMARK '+seq[0:70]
        seq=seq[70:]
        if seq: line=line+'\n'
    return line

#default value d0=0 means we will take d0 as given by scaling formula
def getVerbose(x1,x2,n,pdb1,tx='',ty='',U='',d0=0.):
    if not tx: tx=allocD1(3)
    if not ty: ty=allocD1(3)
    if not U :  U=allocD2(3,3)
    rmsd=getrmsd(x1,x2,n) #; print rmsd
    tm=gettm(x1,x2,n,d0,tx,ty,U) #;print tm
    alg=allocI1(n); 
    cov=coverage(x1,x2,n,tx,ty,U,alg,d0) #;print cov
    seq=processSeq(pdb1,alg) #return seqence, uppercase characters only for aligned residues
    rmsdcov=rmsd_cov(x1,x2,n,tx,ty,U,d0)
    line='REMARK lseq= %4d rmsd= %5.2f tm= %5.3f cov= %4.2f rmsd_cov= %5.2f\n' \
          %(n,rmsd,tm,cov,rmsdcov)
    line=line+seq #add the unaligned/aligned residues sequence
    return [tm,line,tx,ty,U]

#d0=0.0 means we will take d0 from the scaling formula
def gettmCA(pdb1,pdb2,sup='',verbose='',d0=0.0,seqd='',aligned_pairs=[]):
    for file in (pdb1,pdb2):
        if not os.path.exists(file):
            sys.stderr.write('ERROR: file '+file+' does not exists\n')
            return 0
    #allocate big arrays to hold coordinates
    x1=allocD2(N,3) 
    x2=allocD2(N,3)
    #if seq.dat passed, do TM-score only of residues assigned as helix
    #and/or strand. The d0 correspond to all residues, not just those
    #of the helix/strand assigned. The reason is that the expected
    #distance between two corresponding residues after a random
    #alignment is determined by the whole sequence length
    if seqd:
        if not os.path.exists(seqd):
            sys.stderr.write('ERROR: '+seqd+' does not exist!\n')
            sys.exit(1)
        #j1,j2 temporary files holding only those residues assigned as helix or strand
        j1='junk'+`randint(0,999)`  ; j2=j1
        while(j1==j2): j2='junk'+`randint(0,999)` #be sure j2 != j1
        pj1=open(j1,'w')    ;  pj2=open(j2,'w')
        pd1=open(pdb1,'r')  ;  pd2=open(pdb2,'r')
        #we assume that residue order for pdb1,pdb2, and seqd is THE SAME
        n=0 ; m=0
        for line in open(seqd,'r').readlines():
            n+=1
            l1=pd1.readline()
            while l1.find(' CA ') < 0: l1=pd1.readline() #care only for CA
            l2=pd2.readline()  ;
            while l2.find(' CA ') < 0: l2=pd2.readline()
            items=line.split()
            sec_code=int(items[2])
            if sec_code != 1: #residue is not assigned as coil
                m+=1
                pj1.write(l1)
                pj2.write(l2)
        pj1.close()
        pj2.close()
        d0=1.24*(n-15)**(1./3)-1.8  #; print n,m
        #store CA coordinates
        n =getCAcoords2(x1,j1)  
        n2=getCAcoords2(x2,j2) #; print n,n2
        os.system('/bin/rm '+j1+' '+j2)
    else:
        n =getCAcoords2(x1,pdb1)  
        n2=getCAcoords2(x2,pdb2) #; print n,n2
    if n!=n2:
        sys.stderr.write("ERROR: pdb files don't have same number of residues")
        sys.exit(1)
    if n>N:
        sys.stderr.write("ERROR from tm_score.py: chain longer than maximum!")
        sys.exit(1)
    #allocate extra stuff if we want superposition
    if sup:
        os.system('/bin/rm '+sup+' 2>/dev/null')
        [tm,line,tx,ty,U]=getVerbose(x1,x2,n,pdb1,d0=d0)
        supCAs(pdb1,sup,x1,x2,n,tx,ty,U) #pdb1 is used for amino acid composition
        open(sup,'a').write(line+'\n') #append REMARKs at end of file
        if verbose=='y' or verbose=='Y' or verbose=='yes' or verbose=='YES':
            return line
        else: return tm
    if verbose=='y' or verbose=='Y' or verbose=='yes' or verbose=='YES':
        return getVerbose(x1,x2,n,pdb1,d0=d0)[1]
    return gettm(x1,x2,n,d0)

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: tm_score.py [options]',
               '  -a _AR_pdb1 first structure',
               '  -b _AR_pdb2 second structure',
               '  -c _A_sup output file with superposition of the two structures',
               '  -d __verbose output TM,rmsd,cov to standar output',
               '  -e __d0 d0 parameter (default: d0=0.0, which means d0 will be taken from the scaling formula',
               '  -f _A_seqd do TM-score of predicted sec.str. residues only, given by seq.dat',
               '  -g __chunkl list of chunks on which to do TM-score (syntax: 1-34,67-90)',
               '  -i _A_alignf file with alignment in Yang\'s format ',
               '  -j _A_identflag if passed file with alignment, specify which residues to extract: aligned residues("aligned"), or identical residues("identical:) (def: "aligned")',)
    ih.parse(locals(),sys.argv)
    if not d0:d0=0.0
    else:  d0=float(d0)
    if chunkl: #only pass the chunks to gettmCA
        chunks=chunkl.split(',')
        #dump the chunks into junk pdb file
        cas=os.popen('grep -e " CA " '+pdb1+' | grep "ATOM "').readlines()
        pdb3=junkName() ; pt=open(pdb3,'w')
        begin,end=chunks[0].split('-') ; begin=int(begin)-1 ; end=int(end)
        pt.writelines(cas[begin:end]) ; pt.close()
        #dump the chunks into junk pdb file        
        cas=os.popen('grep -e " CA " '+pdb2+' | grep "ATOM "').readlines()
        pdb4=junkName() ; pt=open(pdb4,'w')
        begin,end=chunks[1].split('-') ; begin=int(begin)-1 ; end=int(end)
        pt.writelines(cas[begin:end]) ; pt.close()

        print gettmCA(pdb3,pdb4,sup,verbose,d0=d0,seqd=seqd)
        os.system('/bin/rm '+pdb3+' '+pdb4) #some cleanup
    elif alignf:
        cas1=os.popen('grep -e " CA " '+pdb1+' | grep "ATOM "').readlines()
        cas2=os.popen('grep -e " CA " '+pdb2+' | grep "ATOM "').readlines()
        alg=alignOut(alignf).alg ; gs1=alg.gs[1]  ; gs2=alg.gs[2]
        pdb3=junkName() ; pt3=open(pdb3,'w')
        pdb4=junkName() ; pt4=open(pdb4,'w')
        #find aligned residues
        for i in range(gs2.gsL):
            if gs1[i] in valid_A and gs2[i] in valid_A:
                if identflag=='identical' and gs1[i]!=gs2[i]: continue
                j=gs1.gsi2si[i] ; pt3.write(cas1[j])
                k=gs2.gsi2si[i] ; pt4.write(cas2[k])
        pt3.close() ; pt4.close()
        print gettmCA(pdb3,pdb4,sup,verbose,d0=d0,seqd=seqd)
        os.system('/bin/rm '+pdb3+' '+pdb4) #some cleanup
    else: print gettmCA(pdb1,pdb2,sup,verbose,d0=d0,seqd=seqd) #sup='' by default

    
    sys.exit(0)
