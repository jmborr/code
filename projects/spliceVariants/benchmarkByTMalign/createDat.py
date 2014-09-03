#!/usr/bin/python

'''Go through every line of findCommonTemplate.dat, extract short,
long, templ

Align short and templ, extract isAligned line for short

Eliminate singlet alignet residues, ie, 010 goes to 000 (1:aligned,
0:unaligned) for internal, 10 goes to 00 for N-terminal, and 01 goes
to 00 for C-terminal

Force align of singlet unaligned residues, ie, 11011 goes to 11111 for
internal, 011 goes to 111 for N-terminal, 110 goes to 111 for
C-terminal

Create xxxxxrap3orienrev5s.pdb file from isaligned array, short
sequence, and template coordinates

Apply dat.x to create comb.dat combCA.dat dist.dat distL.dat files'''

import os,sys
from seq.letters import valid_A
from utilities.small_utilities import chomp,junkName,Bye
from utilities.codedir import codedir
from TMalign.TMalignManageResults import TMalignOut
from seq.alignYangManageResults import gappedSeq
from jobs.job import pastry
from mkinput.mkinput2.createDat import formatTempl

SCroot='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
templd='/gpfs1/archive/adrian/dat/pdb_mar07' #template directory
pdbd=templd+'/CA'
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'

def smooth_isAligned(idents):
    l=len(idents)
    #remove singlet internal aligned residues
    for i in range(1,l-1):
        if idents[i-1]==0 and idents[i+1]==0: idents[i]=0
    #remove singlet termini aligned residues
    if idents[0]==1 and idents[0]==0: idents[0]=0
    if idents[l-1]==1 and idents[l-2]==0: idents[l-1]=0
    #Error:I cannot force alignment of singlet unalignedresidues
    #because I don't have template coverage    
    #    #remove singlet unaligned residues
    #    for i in range(1,l-1):
    #        if idents[i-1]==1 and idents[i+1]==1: idents[i]=1
    #    #remove singlet termini aligned residues
    #    if idents[0]==0 and idents[0]==1: idents[0]=1
    #    if idents[l-1]==0 and idents[l-2]==1: idents[l-1]=1
    
junk=junkName()
pastry('/bin/cp '+codedir+'/f77/mkinput/mkinput2.0/dat.x .')
m=0 #variable m for debbuging purposes
#list=open(scd+'/toy.dat','r').readlines() #line for debugging purposes
list=open(scd+'/findCommonTemplate.dat','r').readlines()
for line in list:
    sys.stdout.write(line)
    if line[0]=='#': continue #line is a comment
    (short,long,templ)=line.split()

    #create subdirectory to store files
    subdir=scd+'/out/'+short+'_'+long ; os.system('mkdir -p '+subdir)
    
    spdb=pdbd+'/'+short+'.pdb' ; lpdb=pdbd+'/'+long+'.pdb' ; tpdb=pdbd+'/'+templ+'.pdb'
    pdbs={short:spdb,long:lpdb}

    #gather info on the alignment between short and long 
    cmd='TMalign '+spdb+' '+lpdb+' > '+junk  ;  os.system(cmd)
    
    alg=TMalignOut(junk)
    constitutive=alg.nAligned #number of aligned residues between long and short

    #read ATOM CA entries from template pdb file
    tpdbCA=[]
    for line in open(tpdb,'r').readlines():
        if 'ATOM' in line and ' CA ' in line: tpdbCA.append(line)

    for header in (short,long):
        #read 3-letter code for sequence from ATOM CA entries from short pdb file
        seqh=[]
        pdb=pdbs[header]
        for line in open(pdb,'r').readlines():
            if 'ATOM' in line and ' CA ' in line: seqh.append(line[17:20])
        #create subdirectory to store files
        subsubdir=subdir+'/'+header ; os.system('mkdir -p '+subsubdir)
        #align header to templ
        cmd='TMalign '+pdb+' '+tpdb+'  -L '+`constitutive`+' > '+junk  ;  os.system(cmd)
        alg2=TMalignOut(junk)
        alg2.initIsAlignedArray()
        nres=len(alg2.isAligned[1]) #number of residues in header sequence
        if len(seqh)!=nres:
            sys.stderr.write('ERROR: unequal number of residues for '+header+' in '+line)
        #smooth isAligned line
        smooth_isAligned(alg2.isAligned[1]) #pass alg2.isAligned[1] by reference
        #create headerrap3orienrev5s.pdb like-file
        iress=-1 #current residue in short sequence
        irest=-1 #current residue in templ sequence
        #navigate the gapped sequence for short, look if the position
        #is aligned according to the smoothed isAligned line, and if
        #so, retrieve the coordinates of the corresponding template
        #residue
        buf=''
        for i in range(len(alg2.seq[1])):
            if alg2.seq[2][i] in valid_A: #position in gapped-templ occupied by one aa
                irest+=1                
            if alg2.seq[1][i] in valid_A: #position in gapped-short occupied by one aa
                iress+=1
                if alg2.isAligned[1][iress]==1: #aa in short aligned to some aa in templ
                    line=tpdbCA[irest][0:54] #corresponding aa in templ
                    #thread aa from short into aa from templ
                    aa=line[17:20] ; n=line[22:26] ; xyz=line[30:54]
                    line='ATOM     1   CA  '+seqh[iress]+'%6d'%(iress+1)+'    '+xyz+' '+n+' '+aa+'\n'
                    buf+=line
	buf+='TER\n'
        open(subsubdir+'/'+header+'rap3orienrev5s.pdb','w').write(buf)
        #generate *.dat files with the created headerrap3orienrev5s.pdb file
        cmd='./dat.x -a 1 -b '+`nres`+' -c '+subsubdir+'/'+header+'rap3orienrev5s.pdb'
        pastry(cmd)
        #generate chain.dat with the created headerrap3orienrev5s.pdb file
        formatTempl(subsubdir+'/'+header+'rap3orienrev5s.pdb')
        #move output
        cmd='/bin/mv chain.dat comb.dat combCA.dat dist.dat distL.dat par.dat '+subsubdir
        pastry(cmd)

#    if m>10:Bye('') #for debbuging purposes
    m+=1 #for debbuging purposes

os.system('/bin/rm dat.x '+junk)
sys.exit(0)
