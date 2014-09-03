#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq,pdbnonhom,pdbhomolrestricted
from seq.length import length_header
from random import randint

ih=inpHand('Usage: prospector.py',
           ' -i _R_header five letter pdb code',
           ' -a _A_closemsa 3590.aln type file',
           ' -b _A_distmsa e10.aln type file',
           ' -c _A_secon header.SEQ type file for predicted secondary struct',
           ' -d _A_homol header.homol type file, to exclude homologous',
           ' -f __homoltype homology flavor, either _RESTRICTED_ or _FASTA_ (def=_FASTA_). Restricted homology definition is with prospector, while normal homology definition if with psi-blast',
           ' -e __gendat either of _YES_ or _NO_, to generage contact.dat and template.dat, suitable for mkdat input (def=_NO_)',
           ' -g _A_fastaf fasta file, needed to generate msa if not found/provided, and if fasta file cannot be generated from header file',
           ' -h __help display this message',
           ' -o _AR_outd output directory',
           )
ih.parse(locals(),sys.argv)
#create and go to temporary working directory
workd=os.getcwd()+'/tmp_prospector'+`randint(0,1000000)` #;print workd
os.system('mkdir -p '+workd)
os.chdir(workd)

#create other neccessary directories, and bring common files
os.system('mkdir -p inputs/3590; mkdir -p inputs/e10; mkdir -p inputs/secon')
os.system('/bin/cp -rf '+codedir+'/f77/prospector/prospector4.1/common/common.tbz2 .')
os.system('tar jxf common.tbz2')
os.system('mkdir homologousListd')

os.system('echo "'+workd+'">./bin/homedir')
closemsad=workd+'/inputs/3590/'
os.system('echo "'+closemsad+'" >./bin/3590')
distmsad=workd+'/inputs/e10/'
os.system('echo "'+distmsad+'" > ./bin/e10')
second=workd+'/inputs/secon/'
os.system('echo "'+second+'" > ./bin/secon')

#generate file LIST
os.system('echo "1" > ./bin/LIST')
os.system('echo "'+header+' '+`length_header(header)`+'" >> ./bin/LIST')

#generate/obtain close multiple alignment
finalf=closemsad+header+'.aln'
if closemsa: os.system('/bin/cp -f '+closemsa+' '+finalf)
else:
    query=pdbnonhom+'/3590/'+header+'.aln' #;print query;sys.exit(0)
    if os.path.exists(query):
        os.system('/bin/cp -f '+query+' '+finalf)
    else: #generate close multiple structure alignment
        if not fastaf: #generate fasta file from header name
            if os.system(codedir+'/seq/mkfasta.py -i '+header+\
                         ' -o '+header+'.fasta'):
                ih.abort('ERROR: can\'t generate fasta file from header')
            os.system(codedir+'/bash/msa/sjfasta '+header+'.fasta')
            os.system('/bin/mv -f '+header+'close.aln '+finalf)
            finalf2=distmsad+header+'.aln'
            os.system('/bin/mv -f '+header+'dist.aln '+finalf2)
os.system('/bin/mv '+finalf+' ./bin/3590/')            
print 'Bye!';sys.exit(0)
#generate/obtain distant multiple alignment
finalf=distmsad+header+'.aln'
if distmsa:
    os.system('/bin/cp -f '+distmsa+' '+finalf)
elif os.path.exists(finalf): #finalf generated in the previous section
     pass #nothing to do
else:
    query=pdbnonhom+'/e10/'+header+'.aln'
    if os.path.exists(query):
        os.system('/bin/cp -f '+query+' '+finalf)
    else: #generate distant multiple structure alignment
        if not fastaf: #generate fasta file from header name
            if os.system(codedir+'/seq/mkfasta.py -i '+header+\
                         ' -o '+header+'.fasta'):
                ih.abort('ERROR: can\'t generate fasta file from header')
            os.system(codedir+'/bash/msa/sjfasta '+header+'.fasta')
            os.system('/bin/mv -f '+header+'dist.aln '+finalf)
os.system('/bin/mv '+finalf+' ./bin/e10/')

#generate/obtain secondary structure alignment
finalf=second+header+'.dat'  #is it really .dat extension ??
if secon: os.system('/bin/cp -f '+secon+' '+finalf)
else:
    query=pdbnonhom+'/seq/'+header+'.SEQ'
    if os.path.exists(query):
        os.system('/bin/cp -f '+query+' '+finalf)
    else: #generate secondary structure assignment
        cmd=codedir+'/python/tasser/tasser1.0/mkinput/mkseq.py -i '+header
        if fastaf: cmd=cmd+' -f '+fastaf
        os.system(cmd)
os.system('/bin/mv '+finalf+' ./bin/secon')

#generate/obtain file with list of homologous sequences
finalf=workd+'/homologousListd/'+header+'.homol'
if homol:
    os.system('/bin/cp -f '+homol+' '+finalf)
else:
    if homoltype=='_RESTRICTED_':
        query=pdbhomolrestricted+'/'+header+'.homol'
        if os.path.exists(query):
            os.system('/bin/cp -f '+query+' '+finalf)
        else:
            homoltype=='_FASTA_'
    if homoltype=='_FASTA_': #generate homologous<
        system(codedir+'perl/seq/find_homologs.pl -i '+header)
os.system('/bin/mv '+finalf+' ./homologousListd/')
homol=workd+'/homologousListd/'+header+'.homol'

#run prospector
os.chdir('./bin')
os.system('./direct')
os.system('./list3')
os.system('./run_prospector')
os.system('./summary')

#generate mkdat input
if gendat=='_YES_': #generate contact.dat and template.dat for mkdat input
    os.system(codedir+'/perl/tasser/tasser1.0/mkinput/create_contact.dat.pl -i '+header)
    os.system(codedir+'/perl/tasser/tasser1.0/mkinput/create_template.dat.pl -i '+header)
    os.system('/bin/mv -f contact.dat template.dat '+outd)

#move output to output directory
os.system('/bin/mv -f '+header+'.predictedrap3orienrev '+\
          header+'.potcontactexpr1ap3orienrev '    +\
          header+'rap3orienrev1.pdb '              +\
          header+'.threadrap3orienrev '            +\
          'finish_* list_* stat_* '                +\
          homol+' '                                +\
          outd)
