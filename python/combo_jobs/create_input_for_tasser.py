#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq
from jobs.job import gen_lf,pastry

inpHand('Usage: create_input_for_tasser.py',
        ' -i _R_header (five-letter code)',
        ' -o _AR_outd output directory (will append x/xxxxx/fasta, x/xxxxx/prospector and x/xxxxx/abinitio)',
        ' -a __co sequence identity cutoff (default=0.35)',
        ' -y __logrequest either _ASK_, or _LOG_ (def=_LOG_). _ASK_ and _LOG_ cause to automatically generate a logfile name with the identifier (flag -x), but _ASK_ will only output the logfile to STDOUT and exit. logfile will reside in the output directory',
        ' -x _R_identifier identifier string to generate the logfile', 
        ' -h __help this script will create the neccessary *.dat input files for tasser, and store them in appropriate directories'
        ).parse(locals(),sys.argv)

#generate a logfile, local to the working directory, but with absolute path
logfile=gen_lf(logrequest,outd+'/'+header[1]+'/'+header,identifier)

print '...running create_input_for_tasser.py...'

#set up default values for selected variables
if not co: co=35
else: co=100*eval(co)

currd=os.getcwd()

#create output directories
outd0=outd+'/'+header[1]+'/'+header
if not os.path.exists(outd0): os.system('mkdir -p '+outd0)

outdF=outd0+'/fasta'
if not os.path.exists(outdF): os.system('mkdir -p '+outdF)
outdR=outd0+'/prospector'
if not os.path.exists(outdR): os.system('mkdir -p '+outdR)
outdI=outd0+'/abinitio'
if not os.path.exists(outdI): os.system('mkdir -p '+outdI)

#create pair.dat and pair3.dat, remove unneccessary files
pastry(codedir+'/perl/tasser/tasser1.0/mkinput/mkpair.pl -i ' +header+' -o .')
pastry('/bin/mv pair1.dat pair3.dat '+outd0)
#create input files from FASTA homology
#..run prospector with only header as input (will use list of homologous files
#..as given by fasta filtering)
#os.system(codedir+'/python/prospector/prospector4.0/prospector.py -i '+\
#          header+' -o '+outdF+' -e _YES_')
#..obtain the restricted homology file
#os.system('/bin/cp '+outdF+'/'+header+'.homol .')
#os.system('/bin/cp '+outdF+'/'+header+'.threadrap3orienrev .')
#part=open(header+'.homol','r').readlines()
#threadf=open(header+'.threadrap3orienrev','r')
#threadf.readline()#discard first line
#for line in threadf.readlines():
#    values=line.split()
#    if eval(values[1]) > co:
#        part.append(values[0]+'\n')
#    else:
#        threadf.close()
#        break
#handle=open(header+'2.homol','w');
#handle.write(`len(part)-1`+'\n')
#handle.writelines(part[1:])
#handle.close()
#..purge unneccessary output files
#os.system('/bin/mv '+header+'2.homol '+outdR+header+'.homol')
#os.chdir(outdF)
#os.system('/bin/rm -f '+header+'.predictedrap3orienrev '+\
#          header+'.potcontactexpr1ap3orienrev '    +\
#          header+'rap3orienrev1.pdb '              +\
#          header+'.threadrap3orienrev '            +\
#          'finish_* list_* stat_* '                  )
#os.chdir(currd)
#..run mkdat, remove later uneccessary files
#os.system(codedir+'/python/tasser/tasser1.0/mkinput/mkdat.py -i '+\
#          header+' -o '+outdF)
#os.chdir(outdF)
#os.system('/bin/rm -f contact.dat template.dat')
#os.chdir(currd)

#create input files from PROSPECTOR homology
#..run prospector with homology file the restricted one recently created
#..also remove unneccessary files
#os.chdir(outdR)
#homolf=header+'.homol'
#os.system(codedir+'/python/prospector/prospector4.0/prospector.py -i '+\
#          header+' -o . -d '+homolf)
#os.system('/bin/rm -f '+header+'.predictedrap3orienrev '+\
#          header+'.potcontactexpr1ap3orienrev '    +\
#          header+'rap3orienrev1.pdb '              +\
 #          header+'.threadrap3orienrev '            +\
#          'finish_* list_* stat_* '                  )
#..run mkdat, remove later uneccessary files
#os.system(codedir+'/python/tasser/tasser1.0/mkinput/mkdat.py -i '+\
#          header+' -o .')
#os.system('/bin/rm -f contact.dat template.dat')
#os.chdir(currd)

#create input files from ABINITIO
#..run mkdat, remove later uneccessary files
os.system('echo "0" > contact.dat')
os.system('echo "0" > template.dat')
pastry(codedir+'/python/tasser/tasser1.0/mkinput/mkdat.py -a ./template.dat -b ./contact.dat -i '+header+' -o .')
#..remove all info except for exp.dat and par.dat
os.system('echo "0" >comb.dat;echo "0" >combCA.dat') #remove contact info
os.system('echo "0" >dist.dat;echo "0" >distL.dat') #remove dist-restraint info
os.system('echo "  0 hard" > chain.dat')            #remove template info
pastry('/bin/rm -f contact.dat template.dat seq.dat')
pastry('/bin/mv *dat '+outdI)

print '...finished create_input_for_tasser.py...'
sys.exit(0)
