#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand,does_exists
from utilities.codedir import codedir,libseq
from utilities.small_utilities import chomp,junkName,Bye
from random import randint
from jobs.job import pastry
from seq.alignYangManageResults import alignOut

ih=inpHand('Usage: prospector.py [options]',
           ' -a _R_header header of the query',
           ' -b _A_profdir root input for query. Subdirectories assumed are 3590, e10, profile, and seq (def:/gpfs1/archive/adrian/dat/splice)',
           ' -c _A_templlib future directory in node where library will be located (def: /tmp/jose/templatedir)',
           ' -d _A_execZ compressed files containing executables and scripts used during prospector run (def: codedir/f77/prospector/prospector_genomesmay06/modified/exec.tbz2)',
           ' -e _A_outd outdir (default: current dir)',
           ' -f _A_benchlist list of templates included in execZ (def: LIST.may2006good35 )',
           ' -g __co remove templates with higher global pairwise sequence identity (def: co=1.0)',
           ' -i _A_fastaf fasta file, needed if we pass co',
           ' -j __tbz2f basename of .tbz2 file, if you want to compress results'
           )
ih.parse(locals(),sys.argv)

if not profdir : profdir ='/gpfs1/scratch/adrian/dat/splice'
if not templlib: templlib='/tmp/jose/templatedir'
if not execZ : execZ='/gpfs1/active/jose/code/f77/prospector/prospector_genomesmay06/modified/exec.tbz2'
if not outd  : outd =os.getcwd()
if not os.path.exists(outd): os.system('/bin/mkdir -p '+outd)
for file in (profdir,templlib,execZ): does_exists(file) #check if input exists
if not co: co=1.0

if not os.path.exists(templlib):
    Bye('ERROR from prospector.py: '+templlib+' does not exists')
    
strlist='/gpfs1/active/jose/code/f77/prospector/prospector_genomesmay06/strlist/'
homedir=junkName()  ;  pastry('/bin/mkdir '+homedir) #working dir

#current, workding directories
os.chdir(homedir)  #current dir is home dir
homedir=os.getcwd()#absolute path
wd='workdir'
pastry('/bin/mkdir '+wd)
os.chdir(wd)       #current dir is working directory

print 'create directories to store profiles and secondary structure files.'
querydir={}
for key in ('3590','e10','profile','seq'):
    querydir[key]=homedir+'/'+key
    pastry('/bin/mkdir -p '+querydir[key])
querydir['templatedir']=templlib

print 'create files containing path to previous directories'
dict={'3590' :querydir['3590'],'e10':querydir['e10'], 'mtx':querydir['profile'],
      'secon':querydir['seq'] ,'templatedir':querydir['templatedir']}
for key in dict.keys():
    dir=dict[key]
    pastry('echo '+dir+'/ > '+key) #notice last "/"

print 'files containing file extensions'
extensions={'3590aln':'aln','e10aln':'aln','seconend':'SEQ'}
for key in extensions.keys():
    value=extensions[key]
    pastry('echo .'+value+' > '+key)   #files storing extension, prepend '.'

print 'Creage LIST, input, input2, namesize, homedir files'
pastry('echo -e "1\n'+header+'" > LIST')          #create LIST  file. Contains only query header
pastry('echo pdbb > input && echo pdb > input2')  #create input and input2 files
pastry('echo '+`len(header)`+' > namesize')       #create namesize file
pastry('echo '+homedir+'/ > homedir')             #create homedir  file

print 'expand execZ, copy strlist'
pastry('/bin/cp '+execZ+' . && tar jxf '+execZ) #list of template headers, executables, scripts
pastry('/bin/cp -rf '+strlist+' '+homedir)      #strlist directory under homedir

#remove header from the list benchlist
if not benchlist: benchlist='LIST.may2006good35'
pin=open(benchlist,'r') ; pin.readline() ; line=chomp(pin.readline()) ; hs={}
while line:
    h,l=line.split() #template header and corresponding sequence length
    if header!=h: hs[h]=l #do not include the query header
    line=chomp(pin.readline())
#remove homologs if we pass some sequence cut-off
if co<1.0:
    alignx=junkName() ; pastry('/bin/cp $HOME/bin/align '+alignx) #bring copy of align to node
    co=float(co)
    for h in hs.keys():
        templfasta=templlib+'/input/'+h        
        try:
            os.stat(templfasta)
        except:
            Bye('ERROR in prospectory.py: template fasta file "'+templfasta+'" does not exists') 
        if not os.path.exists(fastaf):
            Bye('ERROR in prospectory.py: query fasta file " '+fastaf+'" does not exists') 
        cmd=alignx+' '+fastaf+' '+templfasta+' 0' #global pairwise alignment
        al=alignOut( os.popen(cmd).readlines() )
        if al.readError!='':
            sys.stderr.write(h+' NOT ALIGNED to query sequence\n'+al.readError)
            continue
        if float(al.seqid) > co: del hs[h] #the template is a homolog
#rewrite benchlist file
buf='%d\n'%len(hs)
for h in hs.keys(): buf+=h+'%5s\n'%(hs[h])
open('./LIST.templates','w').write(buf)

print 'copy query profiles and secondary structure prediction'
pastry('/bin/cp '+profdir+'/seq/'    +header+'.SEQ '+querydir['seq']    )
pastry('/bin/cp '+profdir+'/profile/'+header+'.mtx '+querydir['profile'])
pastry('/bin/cp '+profdir+'/3590/'   +header+'.aln '+querydir['3590']   )
pastry('/bin/cp '+profdir+'/e10/'    +header+'.aln '+querydir['e10']    )

print 'run prospector...'
pastry('./direct && ./list3')
for executable in ('pairorpsigly3','r5saorienrev','r5sap2orienrev','r5sap3orienrev'):
    sys.stdout.write('running '+executable+'...\n')
    pastry('./'+executable+' < nametarg1') #run each at a time
    sys.stdout.flush()
    sys.stderr.flush()

print 'assess results with summary...'
pastry('./summary')

pastry('/bin/mkdir -p '+homedir+'/pdbbsummary')
pastry('/bin/mv stat_* list* STRUCTURE_* '+homedir+'/pdbbsummary')
os.chdir(homedir)

print 'move results to '+outd
if outd!=homedir:
    if(tbz2f):
        pastry('tar jcf '+tbz2f+' pdbb* && /bin/mv '+tbz2f+' '+outd)
    else:
        pastry('/bin/rm -rf '+outd+'/pdbb* && /bin/mv pdbb* '+outd)

sys.exit(0)
