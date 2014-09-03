#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq
from jobs.job import pastry
from random import randint
from utilities.small_utilities import chomp

ih=inpHand('Usage: tasser_data_done.py',
           ' -a __capdb create/import alpha-carbon-only pdb file',
           ' -b _AR_datad directory containing the *.dat files. It may also contain the alpha-carbon-only pdb and file rmsinp',
           ' -l _A_parexpd directory containing exp.dat and par.dat files',
           ' -c __tbz compress trajectories to file.tbz2 or just compress each trajectory separetely (option NOTAR) or compress to header[.iidp].tbz2 (option TAR). Will not compress if flag not specified',
           ' -d _A_indd in.dd file (def=/../tasser_rd/in.dd)',
           ' -e _A_paird directory containing pair1.dat and pair2.dat',
           ' -f _A_seqd directory containing the secondary structure file. If file "seq.dat" exist, it will copy this file, otherwise will look for xxxxx.SEQ',
           ' -g __nidp number of independent runs on the same node(default=1)',
           ' -i _R_header (five-letter code)',
           ' -j __iidp run ID if this run is part of a set of runs on a same target on different nodes',
           ' -k __move name of extra files to move to output directory',
           ' -n __noed directory containing the noe files, if needed',
           ' -o _AR_outd output directory',
           ' -t _A_tasserx tasser flavor (def=codedir/../original/cas)'
           )
ih.parse(locals(),sys.argv)

if not tasserx: tasserx=codedir+'/f77/tasser/tasser.1.0/original/cas'
if not indd: indd=codedir+'/f77/tasser/tasser.1.0/tasser_rd/in.dd'
if tbz:
    if not re.compile('NOTAR|TAR|\.tbz2$').search(tbz):
        ih.abort('compressed file must have .tbz2 extension');
if not nidp: nidp=1
else:
    nidp=int(nidp)
    if nidp<1: nidp=1
#create temporary directory to run tasser
root=os.getcwd(); tmpwd=root+'/tmp_tasser'+`randint(0,99999)`;
if not os.path.exists(tmpwd): os.makedirs(tmpwd)

#bring all input to the temporary working directory
#secondary structure propensities file
if seqd:
    if os.path.exists(seqd):
        if os.path.exists(seqd+'/seq.dat'):
            pastry('/bin/cp -f '+seqd+'/seq.dat '+tmpwd)
        elif os.path.exists(seqd+'/'+header+'.SEQ'):
            pastry('/bin/cp -f '+seqd+'/'+header+'.SEQ '+tmpwd+'/seq.dat')
        else:
            sys.stderr.write('ERROR: no seq.dat nor '+header+'.SEQ secondary propensities found in '+seqd+'\n')
            sys.exit(1)            
    else:
        sys.stderr.write('ERROR: no '+seqd+' found\n')
        sys.exit(1)
        
if paird: pastry('/bin/cp -f '+paird+'/pair*.dat '+tmpwd) #will exit if pair*.dat not found
pastry('/bin/cp -f '+tasserx+' '+tmpwd)
#  remove the directory path from the absolute filename of tasserx
tasserx=os.path.basename(tasserx)
pastry('/bin/cp -f '+indd+' '+tmpwd)
#  remove the directory path from the absolute filename of indd
indd=os.path.basename(indd)
pastry('/bin/cp -f '+datad+'/*.dat '+tmpwd)
if parexpd: pastry('/bin/cp -f '+parexpd+'/par.dat '+parexpd+'/exp.dat '+tmpwd)

#  deal with the noe data, if provided
if noed: pastry('/bin/cp -f '+noed+'/*.dat '+tmpwd)

pastry('/bin/cp -f '+codedir+'/f77/tasser/tasser.1.0/common.tbz2 '+tmpwd)#common input files

#  deal with rmsinp file
if header:
    pastry(codedir+'/python/tasser/tasser1.0/mkrmsinp.py -i '+header+' -s '+tmpwd+'/seq.dat -o '+tmpwd)
elif os.path.exists(datad+'/rmsinp'):
    pastry('/bin/cp -f '+datad+'/rmsinp '+tmpwd+'/rmsinp')
else:
    sys.stderr.write('Error: rmsinp file not in '+datad+' and I can\'t create it because you did not provide a header. Exiting now!')
    sys.exit(1)

#go to temporary working directory and start working!
os.chdir(tmpwd)

#  deal with the CA-only pdb file        
if capdb:
    if header:
        pastry(codedir+'/python/pdb/create_single_chain.py -i '+header+' -f '+capdb)
    elif os.path.exists(datad+'/'+capdb):
        pastry('/bin/cp -f '+datad+'/'+capdb+' '+tmpwd)
    else:
        stderr.write('Error: CA atom file '+capdb+' not in '+datad+' and I can\'t create it because you did not provide a header. Exiting now!')
        sys.exit(1)

#  uncompress *.comm files
pastry('/bin/tar jxf common.tbz2')

#  run tasser
print '...tasser_data_done.py calling',tasserx
print '...running',tasserx
if nidp==1: #do a single run
    runID=''
    if iidp: #we passed an identifier
        tasserx=tasserx+' runID='+iidp
        runID='%02d' % (int(iidp))
    #before running tasser, copy hostname and workding dir to hostfile
    hostf='/gpfs1/active/casp7/Host/jose/'+header+'.ab.'+runID
    open(hostf,'w').write(chomp(os.popen('hostname').readline())+' '+os.getcwd()+'\n')
    #run tasser
    pastry('./'+tasserx)
    pastry("ls -1 rep*tra | wc | awk '{print $1}' > tra"+runID+".in")
    pastry("ls -1 rep*tra >> tra"+runID+".in")
else: #do nidp runs on same node
    indds=open('in.dd','r').readlines() #load whole in.dd file into memory, keep newline chars
    p=re.compile('(\d+)')
    seed=int(p.match(indds[0]).group(1))-1 #find seed
    for n in range(0,nidp):
        seed=seed+1 #update seed
        id='%02d'%(n)
        indds[0]=p.sub(`seed`,indds[0],1) #change seed in the first line of indds
        open('in.dd','w').writelines(indds) #rewrite in.dd file with new seed
        pastry('./'+tasserx)   #run tasser
        reptras2=[]
        reptras=chomp(os.popen('ls -1 rep?.tra').readlines()) #load all created trajectories
        reptras=reptras+chomp(os.popen('ls -1 rep??.tra').readlines())
        for reptra in reptras:
            reptra2=reptra[:-4]+'_'+id+'.tra'
            pastry('mv '+reptra+' '+reptra2)#rename trajectoy files
            reptras2.append(reptra2+'\n')
        #create tra_id.in file that list the rep*_id.tra files
        train=open('tra_'+id+'.in','w')
        train.write(`len(reptras2)`+'\n')
        train.writelines(reptras2)
        train.close()
    pastry("ls -1 rep*tra | wc | awk '{print $1}' > tra.in") #gather all in tra.in
    pastry("ls -1 rep*tra >> tra.in")
print '...finished',tasserx

#save all output
#  maybe we'll compress trajectory files
if tbz:
    if tbz=='NOTAR':
        print '...tasser_data_done.py calling bzip2...'
        print '...running bzip2...'
        pastry('/usr/bin/bzip2 rep*')
        print '...finished bzip2...'
    elif tbz=='TAR':
        if runID: tbz=header+'.'+runID+'.tbz2'
        else:     tbz=header+'.tbz2'
        pastry('/bin/tar jcf '+tbz+' rep*.tra ; /bin/mv '+tbz+' '+outd)
    else: #passed the name of a file
        pastry('/bin/tar jcf '+tbz+' rep*.tra ; /bin/mv '+tbz+' '+outd)


if capdb: pastry('/bin/mv -f '+capdb+' '+outd)

#move some input and also the output
pastry('/bin/mv -f *.dat rmsinp rep* '+indd+' out.d tra*.in '+outd)

#extra files, if required, are moved also
if move: pastry('/bin/mv -f '+move+' '+outd)

#remove temporary directory
os.chdir(root)
pastry('/bin/rm -rf '+tmpwd)

print '...finished tasser_data_done.py...'
sys.exit(0)
