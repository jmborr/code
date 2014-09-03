#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import codedir,libseq
ih=inpHand('Usage: tasser1.0_header.py',
           ' -i _R_header (five-letter code)',
           ' -t _A_tasserx tasser flavor (def=codedir/../tasser_rd/cas_rd)',
           ' -d _A_indd in.dd file (def=/../tasser_rd/in.dd)',
           ' -a __capdb create alpha-carbon-only pdb file',
           ' -o _AR_outd output directory',
           ' -c __tbz compress trajectories to file.tbz2 or just compress each trajectory separetely (option NOTAR)',
           ' -h __help this script runs tasser with header being the only input'
           )
ih.parse(locals(),sys.argv)
if not tasserx: tasserx=codedir+'/f77/tasser/tasser.1.0/tasser_rd/cas_rd'
if not indd: indd=codedir+'/f77/tasser/tasser.1.0/tasser_rd/in.dd'
if tbz:
    if not re.compile('NOTAR|\.tbz2$').search(tbz):
        ih.abort('compressed file must have .tbz2 extension');

root=os.getcwd(); tmpwd=root+'/tmp_tasser';
if not os.path.exists(tmpwd): os.makedirs(tmpwd)
#bring all input to the temporary working directory
os.system('/bin/cp -f '+codedir+'/python/tasser/tasser1.0/mkinput/mkdat.py '+
          tmpwd)
os.system('/bin/cp -f '+codedir+'/perl/tasser/tasser1.0/mkinput/mkpair.pl '+
          tmpwd)
os.system('/bin/cp -f '+codedir+'/f77/tasser/tasser.1.0/original/common.tbz2 '+
          tmpwd)
os.system('/bin/cp -f '+indd+' '+tmpwd) ; indd=os.path.basename(indd)
os.system('/bin/cp -f '+tasserx+' '+tmpwd)
tasserx=os.path.basename(tasserx)

#go to temporary working directory and start working!
os.chdir(tmpwd)
os.system('/bin/tar jxf common.tbz2')

#create CA-only chain, if required
if capdb:
    os.system(codedir+'/python/pdb/create_single_chain.py -i '+header+' -f '+capdb)
#run some scripts
print '...tasser1.0_header.py calling mkdat.py'
os.system('./mkdat.py -i '+header+' -o .')

print '...tasser1.0_header.py calling mkpair.pl'
os.system('./mkpair.pl -d . -o . 2>/dev/null')

os.system('echo 0 > comb.dat; echo 0 > combCA.dat') #remove contact info
os.system('echo 0 > dist.dat; echo 0 > distL.dat')  #remove dist-restraint info
os.system('echo "  0 hard" > chain.dat')            #remove template info

print '...tasser1.0_header.py calling mkrmsinp.py'
os.system(codedir+'/python/tasser/tasser1.0/mkrmsinp.py -i '+header+
          ' -f seq.fasta -o .')
#run tasser
print '...tasser1.0_header.py calling',tasserx
print '...running',tasserx
os.system('./'+tasserx)
print '...finished',tasserx
#create tra.in file
os.system("ls -1 rep*tra | wc | awk '{print $1}' > tra.in")
os.system("ls -1 rep*tra >> tra.in")
#save all output and remove temporary working directory
if tbz:
    if tbz=='NOTAR':
        print '...tasser1.0_header.py calling bzip2...'
        print '...running bzip2...'
        os.system('/usr/bin/bzip2 rep* ; /bin/mv -f rep* '+outd)
        print '...finished bzip2...'
    else:
        print  '...tasser1.0_header.py calling tar jcf...'
        print '...running tar jcf...'
        os.system('/bin/tar jcf '+tbz+' rep*.tra ; /bin/mv '+tbz+' '+outd)
        print '...finished tar jcf...'
else: os.system('/bin/mv rep*.tra '+outd)

if capdb: os.system('/bin/mv -f '+capdb+' '+outd)

os.system('/bin/mv -f seq.fasta *.dat pair.1 pair.3 rmsinp '+indd+' out.d tra.in '+outd+' 2>/dev/null')
os.chdir(root)

os.system('/bin/rm -rf '+tmpwd)

print '...finished tasser1.0_header.py...'
