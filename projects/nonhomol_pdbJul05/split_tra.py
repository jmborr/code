#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from jobs.job import pastry
from random import randint
from utilities.small_utilities import chomp

inpHand('Usage: split_tra.py.py',
        ' -a _A_ind input directory (def=./)',
        ' -b __nsplit split each rep*.tra in nsplit number of times (def=10)'
        ).parse(locals(),sys.argv)

if not ind: ind=os.getcwd();
if not nsplit: nsplit=10
else: nsplit=int(nsplit)

currd=os.getcwd()
os.chdir(ind)
#gather and unzip trajectory files, if zipped
iszip='false'
tras=chomp(os.popen('ls -1 rep*.tra* 2>/dev/null').readlines())
if not tras:
    sys.stderr.write('ERROR: no rep*.tra files in '+ind)
    sys.exit(1)
if re.compile('tbz2').search(tras[0]):
    iszip='true'
    for tra in tras:
        pastry('bunzip2 tra')

prefix='prefix'+`randint(0,1000)`

#create nsplit tra.in files
ntra=len(tras) #number of rep*.tra files
for n in range(0,nsplit):
    id="%02d"%(n)
    pastry('echo "'+`ntra`+'" > tra'+id+'.in')

for tra in tras:
    #count number of snapshots per trajectory
    line=os.popen('head -1 '+tra).readline()
    nres=int(re.compile('\s+(\d+)').match(line).group(1)) #find number of residues
    line=os.popen('grep " '+`nres`+' " '+tra+' | wc ').readline()
    nsnap=int(re.compile('\s+(\d+)').match(line).group(1))   #find number of snapshots
    nlines=(nres+1)*(nsnap/nsplit) #number of lines to output per piece file
    #create temporary piece files
    pastry('split -d -l '+`nlines`+' '+tra+' '+prefix)
    pieces=chomp(os.popen('ls -1 '+prefix+'*').readlines())
    n=0
    prefix2=tra[:-4]
    #rename each piece file and create correspondin tra.in file
    for piece in pieces:
        id="%02d"%(n)
        pastry('mv '+piece+' '+prefix2+'_'+id+'.tra') #prexfix09 --> rep1_09.tra
        pastry('echo "'+prefix2+'_'+id+'.tra" >> tra'+id+'.in')
        n=n+1

#zip trajectories if they were previously unzipped.
if iszip=='true':
    for tra in tras: pastry('bzip2 '+tra)
os.chdir(currd) #go back to initial directory
sys.exit(0)
