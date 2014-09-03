#!/usr/local/python
'''

'''
import os,sys,re,copy
from inputArgs.inputArgs import inpHand
from utilities.codedir import scratchdir,codedir
from utilities.small_utilities import Bye,chomp,junkName
from jobs.job import pastry

#formatdb needs input in "letters", not bytes. byte2letter is an
#estimate of the letter/byte ratio
byte2letter=4

#Parse the command line
ih=inpHand('Usage: filter.py',
           ' -a __volsize volume size, in MegaBytes (default=500)',
           ).parse(locals(),sys.argv)

#initialize defaults
if not volsize: volsize=500*byte2letter
else volsize=int(volsize)*byte2letter

#check we are not in a front node
if os.popen('hostname').readline()[0:4]=='fec0': Bye('do not run from a frontnode')

#clean and work in the database directory
nrdb=os.path.join(scratchdir,'db/nr')
pastry('/bin/rm '+nrdb+'/nr*') #clean old database
os.chdir(nrdb)

#bring database from local computer to database directory
pastry('/usr/bin/ssh fec01 "scpFromLocalComputer /tmp/nr.gz ."')

#mask low complexity, transmembrane, composition bias, and coiled-coil regions
pastry('/bin/gunzip nr.gz')
pfilt.exe=os.path.join(codedir,'c/pfilt/pfilt')
pastry(pfilt.exe+' nr > nr.filter ')
pastry('/bin/gzip nr')

#index the database
formatdb.exe=os.path.join(codedir,'/bin/blast/formatdb')
pastry(formatdb.exe+' -i nr.filter -o T -v '+volsize)

#create appropriate links nr->nr.filter and nr.pal file
pastry('/bin/ln -s nr.filter nr')
pastry('for l in `ls -1 nr.filter.??.p*`;do ll=`echo nr.${l:10:6}` && /bin/ln -s $l $ll;done')
pastry('cat nr.filter.pal | perl -ne \'s/\.filter//g; print;\' > nr.pal')

sys.exit(0)
