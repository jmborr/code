#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from jobs.job import pastry
from random import randint
from utilities.small_utilities import chomp

inpHand('Usage: split_tra.py.py',
        ' -a _A_ind input directory (def=./)',
        ' -b __nsp number of spicker runs (def=10)'
        ).parse(locals(),sys.argv)

if not ind: ind=os.getcwd();
if not nsp: nsp=10
else: nsp=int(nsp)

currd=os.getcwd()
os.chdir(ind)

for n in range(0,nsp):
    id='%02d'%(n)
    tra='tra'+id+'.in'
    centr='centroids'+id+'.pdb'
    clos='closests'+id+'.pdb'
    summ='summary'+id+'.txt'
#    print '/gpfs1/active/jose/code/c/spicker/spicker1.0/spicker.x -l '+tra+' -s seq.dat -a CA -d '+centr+' -e '+clos+' -f '+summ+' -o . -b NOTAR'
    pastry('/gpfs1/active/jose/code/c/spicker/spicker1.0/spicker.x -l '+tra+' -s seq.dat -a CA -d '+centr+' -e '+clos+' -f '+summ+' -o . -b NOTAR')


os.chdir(currd) #go back to initial directory
sys.exit(0)
