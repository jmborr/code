#!/usr/bin/python

import os,sys,re
from jobs.job import pastry
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp
from jobs.job import job

ih=inpHand('Usage: dukka_pulchra.py',
           ' -a _R_header header',
           ' -b __infix infix',
           ' -c __doclosc do we refine closc instead of combo? (def:"no")'           
           ).parse(locals(),sys.argv)

if not doclosc or doclosc in ['n','N','no','NO','No']: type='combo'

pastry('/bin/cp /gpfs1/u/jose/bin/dukka.x /gpfs1/u/jose/bin/pulchra .')
n=0
for i in (1,2,3,4,5,10,25,50,75,99):
    closcf='closc'+`i`+'.pdb'
    if os.path.exists(closcf):
        n+=1
        if type=='combo':
            pastry('pulchra -b -s '+closcf) #only keep the CA's
            combof='combo'+`i`+'.pdb'
            pastry('./dukka.x '+closcf+'.rebuilt '+combof)
            pastry('./pulchra best.pdb')
            pastry('/bin/mv best.pdb.rebuilt '+header+'.'+infix+'.combo'+'%02d'%(i)+'.pdb.rebuilt')
        else:
            pastry('pulchra '+closcf)
            pastry('/bin/mv '+closcf+'.rebuilt '+header+'.'+infix+'.closc'+'%02d'%(i)+'.pdb.rebuilt')
            
            
if n==0: sys.exit(1) #there was no closc*.pdb files

sys.exit(0)
