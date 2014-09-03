#!/usr/bin/python

import os,sys
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import inpHand

inpHand('Usage: findBestPDBHits.py',
        ' -a _AR_model pdb model to align to template library'
        ).parse(locals(),sys.argv)

mintm=0.40
tempd='/tmp/jose/pdbSubset'
list=chomp(open(tempd+'/list').readlines())

l=int(os.popen('grep " CA " '+model+'|wc -l').readline().strip())
pout=open('bestHits.list','w')
pout.write('#template L(template)/L(model) TM(template,model)\n')
for templ in list:
    templ2=tempd+'/'+templ
    cmd='TMalign '+templ2+' '+model+'|grep -a TM-score=|cut -d \'=\' -f 4|cut -d \',\' -f 1'
    tm=float(os.popen(cmd).readline().strip())
    if tm>mintm: #output the hit as a putative significant hit
        L=int(os.popen('grep " CA " '+templ2+'|wc -l').readline().strip())
        #The ration of template and model lengths is important, as we
        #big templates tend to have larger TM scores        
        pout.write('%10s %4.2f %5.2f\n'%(templ,(1.0*L)/l,tm))
pout.close()
sys.exit(0)
