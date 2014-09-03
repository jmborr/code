#!/usr/bin/python

import sys,os
from inputArgs.inputArgs import inpHand
from utilities.codedir import yangpdb
from utilities.small_utilities import chomp
from tmalign_Yang_f77 import tmalign_yang

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: mkpair.py [options]',
               '  -a __templh template header',
               '  -b _A_templistf template list of headers (in place of -a  flag)',
               '  -c _A_templf0 template pdb file (in place of -b or -c)',
               '  -d __targeth target header',
               '  -e _A_targetlistf target list of headers',
               '  -f _A_targetf0 target pdb file (in place of -d or -e)',               
               '  -g _A_pdbDir directory where to find pdb files for headers (def:yangpdb)')
    ih.parse(locals(),sys.argv)

##parse the input
#prepare directory pdb's for templates
if templf0: pdbDirTempl='.'
else:
    if not pdbDir: pdbDirTempl=yangpdb
    elif pbdDir: pdbDirTempl=pdbDir
#prepare list of headers for templates
if templf0: templist=[templf0,]
else:
    if templh: templist=[templh+'.pdb',]
    else: templist=chomp(open(templistf,'r').readlines(),replace='.pdb')
        
#prepare directory pdb's for templates
if targetf0: pdbDirTarget='.'
else:
    if not pdbDir: pdbDirTarget=yangpdb
    elif pdbDir: pdbDirTarget=pdbDir

#prepare list of headers for templates
if targetf0: targetlist=[targetf0,]
else:
    if targeth: targetlist=[targeth+'.pdb',]
    else: targetlist=chomp(open(targetlistf,'r').readlines(),replace='.pdb')

#calculate TM-score of each combination of target and template
scoreMax=0.0
for targeth in targetlist:
    targetf=pdbDirTarget+'/'+targeth     #;print targetf
    for templh in templist:
        templf=pdbDirTempl+'/'+templh    #;print templf
        score=tmalign_yang(templf,targetf)
        if score>scoreMax:
            scoreMax=score
            besttempl=templh
    print 'target=',targeth,'template=',besttempl,'TM-score=',scoreMax

sys.exit(0)
