#!/usr/bin/python

import sys
from spicker.spickerYangResultsManager import spickOut
from TMalign.findBestPDBhitByTMalign import printHeader
from inputArgs.inputArgs import inpHand

inpHand('Usage: findBestPDBHitForYangSpickerResults.py [options]\nFind the pdb structure with maximum TMalign to the top model.',
        ' -a _A_dir directory with output from spicker (def:current dir)',
        ' -c __rebuilt suffix to add to the combo files names, if using the refined files (def:none)',
        ' -d _A_pdbent x-globbed pdb structure (def: CAdir/xxxxx.pdb)',
        ' -e _A_pdblistf list of pdb headers (def:a list generated from CAdir)',
        ' -f _A_nat (x-globbed) native file',
        ' -g _A_outf output file (def:STDOUT)',
        ' -i __header header'
        ).parse(locals(),sys.argv)

if not outf: out=sys.stdout
else:
    out=open(outf,'w')
    out.write(printHeader()+'\n') #some info about the meaning of the output
if not dir: dir='.'
s=spickOut(dir=dir,nat=nat,rebuilt=rebuilt)
if s.readError:
    sys.stderr.write(s.readError+'\n')
    sys.exit(1)
#We can pass empty pdbent and pdblistf, since they will be initialized
line=s.findBestPDBbyTMalign(pdbent=pdbent,pdblistf=pdblistf)
if header: line=header+' '+line
out.write(line+'\n')

sys.exit(0)
