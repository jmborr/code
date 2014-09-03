#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp
from seq.length import length_header

ih=inpHand('Usage: filterListByClass.py [options]\nExample:from pdbJuly05.list, select all entries below 500 amino acids long and of class alpha/beta',
           '  -a _AR_list file list of five letter codes (pdb id plus chain id)',
           '  -b _R_clasi class (_ALPHA_,_BETA_,_AB_,_APLUSB_)',
           '  -c __lco length cutoff (default: no length cutoff)',
           '  -d _A_outf output file (default: stdout)'           
	)
ih.parse(locals(),sys.argv)

#default values
out=sys.stdout
db='/gpfs1/active/jose/projects/db/dir.cla.scop.txt_1.69' #scop release
dict={'_ALPHA_':'a','_BETA_':'b','_AB_':'c','_APLUSB_':'d'} #scop annotation of prot.class
lco=10000000 #veeeeery long protein

#passed arguments
if outf: out=open(outf,'w')
if lco: loc=int(lco)
ldb=open(db,'r').readlines() #load the SCOP file
clas=dict[clasi]
pat0=re.compile('\s'+clas+'\.')

for header in chomp(open(list,'r').readlines()): #iter over entries in list
    header=header.lower() #transformation needed to conform to SCOP
    pat=re.compile(header)
    for l in ldb:          #iter over SCOP entries
        if pat.search(l):              #found entry in SCOP
            if pat0.search(l):
                L=length_header(header)
                if L<lco: #below length cut-off
                    out.write(header+' '+`L`+'\n')
sys.exit(0)
