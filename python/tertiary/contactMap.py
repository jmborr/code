#!/usr/bin/python

import sys,os
from random import randint
from spicker.spickerResultsManager import TMscore
from inputArgs.inputArgs import inpHand
from tertiary.disfit import disfit2 #square of the disfit distances
from utilities.small_utilities import junkName

def genCM(pdbf,co=0):
    if not os.path.exists(pdbf):
        sys.stderr.write('ERROR: '+pdbf+' does not exists\n')
        return []
    junk=junkName() #some junk file
    os.system('/gpfs1/active/jose/code/perl/tertiary/CA2SG.pl '+pdbf+' '+junk)

    #read junk file and store codes and coordinates of sidechains
    aas=[] ; x=[] ; y=[] ; z=[] ; cm=[] ; N=0
    for line in os.popen('grep -e " SG " '+junk+' ').readlines():
        aas.append(line[17:20])
        x.append(float(line[30:38]))
        y.append(float(line[38:46]))
        z.append(float(line[46:54]))
        N+=1

    for i in range(0,N-2-co):
        aai=aas[i] ; xi=x[i] ; yi=y[i] ; zi=z[i]
        for j in range(i+1+co,N):
            dx=xi-x[j] ; dy=yi-y[j] ; dz=zi-z[j]
            if dx*dx+dy*dy+dz*dz < 1.5*disfit2[aai][ aas[j] ]:
                cm.append( [i+1,j+1] )
    os.system('/bin/rm '+junk)
    return cm

if __name__=='__main__':   
    inpHand('Usage: contactMapDisfit.py [options]\nFor a given file, find the contact map according to the distance cut-off of disfit. Will input CA atoms, estimate the SG positions,\ncalculate the contact map based on DISFIT matrix',
            ' -a _AR_inPDBf model in pdb format',
            ' -b _A_outPDBf output file (def: STDOUT)',
            ' -c __cutoff we do not output short-range contacts (def: 3)',
            ).parse(locals(),sys.argv)

    if not outPDBf: out=sys.stdout
    else:
        try: out=open(outPDBf,'w')
        except ValueError:  print outf,'can\'t be opened'
    if not cutoff: cutoff=3
    else: cutoff=int(cutoff)
    cm=genCM(inPDBf,co=cutoff)
    for pair in cm:
        out.write( '%3d %3d\n'%(pair[0],pair[1]) )

    sys.exit(0)
