#!/usr/bin/python

import sys,os
from random import randint
from spicker.spickerResultsManager import TMscore
from inputArgs.inputArgs import inpHand
from tertiary.disfit import disfit2 #square of the disfit distances

def outNdangling(pdbf):
    junk='junk'+`randint(0,999)` #some junk file
    os.system('/gpfs1/active/jose/code/perl/tertiary/CA2SG.pl '+pdbf+' '+junk)

    #read junk file and store codes and coordinates of sidechains
    aas=[] ; x=[] ; y=[] ; z=[]
    N=0
    for line in os.popen('grep -e " SG " '+junk+' ').readlines():
        aas.append(line[17:20])
        x.append(float(line[30:38]))
        y.append(float(line[38:46]))
        z.append(float(line[46:54]))
        N+=1

    #find number of dangling amino acids in the N terminal of the model
    nterdanglingmodel=0
    contact=0 ; n=0
    for i in range(0,N-4):
        aai=aas[i] ; xi=x[i] ; yi=y[i] ; zi=z[i]
        for j in range(i+4,N):
            dx=xi-x[j] ; dy=yi-y[j] ; dz=zi-z[j]
            #stop if i makes a contact with j
            if dx*dx+dy*dy+dz*dz < disfit2[aai][ aas[j] ]:
                contact=1
                break
        if contact==1:
            nterdanglingmodel=n #;print 'nterdanglingmodel=',nterdanglingmodel
            break        
        n+=1
    #find number of dangling amino acids in the C terminal of the model
    cterdanglingmodel=0
    contact=0 ; n=0
    span=range(4,N) ; span.reverse()
    for i in span:
        aai=aas[i] ; xi=x[i] ; yi=y[i] ; zi=z[i]
        for j in range(0,i-4):
            aaj=aas[j] ; dx=xi-x[j] ; dy=yi-y[j] ; dz=zi-z[j]
            #stop if i makes a contact with j
            if dx*dx+dy*dy+dz*dz < disfit2[aai][aaj]:
                contact=1
                break
        if contact==1:
            cterdanglingmodel=n #;print 'cterdanglingmodel=',cterdanglingmodel
            break
        n+=1
    os.system('/bin/rm '+junk)
    return [nterdanglingmodel,cterdanglingmodel,nterdanglingmodel+cterdanglingmodel]

if __name__=='__main__':   
    inpHand('Usage: outNdangling.py [options]\nFor a given file, find the number of dangling residues\nin the N-terminal, the C-terminal, and the sum of both\ntermini. Will input CA atoms, estimate the SG positions,\ncalculate the contact map based on DISFIT matrix',
            ' -a _AR_pdbf model in pdb format',
            ' -b _A_outf output file (def: STDOUT)'
            ).parse(locals(),sys.argv)

    if not outf: out=sys.stdout
    else:
        try: out=open(outf,'w')
        except ValueError:  print outf,'does not exists'

    #print disfit2
    [nN,nC,nA]=outNdangling(pdbf)
    out.write('%3d %3d %3d\n'%(nN,nC,nA))

    sys.exit(0)
