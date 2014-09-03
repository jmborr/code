#!/usr/bin/python

import sys,os
from utilities.small_utilities import chomp
from spicker.spickerYangResultsManager import spickOut,TMscore

#return secondary structure content
def sec_str_content():
    content={'H':0.0,'E':0.0,'C':0.0} #number of residues assigned helix(H), strand(E) or coil(C)
    for line in open('seq.dat','r').readlines():
        state=line.split()[2]
        if state=='1': content['C']+=1.0
        elif state=='2': content['H']+=1.0
        else: content['E']+=1.0
    total=content['C']+content['H']+content['E']
    for key in content: content[key]=content[key]/total
    return content

listf='/gpfs1/scratch/jose/newShortProteins/benchmark/success.list'
#listf='/gpfs1/scratch/jose/newShortProteins/benchmark/toy.list'
myroot='/gpfs1/scratch/jose/newShortProteins/benchmark'
pdbdir='/gpfs1/active/adrian/short_prot/pdb'

print '#./gatherDataForAdrian.py > gatherDataForAdrian.dat'
print '#pdbID, model-number, TM to native, RMSD to native, average cluster RMSD, size/all, average TM to the other models'
        
headers=chomp(open(listf,'r').readlines()) #read listf
for header in headers:
    nat=pdbdir+'/'+header+'.pdb' #native file
    mydir=myroot+'/'+header[1]+'/'+header #; print mydir
    s=spickOut(dir=mydir,nat=nat)
    if s.readError: continue #error reading my output
    s.initTMtoNat()
    s.initRmsdsToNat()
    for id in s.ids:
        print '%s %2d %5.3lf %5.2lf %6.2lf %5.3lf %5.3lf'%(header,id,s.tms[id],s.rmsds[id],s.Rc_in[id],s.numberdens[id],s.avTMid(id=id) )
