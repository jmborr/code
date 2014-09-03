#!/usr/bin/python
import os,sys

from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb
from seq.fastaManager import importFastaEntry

inp=inpHand('Usage: secVar.py [options]\ncalculate percentage of every secondary type for every residue given a bunch of seq.dat files\n',
            ' -a _RA_seqdatlf list of seq.dat files',
            ' -b _A_outf output file (def:STDOUT)',
            ' -c _A_jefseqdatf jeff secondary structure assignments'
            )
inp.parse(locals(),sys.argv)

if not outf: outp=sys.stdout
else: outp=open(outf,'w')

n2t={'1':'C', '2':'A', '4':'E'} #secondary types in of seq.dat
jeff2t={'3':'C', '2':'A', '4':'E'}
    
    
seqdatl=chomp(open(seqdatlf,'r').readlines())
nseqdat=len(seqdatl)

nres=int(os.popen('wc -l '+seqdatl[0]).readline().split()[0]) #number of residues

#percentages for each residue.
pge=[]
for i in range(nres): pge.append({'A':0.0, 'E':0.0, 'C':0.0})

for seqdatf in seqdatl:
    lines=open(seqdatf,'r').readlines()
    for i in range(nres):
        ss= n2t[ lines[i].split()[2] ]
        pge[i][ss]+=1.0

if jefseqdatf:
    jefp=open(jefseqdatf,'r')
    jefp.readline() #fist line of no relevance
    
    jeffss=''.join( chomp(jefp.readlines()) ).replace(' ','')
    for i in range(nres):
        type=jeff2t[ jeffss[i] ]
        pge[i][type]-=nseqdat

buf='#res  %A    %E    %C\n'
for i in range(nres):
    for ss in ('A','E','C'): pge[i][ss]=pge[i][ss]/nseqdat
    buf+='%3d %5.2lf %5.2lf %5.2lf\n'%(i+1,pge[i]['A'],pge[i]['E'],pge[i]['C'])
    
outp.write(buf)
outp.close()

sys.exit(0)
