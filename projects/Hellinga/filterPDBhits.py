#!/usr/bin/python
'''
read each of the xxxxx.inpdb.dat files, and create xxxxx.inpdb200.dat where
we filter out all hits of xxxxx to pdb\'s with more than 200 residues
'''
import sys,os
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb,extractHeader

ih=inpHand('Usage: filterPDBhits.py [options]\nExample:./filterPDBhits.py -a /gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/out_filtered -b 300',
           ' -a _AR_rootd root directory with hits to pdb',
           ' -b __filen0 globbed file name (def: xxxxx.inpdb.dat)',
           ' -c __headersf list of headers (def: root directory/success.list)',
           ' -d __maxPDBl maximum pdb length (def:200)',
           ' -e __filterf file with list of hits to keep (flag in place of -d, def:none)'
           )
ih.parse(locals(),sys.argv)

if not filen0: filen0='xxxxx.inpdb.dat'
if not headersf:headersf=rootd+'/success.list'
filterl=[]
if not filterf:
    if not maxPDBl:maxPDBl=200
    else: maxPDBl=int(maxPDBl)
else:
    filterl=chomp(open(filterf,'r').readlines())
successes=chomp(open(headersf,'r').readlines())
#successes=chomp(open(rootd+'/toy.list','r').readlines())
remnant=len(successes)
for header in successes:
    print 'remaining=%4d '%(remnant)+header
    remnant-=1
    filen=deglobb(filen0,header)
    filenn=rootd+'/'+header+'/'+filen
    if not os.path.exists(filenn): continue
    inf=open(filenn,'r') #;Bye(inf)
    buf=inf.readline()
    line=inf.readline()
    if len(line.split())<4:
        print 'ERROR:'+header
        continue
    while line:
        items=line.split()
        hit=items[2]
        if filterf:
            if hit in filterl: buf+=line
        else:
            M=int(items[3])
            if M<maxPDBl: buf+=line
        line=inf.readline()

    if filterf:
        outf=open(rootd+'/'+header+'/'+filen+'.filtered.dat','w').write(buf)
    else:
        outf=open(rootd+'/'+header+'/'+filen+`maxPDBl`+'.dat','w').write(buf)

sys.exit(0)
