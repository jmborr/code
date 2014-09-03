#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp,Bye,junkName
from utilities.codedir import codedir
from math import sqrt

def distFromStrings(v1,v2):
    x1=float(v1[0:8]) ; y1=float(v1[8:16]) ; z1=float(v1[16:24])
    x2=float(v2[0:8]) ; y2=float(v2[8:16]) ; z2=float(v2[16:24])
    return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) )

#predicted extended and real end-to-end distance of the passed chunk
def end2endDists(chunk,type):
    l=len(chunk)
    #calculate predicted end-to-end distance of a chunk in the extended strand or extended helical
    if type=='2': pe2ed=l*1.5224 #helix
    elif type=='4': pe2ed=l*3.3652 #strand
    #calculate the real end-to-end distance of the chunk
    e2ed=distFromStrings(chunk[0][30:55],chunk[-1][30:55])    
    return [e2ed,pe2ed]

ih=inpHand('Usage: consistencyOfStructures.py [options]',
           ' -a _AR_pdbf pdbf file of the structures to check',
           ' -b _AR_seqdatf seq.dat-like file',
           ' -c _A_canonf pdb file containing the canonical strand and helix (def:codedir/projects/Hellinga/betahelix.pdb',
           ' -d _A_outf outpuf file (def: STDOUT)',
           ' -e __extreme remove flanking residues (def: 0)'
           )
ih.parse(locals(),sys.argv)

if not os.path.exists('./rmsd.x'): Bye('ERROR: please copy rmsd.x to working directory')

if not outf:outp=sys.stdout
else: outp=open(outf,'w')
if not canonf: canonf=codedir+'/projects/Hellinga/betahelix.pdb'
if not extreme:extreme=0
else: extreme=int(extreme)

#read CA's coordinates for pdb file
pdbl=os.popen('grep "ATOM" '+pdbf+'|grep " CA "').readlines() #;Bye(`len(pdbl)`)
basename=os.path.basename(pdbf) #;Bye(''.join(pdbl))
#read secondary structure assignment from seq.dat file
ss=''
for line in open(seqdatf,'r').readlines(): ss+=line.split()[2] #;Bye(ss)

#read CA's coordinates for canonical strand and helix. We know that
#first come the strand, then the helix. We also know that both are 150
#residues long
strandl=os.popen('grep "ATOM" '+canonf+'|grep " CA "').readlines()
helixl=strandl[150:]   #;Bye(''.join(helixl)) #;Bye(`len(helixl)`)
strandl=strandl[0:150] #;Bye(''.join(strandl)) #;Bye(`len(strandl)`)
canond={'2':helixl,'4':strandl}

#two temporary files
chunk1=junkName()
chunk2=junkName()

#iterate over all secondary structure chunks
buf=''
buf+=basename+'\n'
pat=re.compile('2+|4+')
it=pat.finditer(ss) #;Bye(ss)
for match in it:
    #print match.group()
    (b,e)=match.span() #;Bye(`b`+' '+`e`) #limits to beginning and ending of helix
    if e-b-2*extreme>0:
        l=e-b-2*extreme #length of the chunk
        type=ss[b] #'2' or '4'
        #print ''.join(pdbl[b+extreme:e-extreme])
        open(chunk1,'w').writelines( pdbl[b+extreme:e-extreme] )
        [e2ed,pe2ed]=end2endDists(pdbl[b+extreme:e-extreme],type) #end-to-end distances
        open(chunk2,'w').writelines( canond[type][0:e-b-2*extreme] )
        #print ''.join(canond[type][0:e-b-2*extreme])
        rmsd=float(os.popen('./rmsd.x '+chunk1+' '+chunk2).readline()) #;Bye('rmsd='+`rmsd`)
        buf+='%3d %3d %3d %5.2lf %4.1lf %4.1lf %s\n'%(l, 1+b+extreme,e-extreme,rmsd,e2ed,pe2ed,ss[b+extreme:e-extreme])
    
outp.write(buf)
os.system('/bin/rm '+chunk1+' '+chunk2)
sys.exit(0)
