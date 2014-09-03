#!/usr/bin/python
'''

'''
import os,sys,re,copy
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName
from math import sqrt

pat={}

#############################################################
#        FUNCTION DEFINITIONS
#############################################################

def void_header():
    '''
    return a void chain.dat
    '''
    return '    0    hard'
  
#remove from template central residues of strands longer than 8 Angstroms
pat['remc']=re.compile('EEEEEEEEE+')
def removeCentralResFromStrands(sss,template):
    indexestoremove=[] #residue indexes of residues to remove
    it=pat['remc'].finditer(sss) #finds all strands higher than nine residues
    for match in it:
        (b,e)=match.span()
        c=1+(b+e-1)/2 #central residue
        for i in range(c-1,c+2): indexestoremove.append(i)
    template2=[]
    for line in template:
        index=int(line.split()[0]) #residue index
        if index not in indexestoremove: template2.append( line )
    return template2

#retrieve template for the case of all-alpha or all-beta    
def retrieveXYZ(ss,seq,xyz,remc=''):
    L=len(seq)
    template=[]
    for n in range(L):
        if seq[n]==ss:
            x=xyz[ss][n][0]
            y=xyz[ss][n][1]
            z=xyz[ss][n][2]        
            template.append( '%3d  %8.3lf  %8.3lf  %8.3lf'%(1+n,x,y,z) )
    #remove from template central residues of strands longer than 8 Angstroms
    if remc: template=removeCentralResFromStrands(''.join(seq),template)
    return template

#retrieve coordinates from displaced helical or strand chunk
def retrieveXYZH(prev,xyzH):
    template=[]
    for n in range( len(xyzH) ): #retrieve a helix chunk of length last-prev
        x=xyzH[n][0]
        y=xyzH[n][1]
        z=xyzH[n][2]
        template.append( '%3d  %8.3lf  %8.3lf  %8.3lf'%(1+prev+n,x,y,z) )
    return template

#generate coordinates for helical or strand chunk, and position of last placed residue
def translate(ss,length,xyz,u):
    v0=xyz[ss][0] #; print 'ss=',ss,'v0=',v0#origin of helix or strand
    t=[] #translation vector
    for i in range(3): t.append( u[i] - v0[i] ) #translation vector
    xyzH=[]
    for i in range(length):
        d=[0.0, 0.0, 0.0]
        vi=xyz[ss][i] #vector position of residue i of helix or strand
        for j in range(3): d[j]=vi[j]+t[j] #new vector position
        xyzH.append(d)
    return xyzH #helix translated
    
#We lay down every contiguous chunk of h, e, or c residues
def mixtureXYZ(seq,xyz,remc='',nostrand=''):
    L=len(seq)
    template=[]
    #u: position of last placed residue. Let's start at the beginning of the strand
    u=copy.deepcopy(xyz['E'][0])
    #print 'u=',u
    prev=0
    type=seq[prev] #type of current chunk
    last=0
    seq.append('u') #mark end of sequence
    while last<L:
        while seq[last]==type: last+=1
        length=last-prev #; print u,prev,last,`length`,type
        if type=='H':
            #translate helix to coordinates given by beta residue at position given by u
            xyzH=translate('H',length,xyz,u)
            u=copy.deepcopy(xyzH[-1])
            template+=retrieveXYZH(prev,xyzH)
        else: #use the coordinates of the long strand
            if not nostrand:
                xyzH=translate('E',length,xyz,u)
                u=copy.deepcopy(xyzH[-1])
                if type=='E': template+=retrieveXYZH(prev,xyzH) #(don't export coords for coil)
        prev=last
        type=seq[prev]
    seq.remove('u') #remove the mark from seq
    #remove from template central residues of strands longer than 8 Angstroms
    if remc: template=removeCentralResFromStrands(''.join(seq),template)
    return template
        
#############################################################
#        MAIN FUNCTION
#############################################################

#Parse the command line
ih=inpHand('Usage: genLongHelixChainDat.py',
           ' -a _RA_seqdat seq.dat file',
           ' -b _A_hepdb pdb containing the helix/strand superposition (def:./betahelix.pdb)',
           ' -c _A_outf output file (def: chain.dat)',
           ' -d _A_outpdb output pdb file with initial conformation (def: init.pdb)',
           ' -e __remc do not ouput three central residues for strands longer than 8 (def:no)',
           ' -f __nostrand do we output strands (def:yes)'
           ).parse(locals(),sys.argv)

#initialize defaults
if not hepdb: hepdb='./betahelix.pdb'
if not outf: outf='chain.dat'
if not outpdb: outpdb='init.pdb'
if remc and remc[0] in ('n','N'): remc=''
if nostrand and nostrand[0] in ('n','N'): nostrand=True
else: nostrand=False

#read coordinates of the beta (first) and helix (second)
xyz={'E':[],'H':[]}
inpt=open(hepdb,'r')
while 1: #first read beta
    line=inpt.readline()
    if line[0:4]=='ATOM' and line.find(' CA ')>0:
        x=float( line[30:38] )
        y=float( line[38:46] )
        z=float( line[46:54] )
        xyz['E'].append([x,y,z])
    elif line[0:3]=='TER':break
while 1: #second read helix
    line=inpt.readline()
    if line[0:4]=='ATOM' and line.find(' CA ')>0:
        x=float( line[30:38] )
        y=float( line[38:46] )
        z=float( line[46:54] )
        xyz['H'].append([x,y,z])
    if line[0:3]=='TER':break
#print xyz ; Bye('')

#read seq.dat
seq=[]
nss={'H':0,'E':0,'C':0}
n2ss={1:'C',2:'H',4:'E'}
for line in open(seqdat,'r').readlines():
    ss=n2ss[ int(line.split()[2]) ]
    seq.append(ss)
    nss[ss]+=1

#Three cases: all-alpha, all-helix, mixture
template=[]
if not nss['H']:     #protein has not helical content => It's a beta protein
    if nostrand:
        print 'WARNING: you said no strands but this is a beta protein!'
    else:
        template+=retrieveXYZ('E',seq,xyz,remc=remc)
elif not nss['E']:#protein has not beta content => It's a helical protein
    template+=retrieveXYZ('H',seq,xyz)
else:                #protein has both beta and helical content
    template+=mixtureXYZ(seq,xyz,remc=remc,nostrand=nostrand)

if nostrand and not nss['H']:
    template=['    0    hard',]
else:
    template=['    1    hard','%4d'%(len(template))]+template

#Output template and initial configuration
open(outf,'w').write('\n'.join(template))
out=open(outpdb,'w')
for line in template[2:]:
    (i,x,y,z)=line.split()
    i=int(i) ; x=float(x) ; y=float(y) ; z=float(z)
    out.write('ATOM   %4d  CA  VAL  %4d    %8.3lf%8.3lf%8.3lf\n'%(i,i,x,y,z))
out.write('END\n')

#Say bye
sys.exit(0)
