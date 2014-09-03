#!/usr/bin/python
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.codedir import pdbnonhom
from utilities.small_utilities import chomp

"""
Look at seq.dat and estimate the number of helixes/strands, as well as
overall secondary structure content
Typical seq.dat line:
4   SER    1    7
5   GLU    2    9
6   TRP    4    1
Sec.Struc. codes are third column (helix:2,beta:4,coil:1)

"""

minNSS0={'H':5,'B':4}   #minimal lengths for secondary structure elements
maxNSS0={'H':24,'B':11} #maximal lengths for secondary structure elements
rs=[re.compile('HCH'),re.compile('BCB')]
            
def findNsec(header,seqd='',minnss=''):
    if not seqd: seqd=pdbnonhom+'/seq/'+header+'.SEQ'
    minNSS={}
    if not minnss:minNSS=minNSS0
    else: minNSS['H'],minNSS['B']=minnns.split()
    secHBC=''
    #find number of residues assigned as helix, beta, and coil
    lines=chomp(open(seqd,'r').readlines())    
    nrh=0;nrb=0;nrc=0
    for line in lines:
        n=int(line.split()[2])        
        if n==2:
            nrh+=1
            secHBC+='H'
        elif n==4:
            nrb+=1
            secHBC+='B'
        else:
            nrc+=1
            secHBC+='C'
    #print secHBC;sys.exit(1)
    ssct={'H':nrh,'B':nrb,'C':nrc} #secondary structure content    
    l=nrh+nrb+nrc #seqence length
    #calculate number of secondary structures
    #first, eliminate some isolated C: HCH->HHH, BCB->BBB
    secHBC=rs[0].sub('HHH',secHBC)  #;print secHBC
    secHBC=rs[1].sub('BBB',secHBC)  #;print secHBC
    #second, count helixes and strands
    nss={}
    for s in ['H','B']:
        nss[s]=0 #initialize counts
        minssl=minNSS[s] #minimal length of the secondary structure element
        maxl=maxNSS0[s]+minssl #signal to count as two elements of secondary-structure
        regex=s+'+'  #;print regex
        list=re.compile(regex).findall(secHBC) #list of helixes of list of strands
        for hb in list:
            if len(hb)<minssl: continue #too short, go to next element
            if len(hb)>maxl: nss[s]+=2 #count twice if too long
            else: nss[s]+=1
    #return results (length, sec.str.content, number sec.str)
    return l,ssct,nss,secHBC
    
inpHand('Usage: numSecElemSeqDat.py',
        ' -a __header five-letter pdb code',
        ' -b _A_hlist list of headers, mandatory if flag -a not provided',
        ' -c _A_seqd seq.dat file, in place of -a or -b',
        ' -d __minl min length of helix and beta sec.str.elem (def: "6 3")',
        ' -e __minS min sep between two consec sec str elem (def=2)',
        ' -f __secout output string with secondary structure also? (def="no")',
        ' -g _A_outf output file (def:STDOUT)'        
        ).parse(locals(),sys.argv)

if not header and not hlist and not seqd:
    sys.stderr.write('ERROR: either -a or -b or -c flag must be passed')
    sys.exit(1)
if not outf: out=sys.stdout
else:        out=open(outf,'w')
if header: headers=[header,]
if hlist:
    headers=open(hlist,"r").readlines()
    chomp(headers)
if not secout or secout in ['no','N0','n','N']: secout=''
#find number of secondary structures and sec.struct. composition
outputs={}
if seqd:
    if not header: header='xxxxx'
    outputs[header]=findNsec(header,seqd=seqd,minnss=minS)
else:
    for header in headers:
        outputs[header]=findNsec(header,minnss=minS)

#present results
out.write('#header seq-length fraction-unord-res (nrH-nrB)/(nrH+nrB) nsH nsB nsH+nsB\n')
for header in outputs.keys():
    output=outputs[header]
    l=output[0]    #sequence 
    ssct=output[1] #dictionary with secondary structure content
    nss=output[2]  #number of helixes and strands
    nsH=nss['H'] ; nsB=nss['B'] ; nsHB=nsH+nsB
    secHBC=output[3]
    fc=(1.*ssct['C'])/l #fraction of unordered residues
    x=ssct['H']+ssct['B']
    if x==0: hb=111. #no predicted secondary structure content!
    else: hb=(1.*ssct['H']-ssct['B'])/(ssct['H']+ssct['B'])
    out.write('%s %3d %4.2f %5.2f %3d %3d %3d\n'%(header,l,fc,hb,nsH,nsB,nsHB))
    if secout:
        out.write('---------1---------2---------3---------4---------5---------6\n')
        while secHBC:
            out.write(secHBC[0:60]+'\n')
            secHBC=secHBC[60:]


#exit with success
sys.exit(0)

