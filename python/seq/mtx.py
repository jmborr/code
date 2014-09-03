#!/usr/bin/python
"""
 Documentation of the .mtx file format

Format of the .mtx file:

L     = Length of SEQ
SEQ   = Sequence
ka#-* = Karlin/Altschul parameters, block #.  There are three blocks,
        each containing four floating point numbers on separate lines.
pX-Y  = The position specific scores as integers.

The first element of this file format is [L].  This is the sequence
length.  The second line contains the sequence itself, in NCBI AA
notation.  After this, there are three KA blocks (four lines of
floating point numbers each), then the positional scores.

The positional scores are arranged in a grid.  Each line contains 26
elements, corresponding to the 26 elements in the NCBI AA encoding,
and there are L lines where L is the previously mentioned sequence
length.

NCBI AA encoding:
'-','A','B','C','D','E','F','G','H','I','K','L','M',
'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
'O', 'J'

Using the symbols mentioned above, it looks something like this:

<L>
<SEQ>
<ka1-1>
<ka1-2>
<ka1-3>
<ka1-4>
<ka2-1>
<ka2-2>
<ka2-3>
<ka2-4>
<ka3-1>
<ka3-2>
<ka3-3>
<ka3-4>
<p1-1> <p1-2> <p1-3> ... <p1-26>
<p2-1> <p2-2> <p2-3> ... <p2-26>
...
<pL-1> <pL-2> <pL-3> ... <pL-26>

One can find the explanation for the three blocks of KA-parameters in
makemat's source code, lines 188-190:

    putMatrixKbp(checkFile, compactSearch->kbp_gap_std[0], scaleScores, 1/scalingFactor);
    putMatrixKbp(checkFile, compactSearch->kbp_gap_psi[0], scaleScores, 1/scalingFactor);
    putMatrixKbp(checkFile, sbp->kbp_ideal, scaleScores, 1/scalingFactor);

Thus, the first KA block is the standard score, the second is for
PSI-Blast, and the third is the ideal score."""
import sys,re,os
from math import exp,sqrt
from seq.get_seq import get_seq
from inputArgs.inputArgs import inpHand,deglobb,addAbsPath,isglobbed
from jobs.job import pastry
from utilities.codedir import ncbidir
from utilities.small_utilities import abort,junkName,chomp
from seq.blastManager import blastRun
from seq.fastaManager import importFastaEntry

#following parameters derived from  Robinson and Robinson, PNAS 88 8880-8884 (1991)
ncbi_aa_encoding ='-ABCDEFGHIKLMNPQRSTVWXYZU*' #standar ncbi ordering of amino acids
isbiologicalamino='01011111111111111111101000'
isbiolaal=[ int(n) for n in list(isbiologicalamino) ]
lbd=0.3176 #lambda parameter of Robinson and Robinson paper
#average amino acid composition 
background_probs={'A':0.07805,'C':0.01924,'D':0.05364,'E':0.06295,'F':0.03856,'G':0.07377,'H':0.02199,'I':0.05142,'K':0.05744,'L':0.09020,'M':0.02243,'N':0.04487,'P':0.05202,'Q':0.04264,'R':0.05130,'S':0.07120,'T':0.05841,'V':0.06441,'W':0.01330,'Y':0.03216}
bprobs=[0.07805,0.01924,0.05364,0.06295,0.03856,0.07377,0.02199,0.05142,0.05744,0.09020,0.02243,0.04487,0.05202,0.04264,0.05130,0.07120,0.05841,0.06441,0.01330,0.03216]
defaultScalingFactor=1.0/100 #from makemat default settings

class mtx_score:
    
    def __init__(self,pointer,auxpointer=None):
        self.error=''
        self.scoress=[]
        self.scoressnormalized=[]
        #determine what is pointer (filename,file-descriptor,list)
        pin=None
        lines=[]
        if isinstance(pointer,str): #pointer if a file name
            try:
                pin=open(pointer,'r')
            except:
                self.error='ERROR mtx_profile::__init__ could not open '+pointer
                return
        if isinstance(pointer,file): pin=pointer #pointer is a file descriptor
        if pin:
            lines=pin.readlines()
            pin.close()
        else:
            if isinstance(pointer,list): lines=pointer #pointer is a list containing lines of mtx
        self.L=int(lines[0].strip())
        self.seq=lines[1].strip()
        #initialize set of position specific scores
        j=14 ; i=self.L
        while i:
            scores20=[]
            scores26=[int(x) for x in lines[j].split()]
            for k in range(26): #only retain the 20 biological amino acids
                if isbiolaal[k]: scores20.append( scores26[k] ) #squeeze out non-biological aa
            self.scoress.append( scores20 )
            j+=1 ; i-=1
        #if passed, read auxiliary file to determine scalingFactor
        self.scalingFactor=defaultScalingFactor #default initialization
        if auxpointer:
            if isinstance(auxpointer,str): #pointer if a file name
                try:
                    pin=open(auxpointer,'r')
                except:
                    self.error='ERROR mtx_profile::__init__ could not open '+auxpointer
                    return
            if isinstance(auxpointer,file): pin=auxpointer #pointer is a file descriptor
            self.scalingFactor=1/float( pin.readlines()[7].strip() ) #read scaling factor from aux
    
#       *************************************************

    def initNormalScores(self):
        for scores in self.scoress:
            scoresnormalized=[]
            for score in scores: scoresnormalized.append(self.scalingFactor*lbd*score)
            self.scoressnormalized.append( scoresnormalized )
        
#       ************************************************
    def initProfiles(self):
        """
        retrieve profile from position-specific scores

        score(position=i,residue=a)=ln( estimated_prob(i,a)/bprobs(a) ), according to
        Altschul et al, Nucl.Acids.Res 25, 3389-3402 (1997).
        The profile at position=i is the colection of estimated probabilities

        """
        if not self.scoressnormalized: self.initNormalScores() #normalize scores
        self.profiles=[]
        for scores in self.scoressnormalized: #run through all the sequence
            profile=[]
            for i in range(20): profile.append( bprobs[i]*exp(scores[i]) )
            #normalize summ of profile components exactly to one. It
            #should be one by construction, but roundoff errors from
            #mtx file prevent this.            
            t=0.0
            for x in profile: t+=x
            for i in range(20): profile[i]=profile[i]/t
            self.profiles.append(profile)
        
#   ************************************************
def mtx1_mtx2_Z_rank(mtx1,mtx2,alignment=[]):
    """
    calculate the similarity between two BLAST-makemat profiles with a Z-score and rank
    
    mtx1: .mtx file for first profile
    mtx2: .mtx file for second profile

    alignment: optional parameter, empty list by default which
    assumes both profiles have same length. If we pass this
    parameter, it has to be of the form [(i1,j1),(i2,j2),...]
    which indicates position i1 in mtx1 is aligned to position
    j1 in mtx2, and so on. The distance will be calculated only
    on the aligned residues.

    calculate profile-profile distance averaged over all pairs from
    aligment, and compare to all possible distances resulting from
    comparing one profile at a particular position of seq1 with a
    profile at a particular position fo seq2. Comparison is made with
    Z-score and with the percentile rank of distances bigger than the
    calculated average distance
      
    """
    ms1=mtx_score(mtx1) #construc mtx_score object from mtx1 file
    ms2=mtx_score(mtx2)
    if ms1.error or ms2.error: return -1
    for ms in (ms1,ms2): ms.initProfiles() #construct profiles from scores

    #create default alignment
    if not alignment:
        if ms1.L != ms2.L: return -1
        for i in range(ms1.L): alignment.append( (i,i) )
    nalg=len(alignment) #number of residues aligned
    
    #calculate average distance for the alignment
    dx=0.0
    for (i,j) in alignment:
        iprof=ms1.profiles[i]
        jprof=ms2.profiles[j]
        d2=0.0
        for k in range(20):
            dk=iprof[k]-jprof[k]
            d2+=dk*dk
        dx+=sqrt(d2)
    dx=dx/nalg

    #calculate average distances between all possible pair
    #combinations but such that amino acid types coincide. This
    #restriction due to the observation that applying mtx1_mtx2_Z_rank
    #without this restriction to the identical chunk of two
    #evolutionary unrelated sequences also give high Z-score. This
    #means that the typical distance between profile[i] for seq1 and
    #profile[j] for seq2 is much smaller if amino acid at position "i"
    #of seq1 is identical than amino acid at position "j" of seq2    
    av=0.0  ; sigma=0.0; n=0; nbetter=0; rank=0.0
    for i in range(ms1.L): #not neccessary to constraint to the aligment
        iaa=ms1.seq[i]
        iprof=ms1.profiles[i]
        for j in range(ms2.L):
            jaa=ms2.seq[j]
            if jaa==iaa:
#                if (i,j) in alignment: continue #do not include the aligned pairs
                jprof=ms2.profiles[j]
                d2=0.0
                for k in range(20):
                    dk=iprof[k]-jprof[k]
                    d2+=dk*dk
                sigma+=d2
                d=sqrt(d2)
                if d<dx: nbetter+=1               
                av+=d
                n+=1    
    av=av/n
    sigma=sqrt(sigma/n)-av
    rank=(100.0*(n-nbetter))/n

    #return Z-score rank, and number of distances better than the
    #average distance     
    return  (dx-av)/sigma,rank,nbetter

#   ************************************************

def gen_mtx(sequence,outf,header='protein',blastOuts={}):
    outf=addAbsPath(outf)
    #create temporary workding dir and work there
    currd=os.getcwd()  ;  tmpd=currd+'/'+junkName()  ;  pastry('/bin/mkdir -p '+tmpd)
    if blastOuts:
        os.system('/bin/cp '+blastOuts['chk']+' '+blastOuts['fasta']+' '+tmpd)
        os.chdir(tmpd)
    else:
        os.chdir(tmpd)
        blastOuts=blastRun([header,sequence]) #run blast with default options        
    pastry('echo '+blastOuts['chk']+' > '+header+'.pn')
    pastry('echo '+blastOuts['fasta']+' > '+header+'.sn')
    pastry(ncbidir+'/makemat -P '+header) #create the profile
    pastry('/bin/mv '+header+'.mtx '+outf)

    os.chdir(currd)  ;  pastry('/bin/rm -rf '+tmpd)  #remove temporary directory

def gen_mtx_2(blastarf,outf='',outd=''):
    """does same as gen_mtx, but we pass tarred PSI-BLAST results"""
    curd=os.getcwd()
    tmpd=junkName()
    pastry('/bin/mkdir '+tmpd+' && tar xf '+blastarf+' -C '+tmpd)
    os.chdir(tmpd)
    fastaf=os.popen('ls -1 *.fasta').readline().strip()
    header,sequence=importFastaEntry( open(fastaf,'r') )
    chkf=os.popen('ls -1 *.chk').readline().strip()
    pastry('echo '+chkf+' > '+header+'.pn')
    pastry('echo '+fastaf+' > '+header+'.sn')
    pastry(ncbidir+'/makemat -P '+header) #create the profile
    if isglobbed(outf): outf=deglobb(outf,header)
    if not outd: outd=curd
    if isglobbed(outd): outd=deglobb(outd,header)
    outf=os.path.join(outd,outf)
    pastry('/bin/mv '+header+'.mtx '+outf)
    os.chdir(curd)
    pastry('/bin/rm -r '+tmpd)

    
#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: mtx.py [options]',
               '  -a _A_fastaf (globbed) sequence file in fasta format',
               '  -b _A_seqf (globbed) raw sequence file',
               '  -c _A_pdbf (globbed) pdb file',
               '  -j _A_blastarf (globbed) tarred blast output',
               '  -d __pdbh pdb header (one of either -a -b -c -d -j)',
               '  -e __outf (globbed) output file name(default:header.mtx)',
               '  -f _A_listf list of headers, in case we pass some globbed argument',
               '  -o _A_outd (globbed) output directory(default:current dir)')
    ih.parse(locals(),sys.argv)

    blastarf2=''
    listH=[]
    if not listf: listH=['protein',]
    else: listH=chomp(open(listf,'r').readlines())
    for h in listH:        
        #produce the multiple sequence alignment, by passing the target sequence and a header
        #msa={'header':???,'target':???,'alignments':[???,???,...,???]}
        if  fastaf:
            fastaf2=deglobb(fastaf,h)
            header,sequence=importFastaEntry( open(fastaf2,'r') )
        elif seqf:
            seqf2=deglobb(seqf,h)
            sequence=get_seq(seqf2,inputType='seq')
        elif pdbf:
            pdbf2=deglobb(pdbf,h)
            sequence=get_seq(pdbf2,inputType='pdb')
        elif pdbh:            
            sequence=get_seq(pdbh,inputType='header')  ;  header=pdbh
        elif blastarf:
            blastarf2=deglobb(blastarf,h)            
        else: abort('ERROR: get_seq.py needs either one of -a -b -c -d -j options')
      
        if blastarf2:
            gen_mtx_2(blastarf2,outf=outf,outd=outd)
        else:
            if not header: header=h
            print 'header=',header
            #compose absolute path to output file
            if not outf: outf=header+'.mtx'
            else:
                outf=os.path.basename(outf) #remove absolute path
                outf2=deglobb(outf,header)
            outd2=deglobb(outd,header)
            if not outd2: outd2=os.getcwd()     #default is current directory
            if not os.path.exists(outd2): pastry('/bin/mkdir -p '+outd2)
            outf2=outd2+'/'+outf2                #compose absolute output filename            
            gen_mtx(sequence,outf2,header)

    sys.exit(0)
