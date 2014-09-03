#!/usr/bin/python
'''

It turns out that assigning secondary structure based solely on global
distributions for helixes and loops does reproduce such distributions,
but do not reproduce the distributions of number of secondary elements
in a protein of given sequence length. This happens because of
correlations in the lenght of secondary structures, the hairpin being
the epitome of such type of correlation. Thus, we find that assigning
secondary structure based solely on global distributions for helixes
and loops results in an excess of elements, that is, the placed
elements tend to be shorter than usual. This is true for placed
helixes in helical proteins, placed strands in beta proteins, and
placed strands in alpha/beta proteins. Apparently the distribution of
number of placed helixes in alpha/beta proteins agrees with
experimental results.

Create length distributions for helix, strands, and loops, each
derived from a subset of proteins with a given number of elements of
secondary structure.

First, given a target length Ntarget, we derive the distribution of
number of helixes distr[\'H\'] in a compact protein of sequence length
in the [Ntarget-50,Ntarget+50] range. Then, we determine the maximum
of the distribution, distrMax, and from then on we will only consider
producing proteins with a number of helixes nH such that
distr[\'H\'][nh]>0.2*distrMax. The reason is that less a number of
helixes less than 0.2*distrMax will have few cases, and futhermore,
these proteins will tipycally have lengths that deviate most from
Ntarget.

Given a number of helixes, we will scan all proteins with such number
of helixes and a sequence length in the [Ntarget-50,Ntarget+50]
range. The goal is to make a distributions of helix lengths and
intervening loops for such a set of proteins. However, there are not
many proteins that comply with the previous requierements, thus the
only thing we can do is to measure the average length and standard
deviation, then approximate the distribution of lengths to a normal
distribution of said average length and stardard deviation.

We will do the same for beta proteins. For alpha/beta proteins, we
will derive length distributions for a given number of placed strands,
while we will be using the global length distribution of helixes in
strands, since we observed that placement of helixes according to this
distribution does reproduce the distribution of number of placed
helixes. This is great, because we wouldn\'t have enough number of
proteins, if we had to derive length distributions based both on
number of placed strands and helixes.

We can not do the same for alpha/beta proteins because our length
distributions would be two-parameter depending, and we won\'t have
enough population. However, it turns out this is not neccessary
because alpha/beta protein don\'t have correlations between ther
elements of secondary structure, as helix and beta proteins do show.

The output is a file where each line has the number of helixes (or
strands) and the average length and stardard deviation to construct
the normal distribution of lengths. In addition, we also write the
Skewness and excess Kurtosis. There are also lines for strands and for
interior loops in alpha and beta proteins.

'''
import sys,os,re
from utilities.small_utilities import chomp,Bye
from stats.statistics import average,standard_dev,excess_kurtosis,skewness
from inputArgs.inputArgs import inpHand
from random import random,randint

ih=inpHand('Usage: genDistrByL.py',
           ' -a __Ntarget target sequence length (def: 150)',
           ' -b _A_outf output file (def:STDOUT)',
           ' -c __Nseqs number of sequences to output per alpha and beta class (def: 400)',
           ' -d _A_outd output directory (def: ./) will append subdirs a and b',
           ' -e _A_f secondary structure assignment file (def:scratchdir/foldsOrigin/ss.txt'
           ).parse(locals(),sys.argv)

if not Ntarget: Ntarget=150
else: Ntarget=int(Ntarget)
if not outf: outpt=sys.stdout
else: outpt=open(outf,'w')
if not Nseqs: Nseqs=400
else: Nseqs=int(Nseqs)
if not outd: outd=os.getcwd()
if not f: f='/gpfs1/scratch/jose/foldsOrigin/ss.txt' #fasta-like entries for sec. struct.
delta=25
maxLength=40
 
cad='/gpfs1/archive/adrian/dat/pdb_may06.apr06.aux/CA'
rgx='/gpfs1/active/jose/code/cpp/pdb_utils/Rg.x' #executable to calculate radius of gyration
Nx={'H':0.0,'E':0.0,'-':0.0}
nXinY=['nHinH','nEinE','nHinHE','nEinHE']
alph=['H','E','-']  #H:helix E:strand -:(interior)loop
alph2dir={'H':'a','E':'b'}
ssx={'H':2,'E':4,'C':1}  #numeric code for secondary structures
#below are different chunks: helix, strand, N-termini-loop, C-termini-loop, loop between helixes, loop between strands, , loop-preceeded-of-helix-and-followed-by-strand, loop-preceeded-of-strand-and-followed-by-helix
alph2=['(H+)','(E+)','^(-+)','(-+)$','H(-+)H','E(-+)E','H(-+)E','E(-+)H']
px={}
qx={}

N=264 #maximum assumed length for a contiguous stretch
for x in alph2:
    px[x]=re.compile(x) #regular expression to detect contigous stretches
for x in Nx.keys():
    qx[x]=re.compile(x)
    
pf=open(f,'r')
l=pf.readline()

#contain type (H,E), number of secondary struct. elem. (nSS), and a
#list with the lengths of each elemnt (lengths)
protein={}
nstr=0
while(l):
    if l[0]=='>':
        nstr+=1
        #print l,nstr
        header=l[1:6]
        
        #store secondary structure assignment in ll
        ll=''
        l=chomp(pf.readline())
        if l=='': l=chomp(pf.readline())#there's a (stupid) blank line in between entries
        while not l[0]=='>':
            ll=ll+l #glue all lines of the fasta entry
            l=chomp(pf.readline())
            if l=='': l=chomp(pf.readline())
            if not l: break

        #filter according to protein length and radius of gyration
        caf=cad+'/'+header+'.pdb' #CA file
        L=int(os.popen('grep " CA " '+caf+' | wc -l').readline().strip())
        if L<Ntarget-delta or L>Ntarget+delta: continue #protein length not in target length
        rg,rg0=os.popen(rgx+' -a '+caf+' -b ').readline().split()
        rg=float(rg) ; rg0=float(rg0)
        if rg > 1.25 * rg0: continue #protein is not compact

        #Find secondary structure content
        LL=len(ll)
        if LL>0:
            for x in Nx.keys(): Nx[x]=len(''.join(qx[x].findall(ll))) #number of H's, E's, -'s
            for x in Nx.keys(): Nx[x]=(1.0*Nx[x])/LL
            content='%4.2lf %4.2lf %4.2lf'%(Nx['H'],Nx['E'],1-Nx['H']-Nx['E'])
            if   Nx['E']<0.1 and Nx['H']>0.19: type='H'
            elif Nx['H']<0.1 and Nx['E']>0.19: type='E'
            elif Nx['H']>0.1 and Nx['E']>0.1:  type='HB'
            #print '%2s'%(type)+' '+content+' %3d %4.1lf %4.1lf'%(L,rg,rg0)+' '+header
        else: continue #some error forming ll

        if type=='HB':continue
        
        #find number of helices or strands per sequence, and their
        #lengths. Also, renormalize the length of each fragment with
        #Ntarget
        #ratio=float(LL)/Ntarget
        ratio=1
        key='('+type+'+)'
        iterator=px[key].finditer(ll)
        nSS=0
        lengths=[]
        for match in iterator:
            length=int( round(ratio*len(match.group(1))) )
            if length<maxLength: lengths.append(length)  #do not include long elements
            nSS+=1

        #find difference of consecutive lengths between strands
        diffs=[]
        for i in range(len(lengths)-1): diffs.append(lengths[i+1]-lengths[i])
            
        #find length of interior loops
        key=type+'(-+)'+type
        iterator=px[key].finditer(ll)
        lengthsL=[]
        for match in iterator:
            length=len(match.group(1))
            if length<maxLength: lengthsL.append(length) #do not include long loops

        #store protein info
        protein[header]={'type':type,'nSS':nSS,'lengths':lengths,'lengthsL':lengthsL,'diff':diffs}
        #if type=='H' and nSS>8: print header
            
pf.close()
line=''
for type in ('H','E'):
    subdir=outd+'/'+alph2dir[type] #directory to output sequences
    os.system('/bin/mkdir -p '+subdir)
    iseq=0 #current index for output sequence
    line+=type+'\n'
    lengths2headers={} #list of headers with certain number of elements
    for i in range(100): lengths2headers[i]=[]
    
    #find histogram of number of elements
    distrnSS=[0]*100  #initialize histogram of number of elements
    for header in protein.keys():
        if protein[header]['type']==type: #the protein is of appropiate type
            nSS=protein[header]['nSS']
            distrnSS[nSS]+=1
            lengths2headers[nSS].append(header)

    #find maximum of histogram
    distrMax=0
    nSSMax=0
    for nSS in range(len(distrnSS)):
        if distrnSS[nSS] > distrMax:
            distrMax=distrnSS[nSS]
            nSSMax=nSS

    #set to zero any bin with less than 0.2*distrMax
    for nSS in range(len(distrnSS)):
        if distrnSS[nSS] < 0.2* distrMax: distrnSS[nSS]=0

    #normalize distrnSS and output:
    line+='distribution of number of elements\n'
    norm=0
    for nSS in range(len(distrnSS)): norm+=distrnSS[nSS]
    for nSS in range(len(distrnSS)):
        distrnSS[nSS]=float(distrnSS[nSS])/norm
        if distrnSS[nSS]>0: line+='%2d %5.3lf\n'%(nSS,distrnSS[nSS])

    #calculate the number of sequences we have to output per nSS bin:
    nseqPernSS=[]
    for p in distrnSS:nseqPernSS.append( int(p*Nseqs) )
    remainder=Nseqs
    for i in nseqPernSS: remainder -= i
    nseqPernSS[nSSMax]+=remainder #add remainder to most expected, so that we have exactly Nseqs
    
    #for each non-zero bin of distrnSS, find normal distribution
    #parameters for number of elements and interior loops
    totalHeaders=0
    line+='distributions of element, loop length, and length difference between consecutive elements for proteins with fixed number of elements\n'
    for nSS in range(len(distrnSS)):
        print 'nSS=',nSS
        if distrnSS[nSS]>0:
            #lists of all elements and interior loops for proteins of a given nSS
            lengths={'SS':[],'diff':[],'loop':[]}
            for header in lengths2headers[nSS]:
                lengths['SS']=lengths['SS']+protein[header]['lengths']
                lengths['loop']=lengths['loop']+protein[header]['lengthsL']
                lengths['diff']=lengths['diff']+protein[header]['diff']
            #generage histograms of the previous lists. Note that for 'diff',
            #x=lengths['diff'][i] can be negative so that
            #histLength[x]==histLength[2*maxLength+x]
            histLength={'SS':[0]*maxLength,'loop':[0]*maxLength,'diff':[0]*(2*maxLength)}
            for kind in ('SS','loop','diff'):
                for i in lengths[kind]: histLength[kind][i]+=1#histogram of loop length
            #output histogram and some of the distribution moments
            for kind in ('SS','loop'):
                hist=histLength[kind]
                line+='%2d\n'%(nSS)
                for i in range(maxLength): line+=' %1d'%(hist[i])
                line+='\n'
                avL=average(lengths[kind])
                std=standard_dev(lengths[kind])
                kurt=excess_kurtosis(lengths[kind])
                skew=skewness(lengths[kind])
                #line+='%2d %5.2lf %5.2lf %5.2lf %5.2lf %4d\n'%(nSS,avL,std,kurt,skew,len(lengthList))            
            hist=histLength['diff']
            line+='%2d\n'%(nSS)
            for i in range(2*maxLength): line+=' %1d'%(hist[i])
            line+='\n'
            
outpt.write(line)#output distributions

sys.exit(0)
