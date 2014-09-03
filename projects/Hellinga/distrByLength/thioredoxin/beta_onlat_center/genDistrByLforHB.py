#!/usr/bin/python
'''

It turns out that assigning secondary structure based solely on global
distributions for helixes and strands of alpha/beta proteins does
reproduce such distributions, but does not reproduce the distribution
of number of strands in a protein of given sequence length. This
happens because of correlations in the lenght of strands, the hairpin
being the epitome of such type of correlation. Thus, we find that
assigning secondary structure based solely on global distributions for
helixes and loops results in an excess of strands, that is, the placed
strands tend to be shorter than usual. Apparently the distribution of
number of placed helixes in alpha/beta proteins agrees with
experimental results.

Create overall length distributions for helixes and loops, and
particular length distributions of strands derived from the subset of
proteins with a given number of strands.

First, given a target length Ntarget, we derive the distribution of
number of strands distr[\'E\'] in a compact protein of sequence length
in the [Ntarget-50,Ntarget+50] range. Then, we determine the maximum
of the distribution, distrMax, and from then on we will only consider
producing proteins with a number of strands nE such that
distr[\'H\'][nh]>0.2*distrMax. The reason is that less a number of
strands less than 0.2*distrMax will have few cases, and futhermore,
these proteins will tipycally have lengths that deviate most from
Ntarget.

Given a number of strands, we will scan all proteins with such number
of strands and a sequence length in the [Ntarget-50,Ntarget+50]
range. The goal is to make the following distributions:
* distributions of strand lengths
* distribution of helix lengths (remember we\'re dealing with alpha/beta
  proteins
* distribution of intervening loops

In addition, we will compute the ratio of helix to strands that we
placed based on the condition: L_H * N_H = L_E * N_E, that is the
average number of residues in helical state should equal that of
residues in strand state (this is true for experimental compact
alpha/beta proteins)

'''
import sys,os,re
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import inpHand

ih=inpHand('Usage: genDistrByL.py',
           ' -a __Ntarget target sequence length (def: 150)',
           ' -b _A_outf output file (def:STDOUT)',
           ' -c __Nseqs number of sequences to output per alpha and beta class (def: 400)'
           ).parse(locals(),sys.argv)

#initialize default arguments
if not Ntarget: Ntarget=150
else: Ntarget=int(Ntarget)
if not outf: outpt=sys.stdout
else: outpt=open(outf,'w')
if not Nseqs: Nseqs=400
else: Nseqs=int(Nseqs)

#Some global variables
line=''
delta=25
maxLength=40
f='/gpfs1/scratch/jose/foldsOrigin/ss.txt' #fasta-like entries for secondary struct. of each header
cad='/gpfs1/archive/adrian/dat/pdb_may06.apr06.aux/CA'
#cad='/tmp/jose/CA' #find this directory in cng0005 node
rgx='/gpfs1/active/jose/code/cpp/pdb_utils/Rg.x' #executable to calculate radius of gyration
Nx={'H':0.0,'E':0.0,'-':0.0}
alph2dir={'H':'a','E':'b'}
#below are different chunks: helix, strand, N-termini-loop, C-termini-loop, loop between helixes, loop between strands, , loop-preceeded-of-helix-and-followed-by-strand, loop-preceeded-of-strand-and-followed-by-helix
alph2=['(H+)','(E+)','^(-+)','(-+)$','H(-+)H','E(-+)E','H(-+)E','E(-+)H','[HE](-+)[HE]','-+(E+)','-+H+-+(E+)','-+H+-+(H+)','([EH]+)']
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
nstr=0 ;
while(l):
    if l[0]=='>':
        nstr+=1
        header=l[1:6] #;print l,nstr
        
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
            #print type+' '+content+' %3d %4.1lf %4.1lf'%(L,rg,rg0)+' '+header
        else: continue #some error forming ll

        if type=='H' or type=='E':continue
        
        #find number of strands per sequence, and lengths for strands
        #and helixes. Also lengths for difference between consecutive
        #strands and between strands separated by one helix        
        lengths={'E':[],'H':[],'loop':[],'E-E':[],'E-H-E':[]}

        #find length of strands and number of strands
        iterator=px['(E+)'].finditer(ll)
        nSS=0
        for match in iterator:
            length=len(match.group(1))
            if length<maxLength: lengths['E'].append(length)
            nSS+=1

        #find length of strands
        iterator=px['(H+)'].finditer(ll)
        for match in iterator:
            length=len(match.group(1))
            if length<maxLength: lengths['H'].append(length)

        #find length of interior loops
        iterator=px['[HE](-+)[HE]'].finditer(ll)
        for match in iterator:
            length=len(match.group(1))
            if length<maxLength: lengths['loop'].append(length) #do not include long loops

        #find difference of consecutive lengths between
        #strands. finditer() does not report overlapping matches, thus
        #we have to loop over each strand and see if its next element
        #is a strand, or its next two elements are helix followed by
        #strand        
        iterator=px['(E+)'].finditer(ll) #iterate over each strand
        for match in iterator:
            end=match.span()[1]    #where in ll does the strand finish?
            el=end-match.span()[0] #strand length
            match2=px['-+(E+)'].match(ll[end:]) #does a strand follows?
            if match2:
                el2=len( match2.group(1) )    #lenght of consecutive strand
                lengths['E-E'].append(el-el2) #append difference in strand lengths
            match2=px['-+H+-+(E+)'].match(ll[end:]) #does a helix-strand follows?
            if match2:
                el2=len( match2.group(1) )      #lenght of consecutive strand
                lengths['E-H-E'].append(el-el2) #append difference in strand lengths

        #find sequence of secondary elements
        seqSS=''
        iterator=px['([EH]+)'].finditer(ll)
        for match in iterator:
            segment=match.group(1)
            seqSS+=segment[0]
            
        #store protein info
        protein[header]={'nSS':nSS,'lengths':lengths,'seqSS':seqSS}

pf.close()

ntype={'H':0,'E':0}
#Find probabilities of strands and helixes
for header in protein.keys():
    for ss in protein[header]['seqSS']:
        if ss=='E': ntype['E']+=1
        else: ntype['H']+=1

norm=ntype['H']+ntype['E']
ntype['H']=float(ntype['H'])/norm
ntype['E']=float(ntype['E'])/norm

#Find probability of beginning singlet
line+='beginning singlet\n'
singlets={'H':0,'E':0}
for header in protein.keys():
    if protein[header]['seqSS']:
        singlets[ protein[header]['seqSS'][0] ]+=1
norm=sum(singlets.values())
for singlet in singlets.keys():
    singlets[singlet]=float(singlets[singlet])/norm
    line+='%s %5.3lf\n'%(singlet,singlets[singlet])

#Find probability of beginning doublet
line+='beginning doublet\n'
doublets={'HH':0,'EH':0,'HE':0,'EE':0}
for header in protein.keys():
    if len( protein[header]['seqSS'] )>1:
        doublets[ protein[header]['seqSS'][0:2] ]+=1
norm=sum(doublets.values())
for doublet in doublets.keys():
    doublets[doublet]=float(doublets[doublet])/norm
    line+='%s %5.3lf\n'%(doublet,doublets[doublet])    
    
#Find probability of each triplet of secondary structure elements
line+='triplets\n'
triplets={'HHH':0,'HHE':0,'HEH':0,'HEE':0,'EHH':0,'EHE':0,'EEH':0,'EEE':0}
for header in protein.keys():
    seqSS=protein[header]['seqSS']
    for ss in range(0,len(seqSS)-2):
        triplet=seqSS[ss:ss+3]
        triplets[triplet]+=1
norm=sum(triplets.values())
for triplet in triplets.keys():
    triplets[triplet]=float(triplets[triplet])/norm
    line+= '%s %5.3lf\n'%(triplet,triplets[triplet])

#Find histogram of number of strands
distrnSS=[0]*100  #initialize histogram of number of strands
lengths2headers={} #list of headers with certain number of elements
for i in range(100): lengths2headers[i]=[]
for header in protein.keys():
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

#output normalized distribution
line+='HE\ndistribution of number of strands in alpha/beta proteins\n'
norm=sum(distrnSS)
for nSS in range(100):
    if distrnSS[nSS]:
        distrnSS[nSS]=float(distrnSS[nSS])/norm
        line+='%2d %5.3lf\n'%(nSS,distrnSS[nSS])

#for each non-zero bin of distrnSS, find length distributions
line+='distributions of elements and loop length for proteins with fixed number of elements\n'
for nSS in range(len(distrnSS)):
    if distrnSS[nSS]:
        #lists of all elements and interior loops for  a given nSS
        lengths={'H':[],'E':[],'loop':[],'E-E':[],'E-H-E':[]} #list containing all lengths
        for header in lengths2headers[nSS]:
            for x in ('H','E','loop','E-E','E-H-E'):
                lengths[x]+=protein[header]['lengths'][x] #append appropriate list
        #generage histograms of the previous lists. Some elements of
        #lengths['E-E'] and lengths['E-H-E'] are negative, which mean
        #that hist[ lengths['E-E'] ] will be filled from the left        
        histLength={'H':[0]*maxLength,'E':[0]*maxLength,'loop':[0]*maxLength,'E-E':[0]*(2*maxLength),'E-H-E':[0]*(2*maxLength)}
        for x in ('H','E','loop','E-E','E-H-E'):
            hist=histLength[x]
            for L in lengths[x]: hist[L]+=1
        #output histograms
        for x in ('H','E','loop'):
            hist=histLength[x]
            line+='%2d %s\n'%(nSS,x)
            for L in range(maxLength): line+=' %1d'%(hist[L])
            line+='\n'
        for x in ('E-E','E-H-E'):
            hist=histLength[x]
            line+='%2d %s\n'%(nSS,x)
            for L in range(2*maxLength): line+=' %1d'%(hist[L])
            line+='\n'

outpt.write(line)#output distributions

sys.exit(0)
