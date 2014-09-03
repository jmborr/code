#!/usr/bin/python
import sys,os,re
from utilities.small_utilities import chomp,Bye
from stats.statistics import average,standard_dev,excess_kurtosis,skewness
from random import random,randint

#insert remnant residues in the loops
def addToLoops(remnant,loopChunks):
    nloops=len(loopChunks)
    for i in range(remnant):
        n=randint(0,nloops-1)
        #print 'nloops=',nloops,' loopChunks[',n,'][1]=',loopChunks[n][1]
        loopChunks[n][1]+=1

infinity=1000
nall=0
maxLength=40
Nseqs=400
Ntarget=105
outd='.'
alph2dir={'H':'a','E':'b','HE':'c'}
ssx={'H':2,'E':4,'C':1}  #numeric code for secondary structures

inp=open('./HE_HE.dis','r')
type='HE'
subdir=outd+'/c' #directory to output sequences
os.system('/bin/mkdir -p '+subdir)
iseq=0 #current index for output sequence

#read beginning singlet probabilities
singlets={'H':0,'E':0}
singletKeys=singlets.keys()
inp.readline() #(beginning singlet)
for i in range(2):
    singlet,prob=inp.readline().split()
    singlets[singlet]=float(prob)

#read beginning doublet probabilities
doublets={'HH':0,'EH':0,'HE':0,'EE':0}
doubletKeys=doublets.keys()
inp.readline() #(beginning doublet)
for i in range(4):
    doublet,prob=inp.readline().split()
    doublets[doublet]=float(prob)
#fabricate the singlet probability conditional on the previous singlet element
conditionalDoublet={}
for singlet in singletKeys: #previous singlet, the condition
    conditionalDoublet[singlet]={}
    norm=0
    for singlet2 in singletKeys: norm+=doublets[ singlet+singlet2 ]
    for singlet2 in singletKeys:
        conditionalDoublet[singlet][singlet2]=doublets[ singlet+singlet2 ]/norm

#read triplets
triplets={'HHH':0,'HHE':0,'HEH':0,'HEE':0,'EHH':0,'EHE':0,'EEH':0,'EEE':0}
tripletKeys=triplets.keys()
inp.readline() #(triplets)
for i in range(8):
    triplet,prob=inp.readline().split()
    triplets[triplet]=float(prob)
#fabricate the  singlet probability conditional on the previous doublet
conditionalTriplet={}
for doublet in doubletKeys: #previous doublet, the condition
    conditionalTriplet[doublet]={}
    norm=0
    for singlet2 in singletKeys: norm+=triplets[ doublet+singlet2 ]
    for singlet2 in singletKeys:
        conditionalTriplet[doublet][singlet2]=triplets[ doublet+singlet2 ]/norm

#print singlets
#print conditionalDoublet
#print triplets
#print conditionalTriplet
#sys.exit(1)
#read histrogram with number of elements
distrnSS=[0]*100  #initialize histogram of number of elements
inp.readline() #(HE)
inp.readline() #(distribution of number of strands in alpha/beta proteins)
line=inp.readline() #;print 'line=',line
while(line[0]!='d'):
    nSS,y=line.split()
    nSS=int(nSS)
    distrnSS[nSS]=float(y)
    line=inp.readline()
        
#find maximum of histogram
distrMax=0
nSSMax=0
for nSS in range(len(distrnSS)):
    if distrnSS[nSS] > distrMax:
        distrMax=distrnSS[nSS]
        nSSMax=nSS

#posterior distributions, for checking with experimental ones
postDistLeng={'E-E':[0]*(2*maxLength),'E-H-E':[0]*(2*maxLength)}
postN={'H':0,'E':0}
postSinglet={'H':0,'E':0}
postDoublet={'HH':0,'EH':0,'HE':0,'EE':0}
postTriplets={'HHH':0,'HHE':0,'HEH':0,'HEE':0,'EHH':0,'EHE':0,'EEH':0,'EEE':0}

#calculate the number of sequences we have to output per nSS bin:
nseqPernSS=[]
for p in distrnSS: nseqPernSS.append( int(p*Nseqs) )
remainder=Nseqs
for i in nseqPernSS: remainder -= i
nseqPernSS[nSSMax]+=remainder #add remainder to most expected, so that we have exactly Nseqs

#read distributions
for nSS in range(len(distrnSS)):
    if distrnSS[nSS]>0:
        #read length distribution from file
        histLength={}
        inp.readline()
        histLength['H']=inp.readline().split()
        inp.readline()
        histLength['E']=inp.readline().split()
        inp.readline()
        histLength['loop']=inp.readline().split()
        inp.readline()
        histLength['E-E']=inp.readline().split()
        inp.readline()
        histLength['E-H-E']=inp.readline().split()

        #make normalized histogram and cummulative histograms
        for kind in ('H','E','loop'):
            hist=histLength[kind]         #handy reference to histogram
            norm=0
            for i in range(maxLength):
                hist[i]=int(hist[i])
                norm+=hist[i]
            for i in range(maxLength): hist[i]=float(hist[i])/norm
        for kind in ('E-E','E-H-E'):
            hist=histLength[kind]         #handy reference to histogram
            norm=0
            for i in range(2*maxLength):
                hist[i]=int(hist[i])
                norm+=hist[i]
            for i in range(2*maxLength): hist[i]=float(hist[i])/norm
                

        #calculate average lengths
        avLs={'H':0.0,'E':0.0,'loop':0.0}
        for kind in ('H','E','loop'):
            hist=histLength[kind]
            for L in range(maxLength): avLs[kind]+=L*hist[L]
            
        #set expected number of helixes equal average composition
        #of residues in the helix and strand states
        #(number-helixes * average-helix-length=number-strands *
        #average-strand-helix)
        #ratio=avLs['E']/avLs['H'] #prob of choosing a helix over a strand
        #Produce a set of sequences according to the probabilities
        #set by the histogram.
        nToOut=nseqPernSS[nSS] #number of sequences to output with given number of elements
        print 'nSS=',nSS,' nToOut=',nToOut
        ntries=0
        while nToOut:
            remnant=nSS #number of strands yet to place in a sequence
            excess=0    #flag to check if we exceeded the sequence length
            elemL=0      #current sequence length
            elemChunks=[]   #length of each chunk
            nchunks=0   #current number of inserted secondary elements
            minloopL=0 #minimal number of residues in interior loops
            nhelix=0
            while 1: #loop to generate one sequence of secondary elements
                #select strand of helix, based on previouly placed elements
                kind='E' #initialize
                if nchunks==0: #we did not place any element or loop
                    if float(postSinglet['H'])/(1+iseq) < singlets['H']: kind='H' 
                elif nchunks==1: #we placed one element
                    prevSinglet=elemChunks[0][0] #previous element
                    if random()<conditionalDoublet[prevSinglet]['H']: kind='H'
                else:
                    prevDoublet=elemChunks[-2][0]+elemChunks[-1][0] #previous doublet
                    if random()<conditionalTriplet[prevDoublet]['H']: kind='H'                    
                    
                if kind=='E': #we selected strand
                    #select one of E, E-H, and E-H-E distributions
                    prevL=0
                    if nchunks>1 and elemChunks[-1][0]=='E': #prev element also a strand
                        prevL=elemChunks[-1][1] #length of previous strand
                        kind2='E-E' #;print kind2
                    elif nchunks>2 and elemChunks[-1][0]=='H' and elemChunks[-2][0]=='E':
                        kind2='E-H-E'
                        prevL=elemChunks[-2][1] #length of previous strand
                    else:
                        kind2='E'                    
                    hist=histLength[kind2] #select appropriate distribution
                    while 1: #cycle to accept a chunk
                        n=randint(3,11) #putative lenght of a strand
                        d=n-prevL #length difference with previous strand
                        if random()<hist[d]:
                            elemChunks.append( (kind,n) )
                            elemL+=n
                            nchunks+=1
                            break
                else:
                    nhelix+=1
                    hist=histLength[kind] #select appropriate distribution
                    while 1: #cycle to accept a chunk
                        n=randint(4,maxLength-1) #putative lenght of an element
                        if random()<hist[n]:
                            elemChunks.append( (kind,n) )
                            elemL+=n
                            nchunks+=1
                            break
                        
                minloopL=3*(nchunks-1) #current minimum number of residues in interior loops
                if elemL+minloopL>Ntarget: #Oops! we exceeded Ntarget residues
                    excess=1 ; break

                if kind=='E': remnant-=1 #one less strand to place
                if not remnant: break #we placed all nSS strands

            if excess or 2*elemL<Ntarget: #too many or too few residues in the secondary elements
                continue

            if nhelix<=2: continue #there should be at least two helixes
            
            print 'nSS=',nSS,' nchunks=',nchunks,' elemL=',elemL
            loopL=Ntarget-elemL #remanining residues to be placed in the loop conformation
            hist=histLength['loop']
            loopChunks=[]
            norm=0 #number of allocated residues in interior loops
            ntries=infinity
            remnant=nchunks-1 #number of interior loops to generate
            while ntries:
                while remnant:
                    n=randint(3,maxLength-1) #lenght of a loop
                    if random()<hist[n]:
                        loopChunks.append( ['C',n] ) ; norm+=n ; remnant-=1
                if norm>loopL: #try again with other loop
                    ntries-=1
                    loopChunks=[]
                    norm=0
                    remnant=nchunks-1
                    continue
                break

            if not ntries: #we failed in the allocation of the loops
                for i in range(nchunks-1): loopChunks.append(['C',3])
                norm=3*(nchunks-1)
            remnant=loopL-norm
            if remnant<3:
                addToLoops(remnant,loopChunks) ; Nt=0 ; Ct=0
            elif remnant<6:
                Nt=3 ; Ct=0 ; addToLoops(remnant-Nt,loopChunks)
            else:
                Nt=randint(3,remnant/2)
                Ct=randint(3,remnant/2)
                addToLoops(remnant-Nt-Ct,loopChunks)
            print 'nSS=',nSS,'left=',nToOut,' loopL=',loopL,'norm+Nt+Ct=',norm+Nt+Ct
                

            #retrieve sequence of elements
            #seqSS=''
            #for i in range(len(chunks)):
            #    if chunks[i][0]!='C': seqSS+=chunks[i][0]

            #use postSinglet also as general frequence of helix/strand
            #occurrence. We impose this condition because the amino
            #acid length constrain of Ntarget in combination with the
            #constraint of imposing nSS strands favors putting more
            #strands than helixes
            #nh=len(re.compile('H').findall(seqSS))
            #nb=len(re.compile('E').findall(seqSS))
            #if float(nhelix+nh)/(nall+nh+nb)<singlets['H']:
            #    ntries+=1
            #    continue               

            #update number of helixes and strands
            #nhelix+=nh
            #nall+=nh+nb
            #print  float(nhelix)/nall
            
            #update posterior distributions
            for i in range(1,nchunks):
                if elemChunks[i][0]=='E' and elemChunks[i-1][0]=='E':
                    d=elemChunks[i][1]-elemChunks[i-1][1]
                    postDistLeng['E-E'][d]+=1
            for i in range(2,nchunks):
                if elemChunks[i][0]=='E' and elemChunks[i-1][0]=='H' and  elemChunks[i-2][0]=='E':
                        d=elemChunks[i][1]-elemChunks[i-2][1]
                        postDistLeng['E-H-E'][d]+=1
            for i in range(nchunks):
                if elemChunks[i][0]=='E': postN['E']+=1
                elif elemChunks[i][0]=='H': postN['H']+=1

            postSinglet[ elemChunks[0][0] ]+=1 #update post beginning singlet

            postDoublet[ elemChunks[0][0]+elemChunks[1][0] ]+=1 #update post beginning doublet
            
            #for ss in range(0,len(seqSS)-2): #update post triplets
            #    triplet=seqSS[ss:ss+3]
            #    postTriplets[triplet]+=1

            #output the sequence
            if len(loopChunks)!=nchunks-1: Bye('ERROR')
            
            subsubdir=subdir+'/c'+'%04d'%(iseq)
            os.system('/bin/mkdir -p '+subsubdir)
            pf=open(subsubdir+'/seq.dat','w')
            m=1 #current amino acid
            for i in range(Nt): #output the N-terminal
                pf.write('%5d   VAL%5d   10\n'%(m,ssx['C'])) #N-terminal
                m+=1
            for j in range(nchunks-1):
                kind=elemChunks[j][0]
                extent=elemChunks[j][1]
                for i in range(extent):
                    pf.write('%5d   VAL%5d   10\n'%(m,ssx[kind])) ; m+=1
                #print 'loopChunks[',j,'][1]=',loopChunks[j][1]
                extent=loopChunks[j][1]
                for i in range(extent):
                    pf.write('%5d   VAL%5d   10\n'%(m,ssx['C'])) ; m+=1
            kind=elemChunks[nchunks-1][0]
            extent=elemChunks[nchunks-1][1]
            for i in range(extent):
                pf.write('%5d   VAL%5d   10\n'%(m,ssx[kind])) ; m+=1
            for i in range(Ct):
                pf.write('%5d   VAL%5d   10\n'%(m,ssx['C'])) #C-terminal
                m+=1
            pf.close()
            iseq+=1 #one more sequence that has beeen output
            #print alph2dir[type]+'%04d %5d %3d'%(iseq,ntries,nToOut),' nSS=',nSS
            nToOut-=1 #one less sequence to generate

#output posterior distributions
for x in ('E-E','E-H-E'):
    print x+' posterior distribution'
    post=postDistLeng[x]
    norm=sum(post) #add all bins of the histogram
    for i in range(1-maxLength,maxLength-1): post[i]=float(post[i])/norm #normalize distribution
    for i in range(1-maxLength,maxLength-1): print '%3d %5.3lf'%(i,post[i])
print 'singlet probability'
norm=sum(postN.values())
for x in ('E','H'):
    print '%s %5.3lf'%(x,float(postN[x])/norm)
print 'beginning singlet probability'
norm=sum(postSinglet.values())
for x in ('E','H'):
    print '%s %5.3lf'%(x,float(postSinglet[x])/norm)
print 'beginning doublet probability'
norm=sum(postDoublet.values())
for x in doubletKeys:
    print '%s %5.3lf'%(x,float(postDoublet[x])/norm)
print 'triplets  probability'
norm=sum(postTriplets.values())
for x in tripletKeys:
    print '%s %5.3lf'%(x,float(postTriplets[x])/norm)
print nnn,mmm
sys.exit(0)
