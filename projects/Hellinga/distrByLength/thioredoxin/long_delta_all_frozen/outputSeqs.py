#!/usr/bin/python
import sys,os,re
from utilities.small_utilities import chomp,Bye
from stats.statistics import average,standard_dev,excess_kurtosis,skewness
from inputArgs.inputArgs import inpHand
from random import random,randint

infinity=1000
maxLength=40
Nseqs=400
Ntarget=105
outd='.'
alph2dir={'H':'a','E':'b'}
ssx={'H':2,'E':4,'C':1}  #numeric code for secondary structures

#insert remnant residues in the loops
def addToLoops(remnant,loopChunks):
    nloops=len(loopChunks)
    for i in range(remnant):
        n=randint(0,nloops-1)
        loopChunks[n]+=1

inp=open('./HE.dis','r')
for type in ('H','E'):
    subdir=outd+'/'+alph2dir[type] #directory to output sequences
    os.system('/bin/mkdir -p '+subdir)
    iseq=0 #current index for output sequence

    #read histrogram with number of elements
    distrnSS=[0]*100  #initialize histogram of number of elements
    inp.readline() #(H)
    inp.readline() #(distribution of number of elements)
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

    postDistLengh=[0]*maxLength #posterior distr. of element lengths (for all number of elemen)
    distDiffLeng=[0]*(2*maxLength) #distribution of length difference between consecutive elements
    
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
            histLength['SS']  =inp.readline().split()
            inp.readline()
            histLength['loop']=inp.readline().split()
            inp.readline()
            histLength['diff']=inp.readline().split()

            if type=='H':continue #we won't process H, but only E
            
            #make normalized histogram and cummulative histograms
            cummHistLength={'SS':[0]*maxLength,'loop':[0]*maxLength}
            for kind in ('SS','loop'):
                hist=histLength[kind]         #handy reference to histogram
                cummHist=cummHistLength[kind]
                norm=0
                for i in range(maxLength):
                    hist[i]=int(hist[i])
                    norm+=hist[i]
                    cummHist[i]=norm
                for i in range(maxLength): cummHist[i]=float(cummHist[i])/norm
                for i in range(maxLength): hist[i]=float(hist[i])/norm
            hist=histLength['diff']
            norm=0
            for i in range(2*maxLength):
                hist[i]=float(hist[i])
                norm+=hist[i]
            for i in range(2*maxLength): hist[i]=hist[i]/norm
            #print 'histLength'; print histLength; print 'cummHistLength'; print cummHistLength

            #Produce a set of sequences according to the probabilities
            #set by the histogram.
            minloopL=(nSS-1)*3 #minimum number of residues in loop conformation
            nToOut=nseqPernSS[nSS] #number of sequences to output with given number of elements
            #print 'nToOut=',nToOut
            while nToOut:
                elemL=0       #current sequence length
                elemChunks=[]
                prevLSS=0
                ntries=infinity
                remnant=nSS #number of elements yet to place in a sequence
                while remnant and ntries: #loop to generate one sequence
                    #derive fist strand from length distr. Derive any helix in this way too
                    if type=='H' or (type=='E' and len(elemChunks)==0):
                        hist=histLength['SS']
                        while 1: #cycle to accept a chunk
                            n=randint(3,maxLength-1) #lenght of a chunk (element or loop)
                            if random()<hist[n]:
                                elemChunks.append(n) ; elemL+=n ; prevLSS=n ; break

                    #derive strand from length difference distr. from consecutive strands
                    else:
                        hist=histLength['diff']
                        hist2=histLength['SS']
                        while 1: #cycle to accept a chunk
                            n=randint(3,maxLength-1) #no strand smaller than three residues
                            if hist2[n]>0:
                                d=prevLSS-n
                                if random()<hist[d]:
                                    elemChunks.append(n) ; elemL+=n ; prevLSS=n ;break
                    remnant-=1 #one less secondary element

                    if elemL+minloopL>Ntarget:
                        remnant=nSS #begin again
                        elemChunks=[]
                        elemL=0
                        prevLSS=0
                        ntries-=1

                if 2*elemL<Ntarget: continue #too few residues
                print 'nSS=',nSS,'nToOut=',nToOut,'elemL/Ntarget=',float(elemL)/Ntarget
                if not ntries: continue 

                #apportion loops
                loopL=Ntarget-elemL #number of remaining residues
                hist=histLength['loop']
                mtries=infinity
                norm=0
                loopChunks=[]   #length of each chunk
                while mtries:
                    remnant=nSS-1 #number of interior loops to apportion
                    while remnant:
                        n=randint(3,maxLength-1) #lenght of a chunk (element or loop)
                        if random()<hist[n]:
                            loopChunks.append(n)
                            norm+=n
                            remnant-=1
                    if norm>loopL: #allocated residues exceed loopL
                        loopChunks=[]
                        norm=0                
                        mtries-=1
                    else: break
                if not mtries: #we were not successful to allocate the loops from the histogram
                    loopChunks=[3]*(nSS-1) #initialize every loop as three residues
                    norm=3*(nSS-1)                    
                remnant=loopL-norm
                if remnant<3:
                    addToLoops(remnant,loopChunks) ; Nt=0 ; Ct=0
                elif remnant<6:
                    Nt=3 ; Ct=0 ; addToLoops(remnant-Nt,loopChunks)
                else:
                    Nt=randint(3,remnant/2)
                    Ct=randint(3,remnant/2)
                    addToLoops(remnant-Nt-Ct,loopChunks)
                print 'nSS=',nSS,'loopL=',loopL,'norm+Nt+Ct=',norm+Nt+Ct
                #update the posterior distribution of element lengths
                i=0
                while i<len(elemChunks):
                    postDistLengh[ elemChunks[i] ]+=1
                    i+=1

                #update the histogram of the difference in consecutive elements
                i=0
                while i<len(elemChunks)-1:
                    #NOTE: chunks[i+2]-chunks[i] can be negative, but python does handle negative
                    #argument in list indexes. We will take advantage and use
                    #distDiffLeng[0:maxLength] for positive differences and
                    #distDiffLeng[1-maxLength:]==distDiffLeng[maxLength:] for negative differences
                    distDiffLeng[ elemChunks[i+1]-elemChunks[i] ]+=1
                    i+=1
                
                #output the sequence
                subsubdir=subdir+'/'+alph2dir[type]+'%04d'%(iseq)
                os.system('/bin/mkdir -p '+subsubdir)
                pf=open(subsubdir+'/seq.dat','w')
                m=1 #current amino acid
                for i in range(Nt): #output the N-terminal
                    pf.write('%5d   VAL%5d   10\n'%(m,ssx['C'])) #N-terminal
                    m+=1
                #print 'nloops=',nloops,'nSS-1=',nSS-1
                for j in range(nSS-1):
                    for i in range(elemChunks[j]):
                        pf.write('%5d   VAL%5d   10\n'%(m,ssx[type])) ; m+=1
                    for i in range(loopChunks[j]):
                        pf.write('%5d   VAL%5d   10\n'%(m,ssx['C']))  ; m+=1
                for i in range(elemChunks[nSS-1]):
                    pf.write('%5d   VAL%5d   10\n'%(m,ssx[type])) ; m+=1
                for i in range(Ct):
                    pf.write('%5d   VAL%5d   10\n'%(m,ssx['C'])) #C-terminal
                    m+=1
                pf.close()
                iseq+=1 #one more sequence that has beeen output
                #print alph2dir[type]+'%04d %5d %3d'%(iseq,ntries,nToOut),' nSS=',nSS
                nToOut-=1 #one less sequence to generate

    if type=='H':continue #we won't process H, but only E
    
    #normalize the posterior distribution of element lengths
    norm=sum(postDistLengh)
    #Print the posterior distribution of element lengths
    print 'posterior distribution of '+type+' lengths'
    for i in range(maxLength):
        print '%2d %4.2lf'%(i,float(postDistLengh[i])/norm)

    #normalize the distribution of length difference between consecutive elements
    norm=sum(distDiffLeng)
    #Print the distribution of length difference between consecutive elements
    print 'distribution of length difference between consecutive '+type
    for i in range(1-maxLength,maxLength):
        print '%2d %5.3lf'%(i,float(distDiffLeng[i])/norm)
    
sys.exit(0)
