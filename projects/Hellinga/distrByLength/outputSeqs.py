#!/usr/bin/python
import sys,os,re
from utilities.small_utilities import chomp,Bye
from stats.statistics import average,standard_dev,excess_kurtosis,skewness
from inputArgs.inputArgs import inpHand
from random import random,randint

maxLength=40
Nseqs=400
Ntarget=150
outd='.'
alph2dir={'H':'a','E':'b'}
ssx={'H':2,'E':4,'C':1}  #numeric code for secondary structures

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
            nToOut=nseqPernSS[nSS] #number of sequences to output with given number of elements
            #print 'nToOut=',nToOut
            ntries=0
            while nToOut:
                remnant=nSS #number of elements yet to place in a sequence
                excess=0    #flag to check if we exceeded the sequence length
                norm=0       #current sequence length
                chunks=[]   #length of each chunk
                prevLSS=0
                while 1: #loop to generate one sequence
                    #derive fist strand from length distr. Derive any helix in this way too
                    if type=='H' or (type=='E' and len(chunks)==0):
                        hist=histLength['SS']
                        while 1: #cycle to accept a chunk
                            n=randint(0,maxLength-1) #lenght of a chunk (element or loop)
                            if random()<hist[n]:
                                chunks.append(n) ; norm+=n ; prevLSS=n ; break

                    #derive strand from length difference distr. from consecutive strands
                    else:
                        hist=histLength['diff']
                        hist2=histLength['SS']
                        while 1: #cycle to accept a chunk
                            n=randint(2,maxLength-1) #no strand smaller then three residues
                            if hist2[n]>0:
                                d=prevLSS-n
                                if random()<hist[d]:
                                    chunks.append(n) ; norm+=n ; prevLSS=n ;break
                        
                    if norm>Ntarget:
                        excess=1 ; break
                        
                    remnant-=1 #one less element to place
                    if not remnant: break #we placed all nSS elements
                    
                    hist=histLength['loop']
                    while 1: #cycle to accept a chunk
                        n=randint(0,maxLength-1) #lenght of a chunk (element or loop)
                        if random()<hist[n]:
                            chunks.append(n) ; norm+=n ; break
                    if norm>Ntarget:
                        excess=1 ; break

                if excess: #we failed trying to generate a sequence. Start again                 
                    ntries+=1
                    continue 

                #update the posterior distribution of element lengths
                i=0
                while i<len(chunks):
                    postDistLengh[ chunks[i] ]+=1
                    i+=2 #remember chunks lists the SS-loop-SS-loop-...-SS sequence

                #update the histogram of the difference in consecutive elements
                i=0
                while i<len(chunks)-1:
                    #NOTE: chunks[i+2]-chunks[i] can be negative, but python does handle negative
                    #argument in list indexes. We will take advantage and use
                    #distDiffLeng[0:maxLength] for positive differences and
                    #distDiffLeng[1-maxLength:]==distDiffLeng[maxLength:] for negative differences
                    distDiffLeng[ chunks[i+2]-chunks[i] ]+=1
                    i+=2 #remember chunks lists the SS-loop-SS-loop-...-SS sequence
                
                #randomly divide sequence remainder between termini
                remnant=Ntarget-norm
                Nt=randint(0,remnant) ; Ct=remnant-Nt

                #output the sequence
                subsubdir=subdir+'/'+alph2dir[type]+'%04d'%(iseq)
                os.system('/bin/mkdir -p '+subsubdir)
                pf=open(subsubdir+'/seq.dat','w')
                m=1 #current amino acid
                for i in range(Nt): #output the N-terminal
                    pf.write('%5d   VAL%5d   10\n'%(m,ssx['C'])) #N-terminal
                    m+=1
                M=len(chunks) #number of elements plus interior loops
                mch=0         #number of elements plus interior loops placed
                while 1:
                    for i in range(chunks[mch]):
                        pf.write('%5d   VAL%5d   10\n'%(m,ssx[type])) ; m+=1
                    mch+=1
                    if mch==M: break #remember last item in chunk is not a loop, but an element
                    for i in range(chunks[mch]):
                        pf.write('%5d   VAL%5d   10\n'%(m,ssx['C']))  ; m+=1
                    mch+=1 
                for i in range(Ct):
                    pf.write('%5d   VAL%5d   10\n'%(m,ssx['C'])) #C-terminal
                    m+=1
                pf.close()
                iseq+=1 #one more sequence that has beeen output
                #print alph2dir[type]+'%04d %5d %3d'%(iseq,ntries,nToOut),' nSS=',nSS
                ntries=0
                nToOut-=1 #one less sequence to generate

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
