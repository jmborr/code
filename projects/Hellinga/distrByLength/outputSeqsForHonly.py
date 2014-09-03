#!/usr/bin/python
import sys,os,re
from utilities.small_utilities import chomp,Bye
from stats.statistics import average,standard_dev,excess_kurtosis,skewness
from inputArgs.inputArgs import inpHand
from random import random,randint,shuffle

infinity=1E3
maxLength=40
Nseqs=400
Ntarget=150
outd='.'
alph2dir={'H':'a','E':'b'}
ssx={'H':2,'E':4,'C':1}  #numeric code for secondary structures

inp=open('./HE.dis','r')
for type in ('H'):
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

            #calculate expected element content
            expEl=0.0
            hist=histLength['SS']
            for L in range(len(hist)): expEl+=L*hist[L] #average length
            expEl=expEl*nSS/Ntarget #expected element content
            
            #Produce a set of sequences according to the probabilities
            #set by the histogram.
            minloopL=(nSS-1)*4 #minimum number of res in loop conformation for one seq.
            nToOut=nseqPernSS[nSS] #number of sequences to output with given number of elements
            nelements=nToOut*nSS   #number of elements to pull from the length distribution
            hist=histLength['SS']
            minpredEl=expEl+0.02   #minimum helical content
            #maxpredEl=float(Ntarget-minloopL)/Ntarget #maximum helical content
            maxpredEl=minpredEl+0.02 #if too high, then we can't apportion the helixes
            print 'nSS=',nSS,' nToOut=',nToOut,' expEl=',expEl,' minpredEl=',minpredEl,' maxpredEl=',maxpredEl

            while nToOut:
                predEl=0 #initialize predicted element content

                #produce nelements from the distribution, with an
                #average element content 10% above expected
                chunks=[]                            
                remnant=nelements
                while remnant:
                    n=randint(3,maxLength-1)
                    if random()<hist[n]:
                        chunks.append(n+1) #insert an element of length n+1
                        remnant-=1
                predEl=float(sum(chunks))/(nToOut*Ntarget) #compute predicted element content
                #print predEl
                if predEl<minpredEl or predEl>maxpredEl: continue
                print 'predEl=',predEl,' expEl=',expEl,' minpredEl=',minpredEl,' maxpredEl=',maxpredEl
                
                #apportion the nelements among the nToOut sequences 
                ntries=0
                toolong=1
                while toolong:
                    toolong=0 #flag signals if one sequence overflows Ntarget
                    shuffle(chunks) #random shuffle the list of elements
                    for i in range(nToOut):
                        elemL=0 #initialize number of res in element conformation for one sequence
                        for j in range(nSS): elemL+=chunks[i*nSS+j]
                        if elemL+minloopL > Ntarget:
                            toolong=1
                            break
                    if ntries>infinity:break #these nelements are no good
                    ntries+=1
                    if ntries%1000==0: print 'ntries=',ntries
                print 'Ntarget=',Ntarget,' elemL+minloopL=',elemL+minloopL
                if ntries>infinity:
                    print 'ntries above infinity'
                    continue
                #print 'apportion of the nelements: ntries=',ntries
                
                #for every sequence, apportion a set of interior loops, then print sequence
                ntries=0
                hist==histLength['loop']
                for iToOut in range(nToOut):
                    elemL=0
                    for j in range(nSS): elemL+=chunks[iToOut*nSS+j]
                    loopL=Ntarget-elemL #initialize number of residues in loop conformation
                    norm=0 #initialize number of residues in interior loops
                    print 'nSS=',nSS,' left=',nToOut-iToOut,' minloopL=',minloopL,' loopL=',loopL,' elemL=',elemL
                    #print 'seq #',iToOut,' res in loop conf=',loopL
                    ntries=0
                    while 1:
                        if loopL<minloopL+nSS: ntries=infinity #too close to minloopL
                        loopChunks=[]
                        remnant=nSS-1 #number of interior loops to assign
                        while remnant:
                            n=randint(3,maxLength-1)
                            if random()<hist[n]:
                                loopChunks.append(n)
                                remnant-=1
                        ntries+=1
                        if ntries%1000==0: print 'ntries=',ntries
                        if ntries>infinity:
                            loopChunks=[]
                            maxll=int(loopL/(nSS-1))
                            for i in range(nSS-1): loopChunks.append(randint(3,maxll))
                            norm=sum(loopChunks)
                            print 'loop apportion reached infinity. nSS=',nSS,' norm=',norm,' maxll=',maxll
                            break #these nelements are no good
                        norm=sum(loopChunks)
                        if norm>loopL: continue #too many or too few
                        break
                
                    #randomly apportion termini
                    remnant=loopL-norm #;print 'remnant=',remnant
                    Nt=randint(0,remnant)
                    Ct=remnant-Nt
                                         
                    #output the sequence
                    buf='' #buffer
                    subsubdir=subdir+'/'+alph2dir[type]+'%04d'%(iseq)
                    os.system('/bin/mkdir -p '+subsubdir)
                    pf=open(subsubdir+'/seq.dat','w')
                    m=1 #current amino acid
                    for j in range(Nt): #output the N-terminal
                        buf+='%5d   VAL%5d   10\n'%(m,ssx['C']) #N-terminal
                        m+=1
                    for j in range(nSS-1): #output pairs of element/loop
                        elemenL=chunks[iToOut*nSS+j]
                        for k in range(elemenL):
                            buf+='%5d   VAL%5d   10\n'%(m,ssx[type])
                            m+=1
                        loopL=loopChunks[j]
                        for k in range(loopL):
                            buf+='%5d   VAL%5d   10\n'%(m,ssx['C'])
                            m+=1
                    elemenL=chunks[iToOut*nSS+nSS-1] #output last element
                    for k in range(elemenL):
                        buf+='%5d   VAL%5d   10\n'%(m,ssx[type])
                        m+=1
                    for j in range(Ct): #output the C-terminal
                        buf+='%5d   VAL%5d   10\n'%(m,ssx['C']) #C-terminal
                        m+=1                        
                    pf.write(buf) #dump buffer to file
                    pf.close()
                    iseq+=1                    
                nToOut=0
    
sys.exit(0)
