#!/usr/bin/python
import sys,os,re
from utilities.small_utilities import chomp,Bye
from stats.statistics import average,standard_dev,excess_kurtosis,skewness
from inputArgs.inputArgs import inpHand
from random import random,randint,shuffle

infinity=1E3
maxLength=40
Nseqs=600
Ntarget=150
outd='.'
alph2dir={'H':'a','E':'b'}
ssx={'H':2,'E':4,'C':1}  #numeric code for secondary structures
minHelixL=4

#insert remnant residues in the loops
def addToLoops(remnant,loopChunks):
    nloops=len(loopChunks)
    for i in range(remnant):
        n=randint(0,nloops-1)
        loopChunks[n]+=1

#remove a few residues from the secondary elements
def removeFromChunks(excess,iToOut,nSS,chunks):
    n0=iToOut*nSS
    #print 'before '+`sum(chunks[n0:n0+nSS])`
    #print chunks[n0:n0+nSS]
    while excess:
        n=randint(0,nSS-1)
        if chunks[n0+n]>minHelixL:
            #print 'excess=',excess,' sum=',`sum(chunks[n0:n0+nSS])`,' n=',n,'chunks[n0+n]=',chunks[n0+n]
            chunks[n0+n]=chunks[n0+n]-1
            excess-=1
            #print 'chunks[n0+n]=',chunks[n0+n],' excess=',excess
            #print chunks[n0:n0+nSS]                
    #print 'excess=',excess,' sum=',`sum(chunks[n0:n0+nSS])`
    #print chunks[n0:n0+nSS]
    #print 'after '+`sum(chunks[n0:n0+nSS])`
    
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
            minloopL=(nSS-1)*3 #minimum number of res in loop conformation for one seq.
            nToOut=nseqPernSS[nSS] #number of sequences to output with given number of elements
            nelements=nToOut*nSS   #number of elements to pull from the length distribution
            hist=histLength['SS']
            minpredEl=expEl+0.02   #minimum helical content
            #maxpredEl=float(Ntarget-minloopL)/Ntarget #maximum helical content
            maxpredEl=minpredEl+0.02 #if too high, then we can't apportion the helixes
            maxextent=int(maxLength-1-2*nSS) #avoid trombones
            print 'nSS=',nSS,' nToOut=',nToOut,' expEl=',expEl,' minpredEl=',minpredEl,' maxpredEl=',maxpredEl,' maxextent=',maxextent
            while nToOut:
                predEl=0 #initialize predicted element content

                #produce nelements from the length distribution
                chunks=[]                            
                remnant=nelements
                while remnant:
                    n=randint(3,maxextent)
                    if random()<hist[n]:
                        chunks.append(n+1) #insert an element of length n+1
                        remnant-=1
                predEl=float(sum(chunks))/(nToOut*Ntarget) #compute predicted element content
                #print predEl
                if predEl<minpredEl or predEl>maxpredEl: continue
                print 'predEl=',predEl,' expEl=',expEl,' minpredEl=',minpredEl,' maxpredEl=',maxpredEl
                #apportion the nelements among the nToOut sequences
                
                #we need first set of nSS below Ntarget
                if sum(chunks[0:nSS])+minloopL > Ntarget: continue            

                ntries=0
                iToOut=1
                ichunk=nSS
                while iToOut < nToOut:
                    elemL=sum(chunks[ichunk:ichunk+nSS]) #next set of nSS
                    #if this set of nSS is too long, then pick chunks
                    #from previous assigments and interchange
                    while elemL+minloopL > Ntarget and ntries<infinity:
                        #select one previous sequence and one chunk
                        jToOut=randint(0,iToOut-1)    #select previous sequence
                        jchunk=jToOut*nSS             
                        norm=sum(chunks[jchunk:jchunk+nSS]) #length of previous sequence
                        if norm+minloopL+elemL+minloopL>Ntarget: #shuffling will solve nothing
                            ntries+=1
                            continue
                        #try to exchange one chunk of iToOut with one chunk of jToOut
                        infinity2=nSS*nSS
                        mtries=0
                        newelemL=0
                        while mtries<infinity2:
                            m=jchunk+randint(0,nSS-1) #chunk from previous set of nSS chunks
                            lm=chunks[m]             #size of the chunk
                            n=ichunk+randint(0,nSS-1)    #chunk from current set of nSS chunks
                            ln=chunks[n]             #size of the chunk
                            newnorm=norm-lm+ln
                            newelemL=elemL-ln+lm
                            if newnorm+minloopL<=Ntarget and newelemL+minloopL<=Ntarget:
                                chunks[n]=lm
                                chunks[m]=ln
                                break
                            else: mtries+=1
                        if mtries<infinity2: elemL=newelemL
                        else: ntries=+1
                    #print nToOut-iToOut,elemL,elemL+minloopL
                    ichunk+=nSS #advance nSS chunks
                    iToOut+=1   
                    
                if ntries==infinity:
                    print 'Entering putative removal of excess residues'
                    #be permissive if violation is not too big
                    excess=0
                    for iToOut in range(nToOut):
                        helicalL=sum(chunks[iToOut*nSS:(iToOut+1)*nSS])
                        excess=helicalL+minloopL-Ntarget
                        if excess>6: break
                        elif excess>0:
                            #print 'removing '+`excess`+' from sequence '+`iToOut`+' (left='+`nToOut-iToOut`+') which has '+`helicalL`+' helical residues'
                            removeFromChunks(excess,iToOut,nSS,chunks)
                            #print 'now ',`sum(chunks[iToOut*nSS:(iToOut+1)*nSS])+minloopL`
                    if excess>6:
                        print 'ntries above infinity and excess='+`excess`
                        continue

                #for iToOut in range(nToOut):
                #    print nSS,nToOut-iToOut,sum(chunks[iToOut*nSS:(iToOut+1)*nSS])+minloopL
                   
                #for every sequence, apportion a set of interior loops, then print sequence
                ntries=0
                hist==histLength['loop']
                for iToOut in range(nToOut):
                    elemL=sum(chunks[iToOut*nSS:(iToOut+1)*nSS])
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
                            n=randint(3,maxLength-1) #min loop size is 3, TASSER require
                            if random()<hist[n]:
                                loopChunks.append(n)
                                remnant-=1
                        ntries+=1
                        if ntries>infinity:
                            loopChunks=[]
                            maxll=int(loopL/(nSS-1))
                            for i in range(nSS-1): loopChunks.append(randint(3,maxll))
                            print 'loop apportion reached infinity. nSS=',nSS,' norm=',norm,' maxll=',maxll
                        norm=sum(loopChunks)
                        if norm>loopL: continue #too many
                        break
                
                    #apportion termini, taking into account minimum size for loop is 3
                    remnant=loopL-norm #;print 'remnant=',remnant
                    Nt=0
                    Ct=0
                    if remnant<3:
                        addToLoops(remnant,loopChunks)
                    elif remnant<6:
                        Nt=3
                        addToLoops(remnant-Nt,loopChunks)
                    else:
                        Nt=randint(3,remnant/2)
                        Ct=randint(3,remnant/2)
                        addToLoops(remnant-Nt-Ct,loopChunks)
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
