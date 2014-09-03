import sys,numpy,pickle,copy,pdb

"""Gaussian network model. Two nodes interact via a spring
V_{i,j}=1/2 K_{i,j} ( r_{i,j} - <r_{i,j}> )^2, and
        r_{i,j}=relative vector from i to j, and
        K_{i,j} has NxN components
This interaction is a simplification over the more generalized
V_{i,j}=1/2 Sum_m Sum_l K_{i,j,l,m} ( r_{i,l}-r_{j,m}-<r_{i,l}-r_{j,m}> ),
        where r_{i,l} is component 'l' of position vector r_i, and
        K_{i,j,l,m} is a 3x3 spring tensor for a given {i,j} pair, that is,
        a K_{i,j,l,m} has a total of 3Nx3N components
In our model, the interparticle interaction is isotropic

NOTE: THE RESULTING QUASIHARMONIC MODES ARE ONE-DIMENSIONAL:
Example: if at t=0, only  mode 'k' is excited, then
         sqrt(m_i)*r_i=evecs_{i,k}q_k, where q_k is the initial vector
         of the mode, and r_i is the displacement vector of particles 'i'.
         Now notice that there is no implicit sum in 'k', thus vector r_i
         is just proportional to vector q_k. Thus, all particles will
         oscillate along the direction defined by q_k. This is a direct
         consequence of the isotropy of the interparticle potential.
"""

#fit observed histogram to a normal distribution by least squares
def normalmodel(param,s):
    a,b=param
    return a*numpy.exp(-b*s*s)

class gn3D:

    def __init__(self,np):
        self.np=np       #number of particles

    def pickleDump(self,pout):
        """save object with picke"""
        pickle.dump(self,pout,True)

    def pickleLoad(self,pin):
        """load object with pickle"""
        tmp=pickle.load(pin)
        self.__dict__.update(tmp.__dict__)
        
    def loadM(self,M): #load masses and precompute some
        self.M=M
        #Modification when SHAKE performed in the hydrogen atoms
        #print 'Debugging modification in loadM'
        #for i in range(len(self.M)):
        #    if self.M[i]==1.008: self.M[i]=1.42835
        #end
        self.sM=numpy.sqrt(M)  #sqrt(M)
        self.isM=1/self.sM     #1/sqrt(M)
        self.rMM=numpy.outer(self.isM,self.isM) #isM x isM 
        
    def findVarCovar(self,iterator):
        """find variance-covariance matrix of mass-weigthed coordinates
        'iterator' is a function  providing a fresh coordinate snapshot
        in the shape of a numpy.array(N,3)
        """
        assert hasattr(self,'M'), 'Run first loadM'
        avF=numpy.zeros(self.np*3).reshape(self.np,3);
        stdF=numpy.zeros(self.np*self.np).reshape(self.np,self.np)
        frame=iterator(); iframe=0
        while frame.any():
            frame*=self.sM[:,numpy.newaxis] #weight coords with sqrt(m_i)
            #stdF+= <r_i * r_j> vector product for all (i,j) pairs
            stdF+=numpy.inner(frame,frame);avF+=frame
            frame=iterator(); iframe+=1; #load next frame
            #if not iframe%1000: print iframe

            #print 'only for debugging'
            #if iframe==10000: break
            #end of only for debugging
            
        avF/=iframe;
        self.varcov= stdF/iframe - numpy.inner(avF,avF)
        self.avF=avF*self.isM[:,numpy.newaxis]

    def findModes(self):
        """find normal modes given a variance-covariance matrix
        varcovar=E( (x-mu)^t (x-mu) )
        evals,evecs are np, (np x np) matrices, respectively"""
        from numpy.linalg import eigh
        assert hasattr(self,'varcov'), 'Run first findVarCovar'
        #varcov = evecs * evals * evecs^t
        [self.evals,self.evecs]=eigh(self.varcov) #diagonalize varcovar
        #order from higher to smaller correlation
        perm=numpy.argsort(self.evals)[::-1]
        self.evals=self.evals[perm]
        self.evecs=self.evecs[:,perm]

    def findForceConstants(self,kT=1.0):
        """find matrix of force constants between coordinates.
        We divide by '3' because of the isotropic average
        We multiply by self.sM to remove the mass weighting
        sqrt(M)*x=evecs*q
        Last two force constant correspond to translations!"""
        if not hasattr(self,'evals'): self.findModes()
        self.FC=numpy.diag( 1/self.evals )
        self.FC[-1]=0.0; self.FC[-2]=0 #rotation and translation
        self.FC=numpy.dot(self.FC,self.evecs.transpose())
        self.FC=numpy.dot(self.evecs,self.FC)
        self.FC=(kT/3.0)*numpy.outer(self.sM,self.sM)*self.FC
            
    def projectTrajectory(self,iterator,pout):
        """project trajectory frames onto normal modes

         'iterator' is a function  providing a fresh coordinate snapshot
        in the shape of a numpy.array(N,3)

        'pout' is a pointer to file where we store the frames in the
        normal mode system of reference
        """
        from utilities.readingWritingFiles import write_from_numpy
        if not hasattr(self,'evals'): self.findModes()
        qframe=numpy.zeros(self.np*3).reshape(self.np,3)
        frame=iterator(); iframe=0; buf=''; sizeMax=50E6 #50MB
        #q=evecs^t*sqrt(M)*x
        U=self.evecs.transpose()*self.sM[numpy.newaxis:,]
        max=0;min=0
        while frame.any():
            qframe=numpy.dot(U,frame-self.avF)
            if qframe.max()>max: max=qframe.max()
            if qframe.min()<min: min=qframe.min()

            buf+=write_from_numpy(None,qframe,format=' %10.6f',ncols=7)
            if len(buf)>sizeMax: #flush the buffer
                pout.write(buf); buf=''
            frame=iterator(); iframe+=1; #load next frame
            if not iframe%1000: print iframe,max,min
            #if iframe==49999: pdb.set_trace()
            #only for debugging
            #if iframe==3000: break
            #end of only for debugging
        pout.write(buf) #flush remaining frames
        pout.close()
            
    def eflux(self,xiterator,viterator,avF=numpy.zeros(None)):
        """find net energy flux between particles
        'xiterator' returns coordinate frame
        'viterator' returns velocity frame
        flux=(3/4)*k_{ij}*(x_i*v_j-v_i*x_j), and  x_i:displacement vector
        """
        if not avF.any(): avF=self.avF
        if not hasattr(self,'FC'): self.findForceConstants()
        netflux=numpy.zeros(self.np*self.np).reshape(self.np,self.np)
        xframe=xiterator(); vframe=viterator(); iframe=0
        while xframe.any():
            xframe-=avF
            commuter=numpy.inner(xframe,vframe)
            commuter=commuter-commuter.transpose()
            netflux+=self.FC*commuter
            xframe=xiterator(); vframe=viterator(); iframe+=1
            if not iframe%1000: print iframe

            #print 'only for debugging'
            if iframe==10000: break
            #end of only for debugging

        self.netflux=(3.0/4)*netflux
        return self.netflux

    def chisquare(self,iterator,rewind):
        """evaluate the distribution of mode components to the
        theoretical normal distribution.

        The theoretical probability distribution for a single node is
        given by:
        P(\vec{q_i})=[2\pi/k]^(3/2)*exp(-1/2*k*\vec{q_i}^2)=
                    =[2\pi/k]^(3/2)*exp(-1/2*k(q_ix^2+q_iy^2+q_iz^2))

        Each of the q_ix,q_iy,q_iz should be normally distributed with
        zero average and variance 1/k. Thus, we define variable
        q_i=(q_ix+q_iy+q_iz)/sqrt(3) which should also be normally
        distributed and with the same variance.

        We have np q_i variables, one for each mode. We will find the
        pvalue of their histograms with respect to the theoretical
        histograms.

        'iterator' provides us with a new frame
        'rewind' go to beginning of trajectory
        """
        from histograms import histogram
        sphav=1/numpy.sqrt(3) #spherical average
        #find range of q's for each mode
        qmM=numpy.ones(self.np*3).reshape(3,self.np)
        frame=sphav*iterator().sum(axis=1);
        qmM*=frame[numpy.newaxis,:]
        frame=iterator(); iframe=0
        while frame.any():
            qmM[2]=sphav*frame.sum(axis=1)
            qmM[0]=qmM.min(axis=0); qmM[1]=qmM.max(axis=0)
            frame=iterator(); iframe+=1;
            if not iframe%1000: print iframe
        qmM=qmM.transpose()
        #read again the mode trajectory & generate histograms
        rewind();
        cs=1000 #reach trajectory in chunks
        nbins=min(100,iframe/100) #number of bins per histogram
        chunk=numpy.zeros(cs*self.np).reshape(cs,self.np)
        hist=numpy.zeros(self.np*nbins).reshape(self.np,nbins)
        bins=numpy.zeros(self.np*nbins).reshape(self.np,nbins)
        av=numpy.zeros(self.np); s2=numpy.zeros(self.np)
        frame=iterator(); iframe=0; ic=0
        while frame.any():
            frame=sphav*frame.sum(axis=1)#q=(qx+qy+qz)/srqt(3)
            av+=frame; s2+=frame*frame; chunk[ic]=frame
            frame=iterator(); iframe+=1; ic+=1;
            if not iframe%1000: print iframe
            if ic==cs: #update histograms with chunk
                chunk=chunk.transpose()
                for i in range(self.np): #update histograms
                    hs,bins[i]=histogram(chunk[i],bins=nbins,
                                         rg=(qmM[i][0],qmM[i][1]))
                    hist[i]+=hs
                chunk=chunk.transpose(); ic=0   #reset counter
        if ic: #there are some data leftover
            chunk=chunk.transpose()
            for i in range(self.np): #update histograms
                hs,bins[i]=histogram(chunk[i][0:ic],bins=nbins,
                                     rg=(qmM[i][0],qmM[i][1]))
                hist[i]+=hs
        #Calculate Chi-square test of each mode to theoretical
        #distribution 1/sqrt(2*s2) exp(-(q-avq)^2/(2*s2))
        av/=iframe; s2=s2/iframe-av*av; hist/=iframe; b=1/(2*s2)
        theor=numpy.exp(-b[:,numpy.newaxis]*(bins-av[:,numpy.newaxis])**2)
        theor/=theor.sum(axis=1)[:,numpy.newaxis]
        self.chisq=((hist-theor)**2/theor).sum(axis=1)



                
