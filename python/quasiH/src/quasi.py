import os
import numpy
import pickle
#from pdb import set_trace as trace # Uncomment only for debugging purposes

"""Quasi Harmonic Analysis. Two degrees of freedom interact via a spring
V_{i,j}=1/2 Sum_m Sum_l K_{i,j,l,m} ( r_{i,l}-r_{j,m}-<r_{i,l}-r_{j,m}> ),
        where r_{i,l} is component 'l' of position vector r_i, and
        K_{i,j,l,m} is a 3x3 spring tensor for a given {i,j} pair, that is,
        a K_{i,j,l,m} has a total of 3Nx3N components
In our model, the interparticle interaction is isotropic
"""
jk=os.path.join

class modes:

    def __init__(self,nd):
        """
        self.nd        #number of degrees of freedom
        self.M         #mass of each degrees of freedom
        self.sM        #sqrt(M)
        self.isM       #1/sqrt(M)
        self.varcov    #mass weighted variance-covariance matrix
        self.evals     #eigenvalues of varcov, ordered from big to small
        self.evecs     #eigenvectors associated to evals
        self.FC        #force constant acting on coords (not mass weighted)
        self.eflux     #net flux between degrees of freedom
        self.chisq     #Chi-square values for each fit of the observed
                        distribution of values for each mode to a normal
                        distribution
        """
        self.nd=nd

    def pickleDump(self,pout):
        """save object with picke"""
        pickle.dump(self,pout,True)

    def pickleLoad(self,pin):
        """load object with pickle"""
        tmp=pickle.load(pin)
        self.__dict__.update(tmp.__dict__)
        
    def loadM(self,M): #load masses and precompute some
        self.M=M
        self.sM=numpy.sqrt(M)  #sqrt(M)
        self.isM=1/self.sM     #1/sqrt(M)
        
    def findVarCovar(self,iterator):
        """find variance-covariance matrix of mass-weigthed coordinates
        'iterator' is a function  providing a fresh coordinate snapshot
        in the shape of a numpy.array(nd)
        """
        assert hasattr(self,'M'), 'Run first loadM'
        avF=numpy.zeros(self.nd)
        stdF=numpy.zeros(self.nd*self.nd).reshape(self.nd,self.nd)
        #We require a frame.shape==(nd,)
        frame=iterator().ravel(); iframe=0
        while frame.any():
            stdF+=numpy.outer(frame,frame);avF+=frame
            frame=iterator().ravel(); iframe+=1; #load next frame
            if not iframe%1000: print iframe
            #print 'only for debugging'
            #if iframe==2000: break
            #end of only for debugging
        avF/=iframe; stdF/=iframe; self.avF=avF
        self.varcov= stdF - numpy.outer(avF,avF)
        self.varcov*=numpy.outer(self.sM,self.sM) #mass weight it!

    def findModes(self):
        """find normal modes given a variance-covariance matrix
        varcovar=E( (x-mu)^t (x-mu) )    (NOTE: ^t means transpose)
        evals,evecs are nd, (nd x nd) matrices, respectively"""
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
        We multiply by self.sM to remove the mass weighting
        sqrt(M)*x=evecs*q
        Last six force constant correspond to translations!"""
        if not hasattr(self,'evals'): self.findModes()
        self.FC=numpy.diag( 1/self.evals )
        self.FC[-6:]=0.0        #remove rotation and translation
        self.FC=numpy.dot(self.FC,self.evecs.transpose())
        self.FC=numpy.dot(self.evecs,self.FC)
        self.FC=kT*numpy.outer(self.sM,self.sM)*self.FC
            
    def projectTrajectory(self,iterator,pout):
        """project trajectory frames onto normal modes

         'iterator' is a function  providing a fresh coordinate snapshot
        in the shape of a numpy.array(self.nd)

        'pout' is a pointer to file where we store the frames in the
        normal mode system of reference
        """
        from helperfunc import write_from_numpy
        if not hasattr(self,'evals'): self.findModes()
        qframe=numpy.zeros(self.nd)
        frame=iterator(); iframe=0; buf=''; sizeMax=50E6 #50MB
        #q=evecs^t*sqrt(M)*x
        U=self.evecs.transpose()*self.sM[numpy.newaxis:,]
        while frame.any():
            qframe=numpy.dot(U,frame-self.avF)
            buf+=write_from_numpy(None,qframe,fmt=' %10.6f',ncols=7)
            if len(buf)>sizeMax: #flush the buffer
                pout.write(buf); buf=''
            frame=iterator(); iframe+=1; #load next frame
            if not iframe%1000: print iframe
        pout.write(buf) #flush remaining frames
        pout.close()

    def mode2cart(self,q):
        """from mode components q to cartesian coordinates"""
        return self.isM*numpy.dot(self.evecs,q)
    
    def eflux(self,xiterator,viterator,avF=numpy.zeros(None)):
        """find net energy flux between degrees of freedom
        'xiterator' returns coordinate frame
        'viterator' returns velocity frame
        flux=(1/4)*k_{ij}*(x_i*v_j-v_i*x_j), and  x_i:displacement vector
        """
        if not avF.any(): avF=self.avF
        if not hasattr(self,'FC'): self.findForceConstants()
        netflux=numpy.zeros(self.nd*self.nd).reshape(self.nd,self.nd)
        xframe=xiterator(); vframe=viterator(); iframe=0
        while xframe.any():
            xframe-=avF
            commuter=numpy.outer(xframe,vframe)
            commuter=commuter-commuter.transpose()
            netflux+=self.FC*commuter
            xframe=xiterator(); vframe=viterator(); iframe+=1
            if not iframe%1000: print iframe

            #print 'only for debugging'
            #if iframe==10000: break
            #end of only for debugging

        self.netflux=(1.0/4)*netflux
        return self.netflux

    def chisquare(self,iterator,rewind):
        """evaluate the distribution of mode components to the
        theoretical normal distribution.

        The theoretical probability distribution for a single node is
        given by:
        P(q_i)=[2\pi/k]^(1/2)*exp(-1/2*k*q_i^2)

        We have nd q_i variables, one for each mode. We will find the
        pvalue of their histograms with respect to the theoretical
        histograms.

        'iterator' provides us with a new frame
        'rewind' go to beginning of trajectory
        """
        from helperfunc import histogram
        
        #find range of q's for each mode
        qmM=numpy.zeros(self.nd*3).reshape(3,self.nd)
        frame=iterator().ravel(); iframe=0
        while frame.any():
            qmM[2]=frame; qmM[0]=qmM.min(axis=0); qmM[1]=qmM.max(axis=0)
            frame=iterator().ravel(); iframe+=1;
            if not iframe%1000: print iframe
        qmM=qmM.transpose()
        #read again the mode trajectory & generate histograms
        rewind();
        cs=1000 #read trajectory in chunks
        nbins=min(100,iframe/100) #number of bins per histogram
        chunk=numpy.zeros(cs*self.nd).reshape(cs,self.nd)
        hist=numpy.zeros(self.nd*nbins).reshape(self.nd,nbins)
        bins=numpy.zeros(self.nd*nbins).reshape(self.nd,nbins)
        av=numpy.zeros(self.nd); s2=numpy.zeros(self.nd)
        frame=iterator().ravel(); iframe=0; ic=0
        while frame.any():
            av+=frame; s2+=frame*frame; chunk[ic]=frame
            frame=iterator().ravel(); iframe+=1; ic+=1;
            if not iframe%1000: print iframe
            if ic==cs: #update histograms with chunk
                chunk=chunk.transpose()
                for i in range(self.nd): #update histograms
                    hs,bins[i]=histogram(chunk[i],bins=nbins,
                                         rg=(qmM[i][0],qmM[i][1]))
                    hist[i]+=hs
                chunk=chunk.transpose(); ic=0   #reset counter
        if ic: #there are some data leftover
            chunk=chunk.transpose()
            for i in range(self.nd): #update histograms
                hs,bins[i]=histogram(chunk[i][0:ic],bins=nbins,
                                     rg=(qmM[i][0],qmM[i][1]))
                hist[i]+=hs
        #Calculate Chi-square test of each mode to theoretical
        #distribution 1/sqrt(2*s2) exp(-(q-avq)^2/(2*s2))
        #remove those bins with population below 5% of that of the maximum
        av/=iframe; s2=s2/iframe-av*av; hist/=iframe; b=1/(2*s2)
        theor=numpy.exp(-b[:,numpy.newaxis]*(bins-av[:,numpy.newaxis])**2)
        theor/=theor.sum(axis=1)[:,numpy.newaxis]
        self.chisq=((hist-theor)**2/theor).sum(axis=1)
        return self.chisq,hist,bins

    def colorDispl(self,pdbf,rootdir,indexes=None):
        """for every mode, color as b-factor in the provided pdbf the
        displacement of each atom

        indexes: by default we assume that nd=natom*3, but if the
        modes are only calculate for hydrogens, then we need to know
        which atoms in the pdb file are the hydrogens. All others are
        given a zero displacement value"""
        from helperfunc import insertBfact
        #contract in first index. At begginning we have evecs.shape=(nd,nd)
        #      | e_11  .. e_1,nd |
        #evecs=| .          .    |
        #      |e_nd,1 .. e_nd,nd|
        U=self.evecs.transpose().reshape(self.nd,self.nd/3,3)
        U=(U*U).sum(axis=2).transpose()
        #now we have U.shape=(nd/3,nd)
        #  |u_11=(e_11^2+e_21^2+e_31^2) .. u_1nd=(e_1nd^2+e_2nd^2+e_3nd^2)|
        #U=|          .                           .                       |
        #  |u_{nd/3,1}=(e_{nd-2,1}^2+(e_{nd-1,1}^2+(e_{nd,1}^2) ..        |
        #
        #u_ij represents the square of the unit displacement for atom
        #'i' when mode 'j' is excited.
        U=numpy.sqrt(U).transpose()#now a row per mode, instead of a column
        if indexes:
            #guess how many atom in the pdb file
            nat=0; nlen=len(indexes)
            for l in open(pdbf).readlines():
                if l[0:5]=='ATOM ':
                    nat+=1
            for i in range(self.nd):
                bf=numpy.zeros(nat)
                for j in range(nlen): bf[ indexes[j] ]=U[i][j]
                buf=insertBfact(pdbf,bf)
                open(jk(rootdir,'%04d.pdb'%(i+1)),'w').write(buf)
        else:
            nat=(self.nd)/3
            for i in range(self.nd):
                buf=insertBfact(pdbf,U[i])
                open(jk(rootdir,'%04d.pdb'%(i+1)),'w').write(buf)

