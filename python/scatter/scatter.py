#!/usr/bin/python

import os,sys,pdb
import numpy
from utilities.small_utilities import Bye
from copy import deepcopy

#bound incoherent scattering length (http://www.ncnr.nist.gov/resources/n-lengths)
bisl={'H':25.274, '2H':4.04, 'C':0.0, 'N':2.0, 'O':0.0, 'Fe':0.0, 'S':0.0}
def reportScatteringLenthgs(names,masses):
    """ return two numpy arrays containing coherent and incoherent
    scattering lengths"""
    b = { 'coherent':[], 'incoherent':[] }
    for i in range( len(names) ):
        name=names[i][0]
        mass=masses[i]
        if name=='H':
            bi=25.274; bc=-3.7406
            if mass > 1.1:
                bi=4.04; bc=6.671  # deuterium
            if mass > 2.2:
                bi=-1.04; bc=4.792 # tritium
        elif name=='C':
            bi=0.0; bc=6.6511
            if mass > 12.1:
                bi=-0.52; bc=6.19 # C13
            if mass > 13.1:
                bi=0.00; bc=6.19  # C14
        elif name=='O':
            bi=0.0; bc=5.803
            if mass > 16.0:
                bi=0.18; bc=5.78  # O17
            if mass > 17.0:
                bi=0.00; bc=5.84  # O18
        elif 'N':
            bi=2.0; bc=9.37
            if mass > 14.1:
                bi=-0.02; bc=6.44 # N15
        elif 'S':
            bi=0.0; bc=2.804
            if mass > 32.0:
                bi=-0.02; bc=4.74 # S33
            if mass > 33.0:
                bi=0.00; bc=3.48  # S34 S36
        else:
            bi=0.0; bc=0.0
        b['coherent'].append( bc)
        b['incoherent'].append( bi)
    b['coherent'] = numpy.array( b['coherent'] )
    b['incoherent'] = numpy.array( b['incoherent'] )
    return b

def myfourier(ss,conv=numpy.array(None)):
    """this actually gives the same results as numpy.fft.fft"""
    from copy import deepcopy
    tt=deepcopy(numpy.array(ss)); nt=len(tt)
    if isinstance(conv,list): conv=numpy.array(conv)
    if not conv.any(): conv=numpy.ones(nt)
    tt=tt*conv
    a=numpy.arange(nt)
    a=(2*numpy.pi/nt) * (a*a.reshape(nt,1))
    cos=numpy.cos(a)
    sin=numpy.sin(a)
    return numpy.dot(tt,cos), numpy.dot(tt,sin)

def fourier(ss,conv=numpy.array(None)):
    from copy import deepcopy
    #make sure ss and conv are numpy arrays
    tt=deepcopy(numpy.array(ss))
    if isinstance(conv,list): conv=numpy.array(conv)
    if not conv.any(): conv=numpy.ones(len(tt))
    tt=tt*conv
    return numpy.fft.fft( tt )

def fourierRealSymmetric(ss,conv=numpy.array(None)):
    from copy import deepcopy
    #make sure ss and conv are numpy arrays
    tt=deepcopy(numpy.array(ss))
    if isinstance(conv,list): conv=numpy.array(conv)
    if not conv.any(): conv=numpy.ones(len(tt))
    
    tt=tt*conv
    ff = numpy.abs( numpy.fft.fft( tt ).real )
    return numpy.fft.fftshift( ff )

def fourier_abs(ss):
    """modulus of the fourier transform"""
    if type(ss)==type([]): ss=numpy.array(ss)
    ww=numpy.fft.fft(ss) ; ww=abs(ww)
    return ww

def EISF(qq,rr):
    from scipy import random
    qlist=[]
    if type(qq)==type(0.0):  qlist.append(qq)
    elif type(qq)==type([]): qlist=qq
    nq=len(qlist)
    scattering=[]
    nframe=len(rr)  #;print 'nframe=',nframe #number of frames
    NreadFrames=10000
    nxyz=len(rr[0]) #;print 'nxyz=',nxyz #number of scatterers
    #x: ratio between time-between-frames and time-between-time-points
    zeroarray=1E-300+numpy.zeros(nxyz,dtype=float) #avoid dividing by zero
    readFrames=0
    sumI=numpy.zeros(nq*nxyz,dtype=float).reshape(nq,nxyz)
    while readFrames<NreadFrames:
        frame_0=int(random.random()*nframe)
        frame_1=frame_0+int(random.standard_normal()*250)
        while frame_1<0 or frame_1>=nframe:
            frame_1=frame_0+int(random.standard_normal()*250)
        frame_0=rr[frame_0]
        frame_1=rr[frame_1]
        y=frame_1-frame_0
        y=numpy.sqrt( (y*y).sum(axis=1) )#vector of nxyz components
        for iq in range(nq):
            z=qlist[iq]*y+zeroarray
            sumI[iq]+=numpy.sin(z)/z
        readFrames+=1
    sumI=sumI/NreadFrames #; sumI=sumI*sumI
    for iq in range(nq):
        scattering.append( (qlist[iq],sum(sumI[iq])/nxyz) )
    return scattering

def EISF_bak(qq,rr):
    """elasting incoherent structure function
    qq: modulus of momenta exchange, in Angstrom^(-1). Can be a list
    rr=[frame_0,...,frame_N] 
        frame_i=[xyz_0,...,xyz_M] #frame at time 'i'
    EISF(q)=
     (1/3N)*
     \sum_{i=1}^N(       #sum over all scatterers
       \sum_{j=1}^N_fr(       #sum over all frames
         \integral[ d\Omega_q
           exp(i (\vec{q}(\vec{R}_i(j) - \vec{R}_i(frameN)) ))
         )
       )
      )
     )
    The integral in \Omega_q is analytically simple to solve:
    \integral[ d\Omega_q exp(i\vec{q}\vec{z}) ] = 4\pi numpy.sin(qz)/(qz)
    Thus:
    EISF(q)=
     (1/3N)*
     \sum_{j=1}^N_fr(       #sum over all frames
      \sum_{i=1}^N( numpy.sin(qz_i(j))/(qz_i(frameN)) )
     ),
    where z_i(j)=\norm(\vec{R}_i(j) - \vec{R}_i(frameN))
    """
    qlist=[]
    if type(qq)==type(0.0):  qlist.append(qq)
    elif type(qq)==type([]): qlist=qq
    nq=len(qlist)
    scattering=[]
    nframe=len(rr)  #;print 'nframe=',nframe #number of frames
    NreadFrames=10000
    nxyz=len(rr[0]) #;print 'nxyz=',nxyz #number of scatterers
    #x: ratio between time-between-frames and time-between-time-points
    zeroarray=1E-300+numpy.zeros(nxyz,dtype=float) #avoid dividing by zero
    readFrames=0
    sumI=numpy.zeros(nq*nxyz,dtype=float).reshape(nq,nxyz)
    while readFrames<NreadFrames:
        frame_0=rr[int(random()*nframe)]
        frame_1=rr[int(random()*nframe)]
        y=frame_1-frame_0
        y=numpy.sqrt( (y*y).sum(axis=1) )#vector of nxyz components
        for iq in range(nq):
            z=qlist[iq]*y+zeroarray
            sumI[iq]+=numpy.sin(z)/z
        readFrames+=1
    sumI=sumI/NreadFrames #; sumI=sumI*sumI
    for iq in range(nq):
        scattering.append( (qlist[iq],sum(sumI[iq])/nxyz) )
    return scattering

def EISF2(qq,rr):
    """elastic incoherent scattering function calculated according to:
    \vec{dr_i}=\vec{r_i}-<\vec{r_i}>_t,
    1/(4piN)\sum_{i=1}^N ( <  e^{i\vec{q}\vec{dr_i}}  >_{t,\omega_q} )^2=
    1/N\sum_{i=1}^N ( <  numpy.sin(q dr_i)/(q dr_i)  >_t )^2
    """
    nxyz=len(rr[0]) #;print 'nxyz=',nxyz #number of scatterers
    zeroarray=1E-300+numpy.zeros(nxyz,dtype=float) #avoid dividing by zero
    nframe=len(rr) #xyz should be equal to 3
    avframe=rr.sum(axis=0)/nframe #average coordinates
    qlist=[]
    if type(qq)==type(0.0):  qlist.append(qq)
    elif type(qq)==type([]): qlist=qq
    nq=len(qlist)
    scattering=[]
    sumI=numpy.zeros(nq*nxyz,dtype=float).reshape(nq,nxyz)
    for frame in rr:
        y=frame-avframe; y=numpy.sqrt( (y*y).sum(axis=1) )
        for iq in range(nq):
            z=qlist[iq]*y+zeroarray
            sumI[iq]+=numpy.sin(z)/z
    sumI=sumI/nframe
    for iq in range(nq):
        scattering.append( (qlist[iq],sum(sumI[iq])/nxyz) )
    return scattering

    
def msd2_from_EISF2(qq,rr):
    """obtain means square displacements in the gaussian approximation
    S=Fourier(EISF) #elastic incoherent structure factor
    In the gaussian approximation (<r^2>(T)*Q^2 < 2):
    |S|^2(Q,T)=|S_0|^2(T)*e^{-<r^2>(T)*Q^2/3} 
    """
    from scipy.stats import linregress
    scattering=EISF2(qq,rr);    I=[];
    for (q,i) in scattering: I.append(i)
    S=fourier_abs(I)
    a_s,b_s,r,tt,stderr=linregress(qq,S)
    msd2=a_s/6.0
    #second pass restricted to q's such that q^2*msd2<2
    j=0
    while j<len(qq) and qq[j]*qq[j]*msd2 < 2: j+=1
    a_s,b_s,r,tt,stderr=linregress(qq[0:j],S[0:j])
    return a_s/6.0

    
def msd2(rr,delta,tl,nt,norig=1000):
    """average mean square displacement 
    delta:  time in between frames
    rr=[frame_0,...,frame_N] 
        frame_i=[xyz_0,...,xyz_M] #frame at time 'i'
                 xyz_j=[x_j,y_j,z_j] #coords scatterer 'j' at time 'i'
    tl: limit t. Compute SIIScF(q,t) from t=0 up to t=tl
    nt: compute nt different SIIScF(q,t) values
    norig: number of different time-origins we choose to average over
    
    scale: how are the time-points selected ?
          'linear': t_i=tl*i/nt
          'log10' : t_i selected from log_10 scale
    SIIScF: (1D-array)
    Required that len(t)==len(rr)
    SIIScF(q,t)=
     (1/N)*
     \sum_{i=1}^N(
      \integral_0^{tmax-t}( ds
        \norm(\vec{R}_i(t+s) - \vec{R}_i(s))
      )
     )
    """
    nframe=len(rr)  #;print 'nframe=',nframe #number of frames
    nxyz=len(rr[0]) #;print 'nxyz=',nxyz #number of scatterers
    #x: ratio between time-between-frames and time-between-time-points
    scattering=numpy.zeros( (nt,nxyz) )
    
    shifts=[]
    dt=float(tl)/nt #time in between time points
    tratio=dt/delta
    for it in range(nt): shifts.append( int( it*tratio ))

    #print 'shifts=',shifts    
    it=1

    while it < nt:
        shift=shifts[it] ; iframe=0 ; iorig=0
        origjump=max(1,int( float(nframe-shift)/norig ))
        #print 'it=',it,' shift=',shift,' origjump=',origjump
        while iframe+shift<nframe: #time origin shift
            y=rr[iframe+shift]-rr[iframe]
            scattering[it]+=(y*y).sum(axis=1)
            iorig+=1
            iframe+=origjump #next origin frame
        scattering[it]=scattering[it]/iorig
        it+=1 #;print iorig
    #pdb.set_trace()
    scattering=scattering.transpose(1,0)
    return scattering #return scattering[nxyz][nt]

def autoCorr(rr,delta,tl,nt):
    """auto correlation function
    delta:  time in between frames
    rr=[frame_0,...,frame_N] 
        frame_i=[xyz_0,...,xyz_M] #frame at time 'i'
                 xyz_j=[x_j,y_j,z_j] #coords scatterer 'j' at time 'i'
    tl: limit t.
    nt: compute nt different values
    scale: how are the time-points selected ?
          'linear': t_i=tl*i/nt
          'log10' : t_i selected from log_10 scale
    autoCorr(q,t)=
     (1/N)*
     \sum_{i=1}^N(
      1/(tmax-t) * \integral_0^{tmax-t}( ds
        \vec{R}_i(t+s) - \vec{R}_i(s)
      )
     )
    """

    nframe=len(rr)  #;print 'nframe=',nframe #number of frames
    nxyz=len(rr[0]) #;print 'nxyz=',nxyz #number of scatterers
    norig=1000 #number of different origins we choose to average
    scattering=numpy.zeros( (nt,nxyz) )
    
    shifts=[]
    dt=float(tl)/nt #time in between time points
    tratio=dt/delta
    for it in range(nt): shifts.append( int( it*tratio ))
    
    #print 'shifts=',shifts    
    it=0
    while it < nt:
        shift=shifts[it] ; iframe=0 ; iorig=0
        origjump=max(1,int( float(nframe-shift)/norig ))
        #print 'it=',it,' shift=',shift,' origjump=',origjump
        while iframe+shift<nframe: #time origin shift
            y=rr[iframe+shift]*rr[iframe]
            scattering[it]+=y.sum(axis=1)
            iorig+=1
            iframe+=origjump #next origin frame
        scattering[it]=scattering[it]/iorig
        it+=1 #;print iorig
    scattering=scattering.transpose(1,0)
    #normalize to 1 for time=0
    for ixyz in range(nxyz):
        scattering[ixyz]=scattering[ixyz]/scattering[ixyz][0]
    return scattering

def SICScF(qq,rr,b,delta,tl,nt):
    """Self-Intermediate coherent scattering Function (SICSsF)
    q: modulus of momenta exchange, in Angstrom^(-1). Can be a list
    delta:  time in between two consecutive frames
    rr=[frame_0,...,frame_N] 
        frame_i=[xyz_0,...,xyz_M] #frame at time 'i'
                 xyz_j=[x_j,y_j,z_j] #coords scatterer 'j' at time 'i'
    b: coherent scattering factors
    tl: limit t (in psec). Compute SIIScF(q,t) from t=0 up to t=tl
    nt: compute for nt time shifts
    SIIScF: (1D-array)
    Required that len(t)==len(rr)
    SIIScF(q,t)=
     (1/N^2)*
     \sum_{i=1}^N\sum_{j=1}^N(
      \integral_0^{tmax-t}( ds
       \integral[ d\Omega_q
        b_i * b_j * exp(i (\vec{q}(\vec{R}_i(s) - \vec{R}_j(s+t)) ))
       )
      )
     )
    The integral in \Omega_q is analytically simple to solve:
    \integral[ d\Omega_q exp(i\vec{q}\vec{z}) ] = 4\pi numpy.sin(qz)/(qz)
    Thus:
    SIIScF(q,t)=
     (1/N^2)*
     \integral_{0}^{tmax-t}( ds
      b_i * b_j * \sum_{i=1}^N( numpy.sin(qz_i(t+s))/(qz_i(t+s)) )
     ),
    where z_i(t+s)=\norm(\vec{R}_i(t+s) - \vec{R}_i(s))
    """
    #convert to numpy array
    qlist=[]
    if type(qq)==type(0.0):  qlist.append(qq)
    elif type(qq)==type([]): qlist=qq
    qlist=numpy.array(qlist)
    nq=len(qlist)
    
    nframe=len(rr)  #;print 'nframe=',nframe #number of frames
    nxyz=len(rr[0]) #;print 'nxyz=',nxyz #number of scatterers
    scattering=numpy.zeros( (nt,nq,nxyz) )
    scattering[0]+=1
    norig=100 #number of time origins we choose to average over
    shifts=[]
    dt=float(tl)/nt #time in between time points
    tratio=dt/delta
    for it in range(nt): shifts.append( int( it*tratio ))
    zeroarray=1E-300+numpy.zeros(nxyz,dtype=float) #avoid dividing by zero
    it=1
    while it < nt: #iterate over the time points
        shift=shifts[it] ; iframe=0 ; iorig=0
        origjump=max(1,int( float(nframe-shift)/norig ))
        #print 'nframe=',nframe,' norig=',norig
        #print 'it=',it,' shift=',shift,' origjump=',origjump
        while iframe+shift<nframe:
            rr=rr[iframe+shift]-r1[iframe]
            y=numpy.sqrt( (y*y).sum(axis=1) )
            z=(y+zeroarray)*qlist.reshape((nq,1))
            scattering[it]+=numpy.sin(z)/z
            iorig+=1
            iframe+=origjump #next origin frame
        scattering[it]/=iorig
        it+=1 #;print iorig
    scattering=scattering.transpose(2,1,0)

def SIIScF(qq,rr,delta,tl,nt):
    """Self-Intermediate Incoherent Scattering Function (SIIScF)
    q: modulus of momenta exchange, in Angstrom^(-1). Can be a list
    delta:  time in between two consecutive frames
    rr=[frame_0,...,frame_N] 
        frame_i=[xyz_0,...,xyz_M] #frame at time 'i'
                 xyz_j=[x_j,y_j,z_j] #coords scatterer 'j' at time 'i'
    tl: limit t (in psec). Compute SIIScF(q,t) from t=0 up to t=tl
    nt: compute for nt time shifts, number of time points in which to divide tl
    SIIScF: (1D-array)
    Required that len(t)==len(rr)
    SIIScF(q,t)=
     (1/3N)*
     \sum_{i=1}^N(
      \integral_0^{tmax-t}( ds
       \integral[ d\Omega_q
        exp(i (\vec{q}(\vec{R}_i(t+s) - \vec{R}_i(s)) ))
       )
      )
     )
    The integral in \Omega_q is analytically simple to solve:
    \integral[ d\Omega_q exp(i\vec{q}\vec{z}) ] = 4\pi numpy.sin(qz)/(qz)
    Thus:
    SIIScF(q,t)=
     (1/N)*
     \integral_{0}^{tmax-t}( ds
      \sum_{i=1}^N( numpy.sin(qz_i(t+s))/(qz_i(t+s)) )
     ),
    where z_i(t+s)=\norm(\vec{R}_i(t+s) - \vec{R}_i(s))
    """
    #convert to numpy array
    qlist=[]
    if type(qq)==type(0.0):  qlist.append(qq)
    elif type(qq)==type([]): qlist=qq
    qlist=numpy.array(qlist)
    nq=len(qlist)
    nframe=len(rr)  #;print 'nframe=',nframe #number of frames
    nxyz=len(rr[0]) #;print 'nxyz=',nxyz #number of scatterers
    scattering=numpy.zeros( (nt,nq,nxyz) )
    scattering[0]+=1
    norig=1000 #number of time origins we choose to average over
    shifts=[]
    dt=float(tl)/nt #time in between time points
    tratio=dt/delta
    for it in range(nt): shifts.append( int( it*tratio ))
    zeroarray=1E-300+numpy.zeros(nxyz,dtype=float) #avoid dividing by zero
    it=1
    while it < nt: #iterate over the time points
        shift=shifts[it] ; iframe=0 ; iorig=0
        origjump=max(1,int( float(nframe-shift)/norig ))
        #print 'nframe=',nframe,' norig=',norig
        #print 'it=',it,' shift=',shift,' origjump=',origjump
        while iframe+shift<nframe:
            y=rr[iframe+shift]-rr[iframe] #substract element by element
            y=numpy.sqrt( (y*y).sum(axis=1) )
            z=(y+zeroarray)*qlist.reshape((nq,1))
            scattering[it]+=numpy.sin(z)/z
            iorig+=1
            iframe+=origjump #next origin frame
        scattering[it]/=iorig
        it+=1 #;print iorig
    scattering=scattering.transpose(2,1,0)
    return scattering #return scattering[nxyz][nq][nt]

def intSIIScF(qq,rr,delta,tl,nt):

    """Integral of the Self-Intermediate Incoherent Scattering
    Function
    
    q: modulus of momenta exchange, in Angstrom^(-1). Can be a list
    delta:  time in between two consecutive frames
    rr=[frame_0,...,frame_N] 
        frame_i=[xyz_0,...,xyz_M] #frame at time 'i'
                 xyz_j=[x_j,y_j,z_j] #coords scatterer 'j' at time 'i'
    tl: limit t (in psec). Compute SIIScF(q,t) from t=0 up to t=tl
    nt: compute for nt time shifts
    SIIScF: (1D-array)
    Required that len(t)==len(rr)

    The Self-Intermediate Incoherent Scattering Function
    SIIScF(q,t)=
     (1/3N)*
     \sum_{i=1}^N(
      \integral_0^{tmax-t}( ds
       \integral[ d\Omega_q
        exp(i (\vec{q}(\vec{R}_i(t+s) - \vec{R}_i(s)) ))
       )
      )
     )
     
    The integral in \Omega_q is analytically simple to solve:
    \integral[ d\Omega_q exp(i\vec{q}\vec{z}) ] = 4\pi sin(qz)/(qz)
    The spatial integral in q is also analytically simple to solve:
    \integral[dq q^2 sin(qz)/(qz)]=sin(qz)/z^3-q*cos(qz)/z^2

    Thus:
    forward(q,t)=
     (1/N)*
     \integral_{0}^{tmax-t}( ds
      \sum_{i=1}^N( numpy.sin(qz_i(t+s))/(z_i(t+s)^3) )
     ),
    backward(q,t)=
     (1/N)*
     \integral_{0}^{tmax-t}( ds
      \sum_{i=1}^N( numpy.cos(qz_i(t+s))/(z_i(t+s)^2) )
     ),
    where z_i(t+s)=\norm(\vec{R}_i(t+s) - \vec{R}_i(s))

    Note we calculate the indefinite integral. Later we will use the
    results to calculate the integral between two particular q-values
    scattering[nxyz][q_{final}][nt]-scattering[nxyz][q_{initial}][nt]
    """
    #convert to numpy array
    qlist=[]
    if type(qq)==type(0.0):  qlist.append(qq)
    elif type(qq)==type([]): qlist=qq
    nq=len(qlist); qlist=(numpy.array(qlist)).reshape((nq,1))
    nframe=len(rr)  #;print 'nframe=',nframe #number of frames
    nxyz=len(rr[0]) #;print 'nxyz=',nxyz #number of scatterers
    scattering=numpy.zeros( (nt,nq,nxyz) )
    norig=1000 #number of time origins we choose to average over
    shifts=[]
    dt=float(tl)/nt #time in between time points
    tratio=dt/delta
    for it in range(nt): shifts.append( int( it*tratio ))
    zeroarray=1E-8+numpy.zeros(nxyz,dtype=float) #avoid dividing by zero
    scattering[0]=numpy.ones( (nq,nxyz) )
    it=1
    while it < nt: #iterate over the time points
        shift=shifts[it] ; iframe=0 ; iorig=0
        origjump=max(1,int( float(nframe-shift)/norig ))
        #print 'nframe=',nframe,' norig=',norig
        #print 'it=',it,' shift=',shift,' origjump=',origjump
        while iframe+shift<nframe:            
            x=rr[iframe+shift]-rr[iframe] #substract element by element
            y=zeroarray+(x*x).sum(axis=1); z2=1.0/y
            y=zeroarray+numpy.sqrt(y); z=1.0/y; y=qlist*y
            scattering[it]+=(z*z2)*numpy.sin(y)-(qlist*z2)*numpy.cos(y)
            iorig+=1
            iframe+=origjump #next origin frame
        scattering[it]/=iorig;
        it+=1 #;print iorig
    scattering=scattering.transpose(2,1,0)
    return scattering #return scattering[nxyz][nq][nt]

def tt2ww(tt):
    """find energies (in meV) assuming times in picoseconds"""
    return numpy.arange(len(tt))*4.136/(tt[-1]-tt[0])

def storeSIIScF(outfile,scattering,qq,tt):
    """scattering[ len(qq) ][ len(tt) ]
    File format:
    
    #    
    t1 isf(t1,q1) ... isf(t1,qN)
    .
    .
    tM isf(tM,q1) ... isf(tM,qM)

    tt: list of time points. Assume in picoseconds
    """
    
    buf='# t(psec) or w(meV)\n'
    buf+='#'
    for q in qq: buf+=' %8.4f'%(q,)
    buf+='\n'

    for itt in range(len(tt)):
        buf+='%7.3f'%(tt[itt],)
        for iq in range(len(qq)):
            buf+=' %8.6f'%(scattering[iq][itt],)
        buf+='\n'
    open(outfile,'w').write(buf)

def store_intSIIScF(outfile,scattering,tt):
    """scattering[ len(qq) ][ len(tt) ]
    File format:
    
    #    
    t1 isf(t1)
    .
    .
    tM isf(tM)

    tt: list of time points. Assume in picoseconds
    """
    
    buf='# t(psec) or w(meV)\n'

    for itt in range(len(tt)):
        buf+='%7.3f %8.6f\n'%(tt[itt],scattering[itt])
    open(outfile,'w').write(buf)

def FourierSIIScf(scattering,tt=None,resolFun=None,transform=fourier):
    """
    Fourier transform the intermediate scattering function SIIScf
    scattering[dim1][dim2]..[Ntt]  Ntt: number of time points
    errorFun: multiple SIIScf by resolution function in time domain
    """    
    shape=scattering.shape; Ntt=shape[-1]; rf=numpy.ones(Ntt)
    if resolFun: #correct by resolution function
        A=1.0/resolFun(0.0)
        assert tt , 'FourierSIIScf needs tt'
        for itt in range(len(tt)): rf[itt]=A*resolFun(tt[itt])
    sQW=deepcopy(scattering).reshape(numpy.product(shape[:-1]),Ntt)
    for i in range(len(sQW)):
        sQW[i]=transform( sQW[i],conv=rf )
    sQW=sQW.reshape( shape )
    return sQW

def store_sQW(outfile,scattering,qq,ww):
    """scattering[ len(qq) ][ len(tt) ]
    File format:
    
    #    
    w1 sqw(w1,q1) ... sqw(w1,qN)
    .
    .
    wM sqw(wM,q1) ... sqw(wM,qM)
    
    """
    buf='#'+' '*6
    for q in qq: buf+=' %8.4f'%(q,)
    buf+='\n'
    for iww in range(len(ww)):
        buf+='%7.3f'%(ww[iww],)
        for iq in range(len(qq)):
            buf+=' %8.6f'%(scattering[iq][iww],)
        buf+='\n'
    open(outfile,'w').write(buf)
