import numpy,pdb
"""scattering of one-dimensional scatterers.

For instance, assume each scatterer is one quasiharmonic mode with a
single generalized coordinate

"""

def msd2(rr,delta,tl,nt):
    """time dependent mean square displacement 
    delta:  time in between frames
    rr=[frame_0,...,frame_N] 
        frame_i=[q_0,...,q_M] #frame at time 'i'
                 q_j #generalized coord. of scatterer 'j' at time 'i'
    tl: limit t. Compute msd2(t) from t=0 up to t=tl
    nt: compute nt different msd2 values
    returns array of shape=(len(rr[0]),nt)
    """
    nframe=len(rr)  #;print 'nframe=',nframe #number of frames
    nq=len(rr[0])   #;print 'nq=',nq #number of scatterers
    #x: ratio between time-between-frames and time-between-time-points
    norig=1000 #number of different origins we choose to average
    scattering=numpy.zeros( (nt,nq) )
    
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
            scattering[it]+=y*y
            iorig+=1
            iframe+=origjump #next origin frame
        scattering[it]=scattering[it]/iorig
        it+=1 #;print iorig
    scattering=scattering.transpose(1,0)
    return scattering #return scattering[nq][nt]
