import numpy,pdb
from math import pi

co = 0.001 #cut-off for small values

def write_from_numpy(pt,arr,format=' %f',ncols=7,comment=''):
    """write numpy array into file"""
    buf=''
    if comment: buf+=comment+'\n'
    b=arr.ravel(); L=len(b) #reshape to one dimension
    index=0
    while True:
        for jndex in range(ncols):
            if index==L: break
            buf+=format%(b[index]); index+=1
        buf+='\n'
        if index==L: break
    if pt:
        pt.write(buf)
        pt.flush()
    return buf

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

def scatL(U,group):
    """report scattering lenghts
       U: universe
       group: group of atoms"""
    X = U.selectAtoms(group)
    return reportScatteringLenthgs(X.names(), X.masses() )
    
    
def sinappx(x):
    """approximation in the third order to the sin(x) when x is a
    numpy arrray of positive numbers"""
    twopi = 2 * pi    #2pi
    piH   = 0.5 * pi  #pi/2
    minuspiH = -piH   #-pi/2
    piX   = 1.5 * pi  #3pi/2
    
    z = x % twopi  # 0 < z < 2pi
    #shift Z by 2pi if piX < z, since sin(z-2pi)=sin(z) in this range
    z[ numpy.where(z > piX ) ] -= twopi  #z-->z-2pi, piX < z
    #change sign of z and shift by pi if piH < z, since
    # sin(z) = sin(pi-z) in this range
    z[ numpy.where(z > piH ) ] *= 1      #z-->-z, piH < z
    z[ numpy.where(z < minuspiH ) ] += pi  #z-->pi-z, piH < z

    #we have changed the 0<z<2pi range to -pi/2<z<pi/2 range while
    #preserving the value of sin(z), and now we use a Taylor expansion
    #to the third order in sin(z) for the -pi/2<z<pi/2 range. The
    #maximum error is 0.004 for abs(z)=pi/2
    
    z2 = z * z
    return z * ( 1 - 0.166666666 * z2 * ( 1 - 0.02 * z2 ) )

def Qlist(Qmin,Qmax,deltaQ):
    """generate list of Q-values"""
    QQ=[]; Q = Qmin
    while Q<=Qmax:
        QQ.append( Q ); Q += deltaQ
    return numpy.array(QQ)


def getFrames(U,group,skip=1,format='fac'):
    """return frames
       U: universe
       group: group of atoms
       skip: jump between frames
       format: shape of frames array, from permutting 'f', 'a', 'c'
               which are abbreviations for frames, atoms, coordinates.
               Example: frames[#frames][#atoms][3] for 'fac'"""
    X = U.selectAtoms(group)
    frames = U.trajectory.timeseries(X,skip=skip,format=format)
    frames = numpy.array(frames,dtype='float32')
    return frames

def SSF(frames,b,QQ,nsampling=100,
        comment='# static structure factor\n'):
    """static structure factor
    frames: numpy.array with shape=(#frames,3,#atoms)
    b: numpy.array of coherent scattering lenghts
    QQ: numpy.array list of Q values
    nsampling: number of t0

    d_{i,j} = |r_i(t0)-r_j(t0)|
    sum_{i,j} b_i * b_j * sin(Q*d_{i,j})/(Q*d_{i,j})
    """
    from MDAnalysis.analysis.distances import distance_array

    nQ = len(QQ)
    nat = len(b)
    nfr = len(frames)

    buf= comment
    buf+='# Q (A^(-1)) struc-fact\n'

    b_ij = numpy.outer(b,b)
    s_Q = numpy.zeros(nQ)
    delta = max(1,int(nfr/nsampling ))#separat. between t0's
    t0 = 1  #initialize the time point
    isp = 0
    while t0 < nfr: #cycle over all t0's
        fi = frames[t0]
        dij=distance_array(fi,fi,None) #distances between atoms
        for iQ in range(nQ):
            arg = QQ[iQ]*dij
            arg[numpy.where(arg<co)]=co #avoid sin(arg)/arg trouble when arg goes to zero
            s_Q[iQ] += ( b_ij * numpy.sin(arg)/arg ).sum()/(nat*nat)
        isp += 1
        t0 += delta #jump to next t0
    s_Q = s_Q / isp    
    for i in range(nQ):
        buf+= ' %6.3f %10.7f\n'%(QQ[i],s_Q[i])

    return {'buf':buf,}

def SSF_xself(frames,b,QQ,nsampling=100,
        comment='# static structure factor\n'):
    """static structure factor excluded self and normalized by
    the total incoherent scattering lenght, squared

    frames: numpy.array with shape=(#frames,3,#atoms)
    b: numpy.array of coherent scattering lenghts
    QQ: numpy.array list of Q values
    nsampling: number of t0

    d_{i,j} = |r_i(t0)-r_j(t0)|
    sum_{i,j} b_i * b_j * sin(Q*d_{i,j})/(Q*d_{i,j})
    """
    from MDAnalysis.analysis.distances import distance_array

    nQ = len(QQ)
    nat = len(b)
    nfr = len(frames)

    buf= comment
    buf+='# Q (A^(-1)) struc-fact\n'

    b_ij = numpy.outer(b,b)
    bbcoh=(b*b).sum()*(numpy.sin(co)/co) #total self coherent term
    bbinc = b.sum()*b.sum()        #total incoherent term, squared
    s_Q = numpy.zeros(nQ)
    delta = max(1,int(nfr/nsampling ))#separat. between t0's
    t0 = 1  #initialize the time point
    isp = 0
    while t0 < nfr: #cycle over all t0's
        fi = frames[t0]
        dij=distance_array(fi,fi,None) #distances between atoms, exclude self
        for iQ in range(nQ):
            arg = QQ[iQ]*dij
            arg[numpy.where(arg<co)]=co #avoid sin(arg)/arg trouble when arg goes to zero
            s_Q[iQ] += ( (b_ij * numpy.sin(arg)/arg).sum() - bbcoh)/bbinc #remove the self coherent term
        isp += 1
        t0 += delta #jump to next t0
    s_Q = s_Q / isp #average by time points
    for i in range(nQ):
        buf+= ' %6.3f %10.7f\n'%(QQ[i],s_Q[i])

    return {'buf':buf,}

   
def SSFcheck(frames,b,QQ,nsampling=100,
        comment='# static structure factor\n'):
    """print the static structure factor at different t0's to check it
    is invariant with respect to t0.
    frames: numpy.array with shape=(#frames,3,#atoms)
    b: numpy.array of coherent scattering lenghts
    QQ: numpy.array list of Q values
    nsampling: number of t0

    d_{i,j} = |r_i(t0)-r_j(t0)|
    sum_{i,j} b_i * b_j * sin(Q*d_{i,j})/(Q*d_{i,j})
    """
    from MDAnalysis.analysis.distances import distance_array

    nQ = len(QQ)
    nat = len(b)
    nfr = len(frames)

    buf= comment
    buf += '#Q-values (A^(-1)): '
    for Q in QQ: buf += ' %4.2f'%Q
    buf +='\n'
    buf+='# t0 stat-struc-fact for different Q-values\n'

    w_ij = numpy.zeros(nQ*nat*nat).reshape(nQ,nat,nat)
    delta = max(1,int(nfr/nsampling ))#separat. between t0's
    t0 = 1  #initialize the time point
    while t0 < nfr: #cycle over all t0's
        fi = frames[t0]
        dij=distance_array(fi,fi,None) #distances between atoms
        for iQ in range(nQ):
            arg = QQ[iQ]*dij
            arg[numpy.where(arg<co)]=co#avoid sin(arg)/arg troubl
            w_ij[iQ] = numpy.sin(arg)/arg
        #print the structure factor
        w_ij *= numpy.outer(b,b)
        z_ij = w_ij.sum(axis=2).sum(axis=1)/(nat*nat)
        buf += ' %5d'%t0
        for iQ in range(nQ): buf += ' %10.7f'%z_ij[iQ]
        buf += '\n'
        t0 += delta #jump to next t0
    return {'buf':buf,}

    
def IC(frames,b,QQ,nsampling=10,nt=100,dt=1,
       comment='# intermediate coherent scattering factor\n'):
    """Intermediate coherent scattering function
    frames: numpy.array with shape=(#frames,3,#atoms)
    b: numpy.array of coherent scattering lenghts
    QQ: numpy.array list of Q values
    nsampling: number of t0's for a given t
    nt: number of t's
    dt: spacing between consecutive t's, in number of frames

    d_{i,j} = |r_i(t0)-r_j(t0+t)|
    IC = sum_{i,j} b_i * b_j * sin(Q*d_{i,j})/(Q*d_{i,j}) / len(b)^2"""

    from MDAnalysis.analysis.distances import distance_array

    nQ = len(QQ)
    nat = len(b)
    nfr = len(frames)
    
    it = 0 #current number of evaluated t's
    t = 0  #time difference between t0 and t0+t
    
    buf= comment
    buf+='# Q-values=%s' % ' '.join( [str(Q) for Q in QQ] )+'\n'
    buf+='# time(psec) Q-values(A^(-1))\n'

    sf = numpy.zeros((1+nt)*nQ).reshape(1+nt,nQ)
    ts=[]
    while it < 1+nt: #cycle over all t's (+1 for t=0)
        s_ij = numpy.zeros(nQ*nat*nat).reshape(nQ,nat,nat)
        delta = max(1,int((nfr-t)/nsampling ))#separat. between t0's
        isp = 0 #number of sampled time points
        t0 = 1  #initialize the time point
        print 'it=%d t=%d'%(it,t)
        while t0+t < nfr: #cycle over all t0's
            #print 'it=%d t=%d isp=%d t0=%d'%(it,t,isp,t0)
            fi = frames[t0]
            fj = frames[t0+t]
            dij=distance_array(fi,fj,None) #distances between atoms
            for iQ in range(nQ):
                arg = QQ[iQ]*dij
                arg[numpy.where(arg<co)]=co#avoid sin(arg)/arg troubl
                s_ij[iQ] += numpy.sin(arg)/arg
            isp += 1
            t0 += delta #jump to next t0
        s_ij *= numpy.outer(b,b) / isp
        s_ij = s_ij.sum(axis=2).sum(axis=1)/(nat*nat)
        sf[it] = s_ij
        ts.append(t)
        s = write_from_numpy(None,s_ij,ncols=nQ,format=' %10.7f')
        buf+=' %6.3f %s'%(t,s)
        it += 1
        t += dt
    return {'buf':buf,'sf':sf,'ts':numpy.array(ts)}


def II(frames,b,QQ,nsampling=10,nt=100,dt=1,
       comment='# intermediate incoherent scattering factor\n'):
    """Intermediate incoherent scattering function
       frames: numpy.array with shape=(#frames,#atoms,3)
       b: numpy.array of incoherent scattering lenghts
       QQ: numpy.array list of Q values
       nsampling: number of t0's for a given t
       nt: number of t's
       dt: spacing between consecutive t's, in number of frames

       d_{i,j} = |r_i(t0)-r_j(t0+t)|
       II = sum_i b_i * sin(Q*d_{i,j})/(Q*d_{i,j}) / len(b)"""

    nQ = len(QQ)
    nat = len(b)
    nfr = len(frames)
    
    it = 0 #current number of evaluated t's
    t = 0  #time difference between t0 and t0+t

    buf= comment
    buf+='# Q-values=%s' % ' '.join( [str(Q) for Q in QQ] )+'\n'
    buf+='# time(psec) Q-values(A^(-1))\n'

    sf = numpy.zeros((1+nt)*nQ).reshape(1+nt,nQ)
    ts=[]
    while it < 1+nt: #cycle over all t's (+1 for t=0)
        s_ii = numpy.zeros(nQ*nat).reshape(nQ,nat)
        delta = max(1,int((nfr-t)/nsampling ))#separat. between t0's
        isp = 0 #number of sampled time points
        t0 = 1  #initialize the time point
        while t0+t < nfr: #cycle over all t0's
            dii = frames[t0] - frames[t0+t]
            dii *= dii
            dii = numpy.sqrt( dii.sum(axis=1) )
            for iQ in range(nQ):
                arg = QQ[iQ]*dii
                arg[numpy.where(arg<co)]=co#avoid sin(arg)/arg troubl
                s_ii[iQ] += numpy.sin(arg)/arg
            isp += 1
            t0 += delta #jump to next t0
        s_ii *= b / isp
        s_ii = s_ii.sum(axis=1)/nat
        sf[it] = s_ii
        ts.append(t)       
        s = write_from_numpy(None,s_ii,ncols=nQ,format=' %10.7f')
        buf+=' %6.3f %s'%(t,s)
        it += 1
        t += dt
    return {'buf':buf,'sf':sf,'ts':numpy.array(ts)}
