import numpy as np
from math import pi
import h5py
from pdb import set_trace as tr

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
    scattering lengths
    Arguments:
      names: list specifying the atomic element name for each atom.
      masses: list specifying the mass of each atom. Discriminates isotope
    Returns:
      b: dictionary
         b['coherent']: 1D numpy array of coherent lenghts for each atom
         b['incoherent']: 1D numpy array of incoherent lenghts for each atom
"""
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
    b['coherent'] = np.array( b['coherent'] )
    b['incoherent'] = np.array( b['incoherent'] )
    return b


def scatL(group):
    """report scattering lenghts
       U: universe
       group: group of atoms"""
    return reportScatteringLenthgs(group.names, group.masses )


def sinappx(x):
    """approximation in the third order to the sin(x) when x is a
    numpy arrray of positive numbers"""
    twopi = 2 * pi    #2pi
    piH   = 0.5 * pi  #pi/2
    minuspiH = -piH   #-pi/2
    piX   = 1.5 * pi  #3pi/2

    z = x % twopi  # 0 < z < 2pi
    #shift Z by 2pi if piX < z, since sin(z-2pi)=sin(z) in this range
    z[ np.where(z > piX ) ] -= twopi  #z-->z-2pi, piX < z
    #change sign of z and shift by pi if piH < z, since
    # sin(z) = sin(pi-z) in this range
    z[ np.where(z > piH ) ] *= 1      #z-->-z, piH < z
    z[ np.where(z < minuspiH ) ] += pi  #z-->pi-z, piH < z

    #we have changed the 0<z<2pi range to -pi/2<z<pi/2 range while
    #preserving the value of sin(z), and now we use a Taylor expansion
    #to the third order in sin(z) for the -pi/2<z<pi/2 range. The
    #maximum error is 0.004 for abs(z)=pi/2
    z2 = z * z
    return z * ( 1 - 0.166666666 * z2 * ( 1 - 0.02 * z2 ) )


def Qlist(Qmin,Qmax,deltaQ):
    """generate list of Q-values
    Arguments:
      Qmin: minimum (included) Q-value
      Qmax: maximum (maybe not included) allowed Q-value
      deltaQ: separation between consecutive Q-values
    Returns:
      QQ: numpy array with the Q-values
    """
    QQ=[]; Q = Qmin
    while Q<=Qmax:
        QQ.append( Q ); Q += deltaQ
    return np.array(QQ)


def getFrames(universe, group, skip=1,format='fac'):
    """Extricate coordinates for a selected group of atoms
    Arguments:
      universe: universe object from MDAnalysis module
      group: group of atoms for selection
      skip: jump between frames
      format: shape of frames array, from permutting 'f', 'a', 'c'
              which are abbreviations for frames, atoms, coordinates.
              Example: frames[#frames][#atoms][3] for 'fac'
    Returns:
      frames: numpy array containing coordinate frames for the atoms selection
    """
    frames = universe.trajectory.timeseries(group, skip=skip, format=format)
    frames = np.array(frames,dtype='float32')
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

    b_ij = np.outer(b,b)
    s_Q = np.zeros(nQ)
    delta = max(1,int(nfr/nsampling ))#separat. between t0's
    t0 = 1  #initialize the time point
    isp = 0
    while t0 < nfr: #cycle over all t0's
        fi = frames[t0]
        dij=distance_array(fi,fi,None) #distances between atoms
        for iQ in range(nQ):
            arg = QQ[iQ]*dij
            #avoid sin(arg)/arg trouble when arg goes to zero
            arg[np.where(arg<co)]=co
            s_Q[iQ] += ( b_ij * np.sin(arg)/arg ).sum()/(nat*nat)
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

    b_ij = np.outer(b,b)
    bbcoh=(b*b).sum()*(np.sin(co)/co) #total self coherent term
    bbinc = b.sum()*b.sum()        #total incoherent term, squared
    s_Q = np.zeros(nQ)
    delta = max(1,int(nfr/nsampling ))#separat. between t0's
    t0 = 1  #initialize the time point
    isp = 0
    while t0 < nfr: #cycle over all t0's
        fi = frames[t0]
        dij=distance_array(fi,fi,None) #distances between atoms, exclude self
        for iQ in range(nQ):
            arg = QQ[iQ]*dij
            #avoid sin(arg)/arg trouble when arg goes to zero
            arg[np.where(arg<co)]=co
            #remove the self coherent term
            s_Q[iQ] += ( (b_ij * np.sin(arg)/arg).sum() - bbcoh)/bbinc
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

    w_ij = np.zeros(nQ*nat*nat).reshape(nQ,nat,nat)
    delta = max(1,int(nfr/nsampling ))#separat. between t0's
    t0 = 1  #initialize the time point
    while t0 < nfr: #cycle over all t0's
        fi = frames[t0]
        dij=distance_array(fi,fi,None) #distances between atoms
        for iQ in range(nQ):
            arg = QQ[iQ]*dij
            arg[np.where(arg<co)]=co#avoid sin(arg)/arg troubl
            w_ij[iQ] = np.sin(arg)/arg
        #print the structure factor
        w_ij *= np.outer(b,b)
        z_ij = w_ij.sum(axis=2).sum(axis=1)/(nat*nat)
        buf += ' %5d'%t0
        for iQ in range(nQ): buf += ' %10.7f'%z_ij[iQ]
        buf += '\n'
        t0 += delta #jump to next t0
    return {'buf':buf,}


def IC(frames,b,QQ,nsampling=200,nt=100,dt=1,
       comment='# intermediate coherent scattering factor\n'):
    """Intermediate coherent scattering function
    d_{i,j} = |r_i(t0)-r_j(t0+t)|
    IC = sum_{i,j} b_i * b_j * sin(Q*d_{i,j})/(Q*d_{i,j}) / len(b)^2

    frames: numpy.array with shape=(#frames,3,#atoms)
    b: numpy.array of coherent scattering lenghts
    QQ: numpy.array list of Q values
    nsampling: number of t0's for a given t
    nt: number of t's (number time lapses)
    dt: spacing between consecutive t's, in number of frames
    """

    from MDAnalysis.analysis.distances import distance_array

    nQ = len(QQ)
    nat = len(b)
    nfr = len(frames)
    bb =  np.outer(b,b)

    it = 0 #current number of evaluated t's
    t = 0

    buf= comment
    buf+='# Q-values=%s' % ' '.join( [str(Q) for Q in QQ] )+'\n'
    buf+='# time(psec) Q-values(A^(-1))\n'

    sf = np.zeros((1+nt)*nQ).reshape(1+nt,nQ)
    ts=[]
    while it < 1+nt: #cycle over all t's (+1.0 for t=0)
        s_ij = np.zeros(nQ*nat*nat).reshape(nQ,nat,nat) #store signal
        delta = max(1,int((nfr-t)/nsampling ))#separation between t0's
        isp = 0 #actual number of sampled t0 time points (isp<=nsampling)
        t0 = 1  #initialize the time point
        print 'it=%d t=%d'%(it,t)
        while t0+t < nfr: #cycle over all t0's
            #print 'it=%d t=%d isp=%d t0=%d'%(it,t,isp,t0)
            fi = frames[t0]
            fj = frames[t0+t]
            dij=distance_array(fi,fj,None) #distances between atoms
            for iQ in range(nQ):
                arg = QQ[iQ]*dij
                arg[np.where(arg<co)]=co#avoid sin(arg)/arg troubl
                s_ij[iQ] += np.sin(arg)/arg
            isp += 1 # update the number of sampled t0 time points
            t0 += delta #jump to next t0
        s_ij *= bb / isp
        s_ij = s_ij.sum(axis=2).sum(axis=1)/(nat*nat)
        sf[it] = s_ij
        ts.append(t)
        s = write_from_numpy(None,s_ij,ncols=nQ,format=' %10.7f')
        buf+=' %6.3f %s'%(t,s)
        it += 1
        t += dt
    return {'buf':buf,'sf':sf,'ts':np.array(ts)}

def II(frames,b,QQ,nsampling=200,nt=100,dt=1,begt=0,c1=np.empty(0),
       comment='# intermediate incoherent scattering factor\n'):
    """Intermediate incoherent scattering function
    Calculations:
      d_i = |r_i(t0)-r_(t0+t)|
      II = sum_i b_i * sin(Q*d_i)/(Q*d_i) / len(b)
      sin(Q*d_i)/(Q*d_i) denotes average over all orientations of Q-vector
    Arguments:
      frames: numpy.array with shape=(#frames,#atoms,3)
      b: numpy.array of incoherent scattering lenghts
      QQ: numpy.array list of Q values
      nsampling: number of t0's for a given t
      nt: number of t's
      dt: spacing between consecutive t's, in number of frames
      begt: smallest starting time lapse, in number of
            frames. Biggest time lapse = begt+nt*dt
      c1: correlation array of shape (#frames, #atoms). Will multiply each
          self-interference term by c1(t0)*c1(t0+t)
    Returns:
      Dictionary containing the following keys:
        buf: string containing II(Q,t)
        sf: II(Q,t) numpy array of shape (1+nt, nQ)
        ts: list of times t where II(Q,t) is calculated
      If c1 is passed as argument, additional keys:
        sfYY: II(Q,t) for correlation of c1(t0) and c1(t0+t)
        sfYN: II(Q,t) for correlation between c1(t0) and 1-c1(t0+t)
        sfNY: II(Q,t) for correlation between 1-c1(t0) and c1(t0+t)
        sfNN: II(Q,t) for correlation between 1-c1(t0) and 1-c1(t0+t)
    """
    isc1 = True if len(c1) else False
    nQ = len(QQ)
    nat = len(b)
    nfr = len(frames)
    bb = b*b # square of the incoherent scattering lengths
    it = 0   #current number of evaluated t's
    t = begt #time difference between t0 and t0+t
    buf= comment
    buf+='# Q-values=%s' % ' '.join( [str(Q) for Q in QQ] )+'\n'
    buf+='# time(psec) Q-values(A^(-1))\n'
    sf = np.zeros((1+nt)*nQ).reshape(1+nt,nQ) # stores II(t,Q)
    if isc1:
        sfYY = np.zeros((1+nt)*nQ).reshape(1+nt,nQ)
        sfYN = np.zeros((1+nt)*nQ).reshape(1+nt,nQ)
        sfNY = np.zeros((1+nt)*nQ).reshape(1+nt,nQ)
        sfNN = np.zeros((1+nt)*nQ).reshape(1+nt,nQ)
    ## Now cycle over all target time lapses (+1 for t=0)
    ts=[] # list of times
    while it < 1+nt:
        s_ii = np.zeros(nQ*nat).reshape(nQ,nat)
        if isc1:
            s_iiYY = np.zeros(nQ*nat).reshape(nQ,nat)
            s_iiYN = np.zeros(nQ*nat).reshape(nQ,nat)
            s_iiNY = np.zeros(nQ*nat).reshape(nQ,nat)
            s_iiNN = np.zeros(nQ*nat).reshape(nQ,nat)
        delta = max(1,int((nfr-t)/nsampling ))#separation between t0's
        isp = 0 #number of sampled t0 time points (isp<=nsampling)
        ## Now cycle over the sampling t0 points
        t0 = 1  #initialize the time point
        while t0+t < nfr:
            dii = frames[t0] - frames[t0+t]
            dii *= dii
            #dii distance traveled by each atom in the lapse t0 to to+t
            dii = np.sqrt(dii.sum(axis=1)) #sum over x, y,a nd z
            for iQ in range(nQ):
                arg = QQ[iQ]*dii
                #avoid sin(arg)/arg trouble when arg goes to zero
                arg[np.where(arg<co)]=co
                selfInterference = np.sin(arg)/arg # shape=(natoms,)
                s_ii[iQ] += selfInterference
                if isc1:
                     s_iiYY[iQ] += selfInterference*c1[t0]*c1[t0+t]
                     s_iiYN[iQ] += selfInterference*c1[t0]*(1-c1[t0+t])
                     s_iiNY[iQ] += selfInterference*(1-c1[t0])*c1[t0+t]
                     s_iiNN[iQ] += selfInterference*(1-c1[t0])*(1-c1[t0+t])
            isp += 1
            t0 += delta #jump to next t0
        s_ii *= bb/isp  # Normalize by number of sample t0 points (isp)
        s_ii = s_ii.sum(axis=1) # Incoherent means sum, i.e., no interferece
        if isc1:
            s_iiYY *= bb/isp
            s_iiYY = s_iiYY.sum(axis=1)
            s_iiYN *= bb/isp
            s_iiYN = s_iiYN.sum(axis=1)
            s_iiNY *= bb/isp
            s_iiNY = s_iiNY.sum(axis=1)
            s_iiNN *= bb/isp
            s_iiNN = s_iiNN.sum(axis=1)
        sf[it] = s_ii # stores II(t,Q) for the current t
        if isc1:
            sfYY[it] = s_iiYY
            sfYN[it] = s_iiYN
            sfNY[it] = s_iiNY
            sfNN[it] = s_iiNN
        ts.append(t)
        s = write_from_numpy(None,s_ii,ncols=nQ,format=' %10.7f')
        buf+=' %6.3f %s'%(t,s) # write additional I(t,Q) line to buffer
        it += 1 # next time t count
        t += dt # next time t
    results = {'buf':buf,'sf':sf,'ts':np.array(ts)}
    if isc1:
        results.update({'sfYY':sfYY, 'sfYN':sfYN, 'sfNY':sfNY, 'sfNN':sfNN})
    return results

def IIv2(frames,b,QQ,c1,nsampling=200,nt=100,dt=1,begt=0,
       comment='# intermediate incoherent scattering factor\n'):
    """Intermediate incoherent scattering function
    Calculations:
      d_i = |r_i(t0)-r_(t0+t)|
      II = sum_i b_i * sin(Q*d_i)/(Q*d_i) / len(b)
      sin(Q*d_i)/(Q*d_i) denotes average over all orientations of Q-vector
    Arguments:
      frames: numpy.array with shape=(#frames,#atoms,3)
      b: numpy.array of incoherent scattering lenghts
      QQ: numpy.array list of Q values
      nsampling: number of t0's for a given t
      nt: number of t's
      dt: spacing between consecutive t's, in number of frames
      begt: smallest starting time lapse, in number of
            frames. Biggest time lapse = begt+nt*dt
      c1: correlation array of shape (#frames, #atoms). Will multiply each
          self-interference term by c1(t0)*c1(to+dt)*..*c1(t0+t+dt)
    Returns:
      Dictionary containing the following keys:
        buf: string containing II(Q,t)
        sf: II(Q,t) numpy array of shape (1+nt, nQ)
        ts: list of times t where II(Q,t) is calculated
    """
    nQ = len(QQ)
    nat = len(b)
    nfr = len(frames)
    bb = b*b # square of the incoherent scattering lengths
    it = 0   #current number of evaluated t's
    t = begt #time difference between t0 and t0+t
    buf= comment
    buf+='# Q-values=%s' % ' '.join( [str(Q) for Q in QQ] )+'\n'
    buf+='# time(psec) Q-values(A^(-1))\n'
    sf = np.zeros((1+nt)*nQ).reshape(1+nt,nQ) # stores II(t,Q)
    ## Now cycle over all target time lapses (+1 for t=0)
    ts=[] # list of times
    while it < 1+nt:
        s_ii = np.zeros(nQ*nat).reshape(nQ,nat)
        delta = max(1,int((nfr-t)/nsampling ))#separation between t0's
        # isp is number of contributing sampled t0 time points for each atom
        isp = np.zeros(nat)
        ## Now cycle over the sampling t0 points
        t0 = 1  #initialize the time point
        while t0+t < nfr:
            dii = frames[t0] - frames[t0+t]
            dii *= dii
            #dii distance traveled by each atom in the lapse t0 to to+t
            dii = np.sqrt(dii.sum(axis=1)) #sum over x, y,a nd z
            # c is the product of elementary correlations
            c = np.prod(c1[t0:t0+t+1,],axis=0) # shape = (natoms,)
            for iQ in range(nQ):
                arg = QQ[iQ]*dii
                #avoid sin(arg)/arg trouble when arg goes to zero
                arg[np.where(arg<co)]=co
                selfInterference = np.sin(arg)/arg # shape=(natoms,)
                s_ii[iQ] += selfInterference * c
            isp += c
            t0 += delta #jump to next t0
        # if there are components i for which isp[i] is zero, then
        # s_ii[i]/isp[i] returns zero because of the previous np.seterr
        for i in range(nat):
            if isp[i]>0:
                s_ii[:,i] *= bb[i]/isp[i]  # Normalize by number of sample t0 points (isp)
            else:
                s_ii[:,i] = 0.0
        s_ii = s_ii.sum(axis=1) # Incoherent means sum, i.e., no interferece
        sf[it] = s_ii # stores II(t,Q) for the current t
        ts.append(t)
        s = write_from_numpy(None,s_ii,ncols=nQ,format=' %10.7f')
        buf+=' %6.3f %s'%(t,s) # write additional I(t,Q) line to buffer
        it += 1 # next time t count
        t += dt # next time t
    results = {'buf':buf,'sf':sf,'ts':np.array(ts)}
    return results

def saveIISassenaFormat(II, QQ, filename):
    """Save the scattering factor I(Q,t) computed with function II in as HDF5
    file with same entries as output from Sassena
    Arguments:
      II: I(Q,t) computed with function scatterMDA.II. Shape is (nt, nQ)
      QQ: associated list of Q values
      filename: string containing the output HDF5 file
    Returns:
      h5py._hl.files.File object, handle to the HDF5 file
    Except:
      ValueError when II.shape[1] != len(QQ)
    """
    # Validate
    ntimes, nQ = II.shape
    if nQ != len(QQ):
        raise ValueError("II.shape[1] != len(QQ)")
    f = h5py.File(filename, "w")
    # Save Q-values. Each Q is saved as a vector (Q,0,0)
    qSet = f.create_dataset("qvectors", (len(QQ),3), dtype='f')
    q3 = np.array(QQ).reshape([len(QQ),1])  # column-vector
    q3 = np.insert(q3, 1, 0, axis=1) # column of zeros
    q3 = np.insert(q3, 1, 0, axis=1)
    qSet[...] = q3
    # Save I(t=0,Q). Each has a real part and imaginary part (set to zero)
    fq0Set = f.create_dataset("fq0", (len(QQ),2), dtype='f')
    fq0 = II[0].reshape([len(QQ),1])  # column-vector
    fq0 = np.insert(fq0, 1, 0, axis=1) # column of zeros
    fq0Set[...] = fq0
    # Save I(t,Q). Each has a real part and imaginary part (set to zero)
    ntimes, nQ = II.shape
    fqtSet = f.create_dataset("fqt", (nQ, ntimes, 2), dtype='f')
    fqt = np.transpose(II) # now shape (nQ, ntimes)
    fqt = fqt.reshape(nQ, ntimes, 1)
    imag = np.zeros([nQ, ntimes, 1]) # imaginary component
    fqt = np.concatenate([fqt,imag], axis=2)
    fqtSet[...] = fqt
    # Temporary quick fix
    fqSet = f.create_dataset("fq", (len(QQ),2), dtype='f')
    fqSet[...] = fq0
    fq2Set = f.create_dataset("fq2", (len(QQ),2), dtype='f')
    fq2Set[...] = fq0
    # Always close the file to make sure we finish writing to file
    f.attrs["sassena_version"] = "1.4.1"
    f.close()
