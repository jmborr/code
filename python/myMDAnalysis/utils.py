import numpy,pdb

def getFrames(U,group,skip=1,format='fac'):
    """return frames
       U: universe
       group: group of atoms (charmm atom selection syntax)
       skip: jump between frames
       format: shape of frames array, from permutting 'f', 'a', 'c'
               which are abbreviations for frames, atoms, coordinates.
               Example: frames[#frames][#atoms][3] for 'fac'
    """
    X = U.selectAtoms(group)
    frames = U.trajectory.timeseries(X,skip=skip,format=format)
    frames = numpy.array(frames,dtype='float32')
    return frames

def shortest_distance(U,group1,group2,skip=1):
    """For every atom in group1, find its shortest distance to any
       atom in group2 for each of the read frames

       group1: first group of atoms (charmm atom selection syntax)

       returns array of dimensions (#frames, #atoms-in-group1)
    """
    from MDAnalysis.analysis.distances import distance_array

    fr1 = getFrames(U,group1,skip=skip,format='fac')
    fr2 = getFrames(U,group2,skip=skip,format='fac')
    sd = []

    for ifr in range( len(fr1) ):
        dij=distance_array(fr1[ifr],fr2[ifr],box=U.coord.dimensions)
        sd.append( numpy.min(dij,axis=1) )
    return numpy.array(sd,dtype='float32')

def residence_times(U,group1,group2,ico,oco,skip=1,dt=1.0,minres=0):

    """For every atom in group1, find the atoms in group2 that are
    within distance ico (internal cut-off) and oco (outer cut-off)

    For a particular atom pair i,j (i from group1 and j from group2),
    we define three tipes of events (A, B, C)according to the
    evolution of their mutual distance at time t (d_{i,j}) and at time
    t-1 (pd_{i,j}):

     d_{i,j}
        ^
        |
        |\                                               _/
        | --\                                           /#C
     oco|----\-----------------------------------------/------
        |     \__                                  /\_/
        |        \                  _/\          _/
        |         \ #A             /   \_ #D    /  
     ico|----------\--------------/------\-----/--------------
        |           \    __     _/ #B     \   / #B
        |            \__/  \  _/           \_/
        |                   \/
    ----|------------------------------------------------------> time
        |

    event A: signal the particle entered the 'contact volume'
    if pd_{i,j}>ico and d_{i,j}<ico and tau0_{i,j}<0: tau0_{i,j}=t

    event B: particle left 'contact volume' and entered  'no-mans land'
    if pd_{i,j}<ico and d_{i,j}>ico: tau_{i,j}=t

    event c: particle left 'no-man land'. Store the residence time
    if pd_{i,j}<oco and d_{i,j}>oco:
       save(tau_{i,j}-tau0_{i,j}); tau0_{i,j}=-1

    The no-man land serves to avoid events when the distance leaves
    the contact volume just to come back after a little while (walk
    from #B to #D)

    skip: jump between consecutively read frames. 1 means read all frames
    dt: time-scale, time between consecutive frames
    minres: log only residence time bigger than minres. Units of minres are 1+skip
    """
    
    from MDAnalysis.analysis.distances import distance_array

    fr1 = getFrames(U,group1,skip=skip,format='fac')
    (nfr,N1,x)=fr1.shape
    X1 = U.selectAtoms(group1)
    fr2 = getFrames(U,group2,skip=skip,format='fac')
    (nfr,N2,x)=fr2.shape
    X2 = U.selectAtoms(group2)
    
    #initialize tau0, tau, and pd
    tau0 = -numpy.ones(N1*N2).reshape(N1,N2)
    tau = numpy.zeros(N1*N2,dtype='float32').reshape(N1,N2)
    pd = distance_array(fr1[0],fr2[0],box=U.coord.dimensions)

    #trajectory started. Treat as if suddenly distances entered the
    #'contact volume'
    for i in range(N1):
        for j in range(N2):
            if pd[i][j]<ico: tau0[i][j]=0
    
    #iterate over all frames
    rt=[] #list of residence times
    for ifr in range(1,nfr):
        if not ifr%100: print ifr
        d = distance_array(fr1[ifr],fr2[ifr],box=U.coord.dimensions)
        for i in range(N1):
            for j in range(N2):
                #event A
                if pd[i][j]>ico and d[i][j]<ico and tau0[i][j]<0: 
                    tau0[i][j]=ifr
                #event B
                elif pd[i][j]<ico and d[i][j]>ico:
                    tau[i][j]=ifr
                #event C
                elif pd[i][j]<oco and d[i][j]>oco and tau0[i][j]>0:
                    H = X1[i]; O = X2[j]
                    delta=tau[i][j]-tau0[i][j]
                    info=[delta,
                          H.resname,H.name,1+H.number,H.charge,
                          O.resname,O.name,1+O.number,O.charge]
                    if delta>minres: rt.append(info)
                    tau0[i][j] = -1
        pd = d #update previous distance

    #trajectory reached end. Treat as if suddenly distance left the
    #no-man land
    for i in range(N1):
        for j in range(N2):
            if tau0[i][j]>0:
                H = X1[i]; O = X2[j]
                info=[nfr-tau0[i][j],
                      H.resname,H.name,1+H.number,H.charge,
                      O.resname,O.name,1+O.number,O.charge]
                rt.append(info)

    #apply dt time scale to the residency times
    for irt in range( len(rt) ):  rt[irt][0] *= dt
        
    return rt
                    
