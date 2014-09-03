#!/usr/bin/python
import numpy,copy

def rmsd( X, Y ):
    """X.shape=(N,3), Y.shape=(N,3)
    rmsd fitting of X onto Y:  X --> U * (X-CMx) + CMy
    """
    L=len(X)
    CMx = numpy.sum(X,axis=0) / float(L)
    CMy = numpy.sum(Y,axis=0) / float(L)

    Xt = copy.copy(X) - CMx
    Yt = copy.copy(Y) - CMy

    # Initial residual, see Kabsch.
    E0 = numpy.sum( numpy.sum(Xt*Xt,axis=0),axis=0) + numpy.sum( numpy.sum(Yt*Yt,axis=0),axis=0)

    # This beautiful step provides the answer.  V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(Yt), Xt) )

    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))

    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = numpy.sqrt(abs(RMSD / L))

    #U is simply V*Wt
    U = numpy.dot(V, Wt)

    return RMSD,U,CMx-CMy

def euler2rot(xyz):
    """return rotation matrix given Eugler angles in the following sequence:
    1- rotate by an angle x counterclockwise about Z-axis
    2- rotate by an angle y counterclockwise about new X-axis
    3- rotate by an angle z counterclockwise about new Z-axis   

    xyz is numpy 1D array with x, y, and z angles
    Angles in radians

    Taken from the book "Classical Mechanics" by Herbet Goldstein, 2nd
    Ed., page 147.

    """
    x=xyz[0]
    a = numpy.cos(x)
    b = numpy.sin(x)
    D = numpy.array([a,b,0,-b,a,0,0,0,1]).reshape(3,3)

    y=xyz[1]
    a = numpy.cos(y)
    b = numpy.sin(y)
    C = numpy.array([1,0,0,0,a,b,0,-b,a]).reshape(3,3)

    z=xyz[2]
    a = numpy.cos(z)
    b = numpy.sin(z)
    B = numpy.array([a,b,0,-b,a,0,0,0,1]).reshape(3,3)

    A = numpy.dot(B, numpy.dot(C,D) )
    return A
