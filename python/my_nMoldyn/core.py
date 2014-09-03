import pdb,time
from nMoldyn.misc import timePrepare, BincohList, logfileInit,logfileUpdate,cleanUp
from Numeric import *
import Numeric; N = Numeric
from Scientific.IO.NetCDF import NetCDFFile
from nMoldyn.calc import AutoCorrelationFunction

def IncoherentScatteringFunction(traj,qVect,ncDataFN,
                                 atoms=None,bincoh=1.0,timeInfo=None,
                                 nsteps=None):
    """notice bincoh=1"""
    
    file, filename = logfileInit('isf')

    if atoms is None: atoms = traj.universe
    tim    = timePrepare(traj,timeInfo)
    if nsteps is None:
        nsteps = len(tim)
    nsteps = min(nsteps, len(tim))

 
    comment = ''
    try:
        comment = comment + 'Trajectory: ' + traj.filename + '\n'
    except AttributeError:
        pass

    #print 'initializing output file'
    file2       = NetCDFFile(ncDataFN, 'w', comment)
    file2.title = 'Incoherent Scattering Function'
    file2.createDimension('Q', None)
    file2.createDimension('TIME', nsteps)
    SF      = file2.createVariable('sf', Float, ('Q','TIME'))
    Time    = file2.createVariable('time', Float, ('TIME',))
    Qlength = file2.createVariable('q', Float, ('Q',))
    Time[:]=tim[:nsteps]


    for j in range(len(qVect[1])):
        qVect[1][j] = transpose(array(map(lambda v: v.array, qVect[1][j])))

    normFact = len(qVect[0])
    Fincoh = zeros((nsteps,len(qVect[1])),Float)
    for j in range(len(qVect[0])):
        for at in atoms.atomList():
            #print 'qi=%d readParticleTrajectory  %s'%(j,at.name)
            sekw        = traj.readParticleTrajectory(at,first=timeInfo[0],\
                               last=timeInfo[1],skip=timeInfo[2]).array
            series      = exp(1j*dot(sekw, qVect[1][j]))
            #print 'AutoCorrelationFunction'
            res         = AutoCorrelationFunction(series)[:nsteps]
            Fincoh[:,j] = Fincoh[:,j] + add.reduce(bincoh*res,1)/\
                                        qVect[1][j].shape[1]
        SF[j,:]    = Fincoh[:,j]
        Qlength[j] = qVect[0][j]
        file2.flush()
        gj = logfileUpdate(file,j,normFact)

    cleanUp(file,filename)
    file2.close()
    return 
