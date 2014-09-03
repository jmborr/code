'''
Created on Feb 13, 2013

@author: jmborr
'''
from pdb import set_trace as trace

def mockSassenaOutput(filename,stamp='1.4.1',kind='Gaussian+i*tGaussian'):
  """Create mock Sassena output.
  
  kind==Gaussian+i*tGaussian creates I(Q,t) with:
     - only one Q-value (Q=1),
     - Real[I(Q,t)] = exp(-t^2/2)
     - Imag[I(Q,t)] = t*exp(-t^2/2)
  The Fourier transform of exp(-t^2/2) is sqrt(2*pi)*exp(-E^2/2)
  The Fourier transform of t*exp(-t^2/2) is -i*sqrt(2*pi)*E*exp(-E^2/2)
  Thus S(Q,E) = sqrt(2*pi)*exp(-E^2/2)*(1+E)
  
  """
  import h5py
  import numpy
  if kind=='Gaussian+i*tGaussian':
    Q=numpy.array((1.,0.,0.))
    Q=Q.reshape(1,3)
    x=numpy.arange(0,10,0.01)
    L=len(x)
    x=x.reshape(1,L,1)
    real=numpy.exp(-x*x/2)
    imag=x*numpy.exp(-x*x/2)
    fqt=numpy.concatenate((real,imag),axis=2)
    f=h5py.File(filename, 'w')
    f.attrs['sassena_version']=stamp
    f.create_dataset('qvectors',data=Q)
    f.create_dataset('fqt',data=fqt)
    f.close()