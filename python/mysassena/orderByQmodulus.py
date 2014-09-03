#from pdb import set_trace as trace # Uncomment only for debugging purposes
import os
from tempfile import mkstemp

def orderByQmodulus(filename,outfile=None):
  """ Sassena does not enforce any ordering of the structure factors.
  Here we order by increasing value of modulus of Q-vectors. """
  import h5py
  import numpy
  f=h5py.File(filename,'r')
  overwrite=False
  if not outfile:
    handle,outfile=mkstemp(dir='/tmp') # temporaty output file
    overwrite=True
  g=h5py.File(outfile,'w')
  ds_q = numpy.array(f["qvectors"]) # shape==(nvectors,3)
  moduli=numpy.square(ds_q).sum(axis=1) # moduli-squared of the Q-vectors
  rank=numpy.argsort(moduli) # rank from smallest to greatest
  for dset in ('qvectors', 'fqt', 'fq', 'fq0', 'fq2'):
    if dset in f.keys(): 
      x=numpy.array(f[dset])
      if not outfile:
        del f[dset]
        f[dset]=x[rank]
      else:
        g[dset]=x[rank]
  for key,val in f.attrs.items(): g.attrs[key]=val
  g.close()
  f.close()
  if overwrite:
    os.system('/bin/mv %s %s'%(outfile,filename))
  return None


if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='Reorder datasets by Q-modulus')
  parser.add_argument('--filename',help='sassena output file, serves as input')
  parser.add_argument('--outfile', help='Optional, write ordered file to outfile. Otherwise rewrite filenam')
  args = parser.parse_args()
  orderByQmodulus(args.filename,outfile=args.outfile)
