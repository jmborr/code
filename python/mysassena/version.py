'''
Created on Feb 8, 2013

@author: jmborr
'''
import h5py
#from pdb import set_trace as trace  # Uncomment only for debugging purposes

def hasVersion(filename):
  """Check filename as sassena version"""
  f = h5py.File(filename,'r')
  value=False
  if 'sassena_version' in f.attrs.keys(): value=True
  f.close()
  return value

def addVersionStamp(filename,stamp):
  """ Insert stamp as version attribute in and HDF5 file. """
  f = h5py.File(filename,'r+')
  f.attrs['sassena_version']=stamp
  f.close()

if __name__=='__main__':
  import argparse
  parser = argparse.ArgumentParser(description='add version to file')
  parser.add_argument('--filename',help='sassena output file, serves as input')
  parser.add_argument('--version',help='optional, sassena version string. Default is "1.4.1"')
  args = parser.parse_args()
  if not args.version: args.version='1.4.1'
  addVersionStamp(args.filename, args.version)

