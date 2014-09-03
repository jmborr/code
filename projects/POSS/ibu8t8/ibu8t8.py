'''
Created on Jun 3, 2013

@author: jbq
'''
#!/usr/bin/python

import argparse, sys, os, numpy, random, re
from math import exp, erfc, sqrt
from pdb import set_trace as tr
from tempfile import mkstemp

projd='/projects/research/POSS/ibu8t8'
coded='/home/jbq/code/projects/POSS/ibu8t8'

def vmdexec(cmds,redirect=False):
  """execute a set of commands by calling vmd"""
  handle,filename=mkstemp(dir='/tmp')
  open(filename,'w').write(cmds)
  cmd='vmd -dispdev text -e %s'%filename
  if redirect:
    return os.popen(cmd)
  else:
    os.system(cmd)
  return None

def single_molecule(*kargs,**kwargs):
  from utilities.readingWritingFiles import write_from_numpy
  job = kwargs['job']
  if job=='remove residual total charge':
    """Remove residual total charge"""
    pt=open(kwargs['prmtop'])
    buf=''
    line=pt.readline()
    while '%FLAG CHARGE' not in line:
      buf+=line
      line=pt.readline()
    buf+=line+pt.readline()
    line=pt.readline()
    charges=[]
    while '%FLAG MASS' not in line:
      charges+=[float(q) for q in line.split()]
      line=pt.readline()
    charges=numpy.array(charges)
    charges=charges-numpy.sum(charges)/len(charges)
    buf+=write_from_numpy(None,charges,format='%16.8e',ncols=5)
    buf+=line
    buf+=''.join(pt.readlines())
  open(kwargs['outprmtop'],'w').write(buf)

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='service provider for dhfr_solv project')
  parser.add_argument('service',help='requested service, the name of a function defined in this module')
  parser.add_argument('--kargs',help='required arguments of service. Ex: "arg1,arg2,arg3"')
  parser.add_argument('--kwargs',help='optional arguments of service. Ex: "arg1=val1,arg2=val2"')
  args = parser.parse_args()
  service=args.service
  reqargs=[]
  if args.kargs: reqargs=args.kargs.split(',')
  optargs={}
  if args.kwargs: optargs=dict( map( lambda x: x.split('='), args.kwargs.split(',') ) )
  exitCode=0
  try:
    locals()[service](*reqargs,**optargs)
  except:
    sys.stderr.write('Error: service '+service+' not found or error in service\n')
    exitCode=1
  sys.exit(exitCode)

