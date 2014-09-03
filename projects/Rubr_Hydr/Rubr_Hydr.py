'''
Created on Jun 3, 2013

@author: jbq
'''
#!/usr/bin/python

import argparse, sys, os, numpy, random, re
from math import exp, erfc, sqrt
from pdb import set_trace as tr
from tempfile import mkstemp

projd='/projects/research/Rubr_Hydr'
coded='/home/jbq/code/projects/Rubr_Hydr'

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


def mantidjobs(*kargs, **kwargs):
  """Jobs carried out in Mantid"""

  job=kwargs['job']
  if job=='create experimental Nexus files':
    pass
  else:
    print 'job not recognized, nothing to calculate'


def sassena(*kargs,**kwargs):
  from utilities.readingWritingFiles import write_from_numpy
  job=kwargs['job']
  if job=='convert ptraj DCD to vmd DCD':
    """sassena cannot read a DCD trajectory output by ptraj. We must load it into VMD and then
    save it as DCD. Now sassena can read the new DCD trajectory"""
    cmds="""mol new %s waitfor all
mol addfile %s waitfor all
animate write dcd %s
quit"""%(kwargs['prmtop'],kwargs['indcd'],kwargs['outdcd'])
    vmdexec(cmds)
  else:
    print 'job not recognized, nothing to calculate'


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

