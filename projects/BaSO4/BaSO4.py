'''
Created on Jun 3, 2013

@author: jbq
'''
#!/usr/bin/python

import argparse, sys, os, numpy, random, re
from math import exp, erfc, sqrt
from pdb import set_trace as tr
from tempfile import mkstemp

projd='/projects/research/BaSO4'
coded='/home/jbq/code/projects/BaSO4'

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
  from mantid.simpleapi import (LoadNexus, SaveNexus, Scale, ScaleX)
  job=kwargs['job']

  if job=='create experimental Nexus files':
    rootd='/projects/research/POSS/me8t8/expdata'
    for T in (100,150,175,200,225,250,300,350):
      LoadDaveGrp(Filename='%s/exp%dK.dat'%(rootd,T),XAxisUnits='DeltaE',YAxisUnits='Empty',IsMicroEV=1,OutputWorkspace='exp%dK'%T)
      Rebin(InputWorkspace='exp%dK' %T,Params=[-0.14,0.0004,0.14],OutputWorkspace='exp%dK_2' %T)
      SaveNexus(InputWorkspace='exp%dK_2' %T,Filename='/projects/research/POSS/me8t8/expdata/exp%dK.nxs'%T)
    T=13
    LoadDaveGrp(Filename='%s/exp%dK.dat'%(rootd,T),XAxisUnits='DeltaE',YAxisUnits='Empty',IsMicroEV=1,OutputWorkspace='elastic')
    Rebin(InputWorkspace='elastic',Params=[-0.14,0.0004,0.14],OutputWorkspace='elastic')
    SaveNexus(InputWorkspace='elastic',Filename='/projects/research/POSS/me8t8/expdata/elastic.nxs')
    Scale(InputWorkspace='elastic',Factor=-1,Operation='Multiply',OutputWorkspace='resolution')
    SaveNexus(InputWorkspace='resolution',Filename='/projects/research/POSS/me8t8/expdata/resolution.nxs')

  else:
    print 'job not recognized, nothing to calculate'

def visualize(*kargs,**kwargs):
  """task related to preparing pictures or animations"""
  job=kwargs['job']

  if job=='short wrapped dcd trajectory':
    """save the trajectory from 1ns to 1.1ns with wrapped coordinates and with the BaSO4 in the center of the simulation box"""
    cmds='''mol new %s waitfor all
    mol addfile %s type lammpstrj first 10000 last 10999 waitfor all
    pbc set {34.4475 28.1422 57.2552} -all
    pbc box_draw
    set sel [atomselect top all]
    set nf 1000
    for {set i 0} {$i < $nf} {incr i} { 
      $sel frame $i 
      $sel moveby {0. 0. 13.5} 
    }
    pbc wrap -all
    pbc join residue -all
    animate write dcd %s
    quit'''%(kwargs['pdb'],kwargs['lammpstrj'],kwargs['dcd'])
    vmdexec(cmds)

  else:
    print 'job not recognized, nothing to calculate'


def preparing(*kargs,**kwargs):
  """tasks related to preparation of simulation system"""
  job=kwargs['job']

  if job=='translate to DCD':
    """translate from .lammpstrj to .dcd, neglecting the first nanosecond (treated as equilibration)"""
    cmds='''mol new %s waitfor all
    mol addfile %s type lammpstrj first 10001 waitfor all
    animate write dcd %s
    quit'''%(kwargs['pdb'],kwargs['lammpstrj'],kwargs['dcd'])
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

