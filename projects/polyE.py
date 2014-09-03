#!/usr/bin/python

import argparse, sys, os, numpy
from math import exp, erfc, sqrt
from pdb import set_trace as trace

projd = '/projects/research/polyE'

def single_chain(*kargs,**kwargs):
  jobList = ['',
             ]

def misc(*kargs,**kwargs):
  jobList = ['',
             ]
  job = kwargs['job']  #what job to do
  if job == '':
    """description"""
  else:
    print 'Job not recognized. Nothing to do'

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

