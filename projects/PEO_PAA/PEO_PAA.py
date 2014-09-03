#!/usr/bin/python

import argparse
import sys
import os
import numpy
from pdb import set_trace as trace

sys.path.append('/home/jmborr/devel/mantidproject/mantid/Code/debug/bin')
PROD='/projects/research/PEO_PAA'

def addSassenaVersion(*kargs, **kwargs):
  """insert Sassena version into an HDF5 file"""
  from mysassena.version import addVersionStamp
  addVersionStamp(kwargs['filename'],'1.4.1')
  return None

def reOrderByQmodulus(*kargs, **kwargs):
  """ Ordering in the H5 by qvector modulus is not enforced by default """
  from mysassena.orderByQmodulus import orderByQmodulus
  orderByQmodulus(kwargs['filename'])
  return None

def saveAsAscii(*kargs, **kwargs):
  """Load a structure factor and save the fqt.Re workspace as ASCII file"""
  from mantid.simpleapi import mtd,LoadSassena,Scale,SaveAscii
  LoadSassena(Filename=kwargs['filename'],OutputWorkspace='PEO-AA',TimeUnit=1)
  ws=mtd['PEO-AA_fq0']
  I=ws.readY(0)[0]
  Scale(InputWorkspace='PEO-AA_fqt.Re',OutputWorkspace='PEO-AA_fqt.Re',Operation='Multiply',Factor=1.0/I)
  SaveAscii(InputWorkspace='PEO-AA_fqt.Re',Filename=kwargs['outfile'],CommentIndicator='#',Separator='Tab')

def createPNGs(*kargs, **kwargs):
  """Load ASCII file containing I(Q,t) and create PNG images for each Q with xmgrace"""
  template='''READ BLOCK "_INFILE_"
BLOCK xy "1:_INDEX_"
s0 line color 1
title "_SYSTEM_ - permittivity=_PERMITTIVITY_ - Q=_Q-VALUE_\cE\C\S-1\N"
xaxis label "time (reduced-units)"
yaxis label "_TYPE_ I(Q,t) (arb-units)"
PRINT TO "_OUTFILE_"
HARDCOPY DEVICE "PNG"
PAGE SIZE 2560, 2048
DEVICE "PNG" FONT ANTIALIASING on
DEVICE "PNG" OP "transparent:off"
DEVICE "PNG" OP "compression:9"
PRINT
'''
  substitutions={'_INFILE_':kwargs['infile'],
                 '_SYSTEM_':kwargs['system'],
                 '_PERMITTIVITY_':kwargs['permittivity'],
                 '_TYPE_':kwargs['type'],
                 }
  tpl=template
  for (key,val) in substitutions.items(): tpl=tpl.replace(key,val)
  for i in range(1,61):
    print i
    tpl2=tpl
    tpl2=tpl2.replace('_INDEX_',str(2*i))
    Q=i*0.05
    tpl2=tpl2.replace('_Q-VALUE_',str(Q))
    tpl2=tpl2.replace('_OUTFILE_',kwargs['outfileroot']+'_Q%s.png'%Q)
    open('junk.createPNGs','w').write(tpl2)
    os.system('xmgrace -batch "junk.createPNGs" -nosafe -hardcopy')

def createHTMLreport(*kargs, **kwargs):
  """Make an HTML page showing the different I(Q,t)"""
  for system in ('cion','charge'):
    for typefactor in ('coherent', 'incoherent'):
      for permittivity in range(2,12,2):
        pstr='%02d'%permittivity
        hs = '<!DOCTYPE html>\n<html>\n<body>\n\n'
        hs += '<h1>'+typefactor+' I(Q,t) for '+system+pstr+'</h1>\n'
        hs += '<p>These plots show '+typefactor+' intermediate structure factors I(Q,t) for 0.05<Q<3 and permitivities 2<=p<=10</p>\n'
        for iQ in range(1,61):
          Qstr=str(iQ*0.05)
          imagefile=system+pstr+'_'+typefactor+'IQt_Q'+Qstr+'.png'
          hs+='<img src="figures/'+imagefile+'" alt="I(Q,t)" width="550" height="400">\n'
        hs += '</body>\n</html>'
        open( system+pstr+'_'+typefactor+'_'+kwargs['outfile'],'w').write(hs)


def prepSys( *kargs, **kwargs ):
    joblist = ['something',
               ]
    JOB = kwargs['job'] #what job to do

    if JOB=='something':
        """do something"""
        pass

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
    if args.kwargs:
        optargs=dict( map( lambda x: x.split('='), args.kwargs.split(',') ) )
    exitCode=0
    try:
        locals()[service](*reqargs,**optargs)
    except:
        sys.stderr.write('Error: service '+service+' not found or error in service\n')
        exitCode=1
    sys.exit(exitCode)

