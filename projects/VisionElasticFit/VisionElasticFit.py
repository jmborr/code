#!/usr/bin/python

import argparse, sys, os, numpy
from math import exp, erfc, sqrt
from pdb import set_trace as trace
import ContactMapAnalysis.ContactMapAnalysisAPI as CMAPI
import ContactMapAnalysis.ContactMapAnalysis as CMA

projd = '/projects/development/VisionElasticFit'

def back2backExp(x,params):
  I=params['I'];  a=params['A'];  b=params['B'];  x0=params['X0'];  s=params['S']
  s2 = s*s
  diff = x-x0
  y1 = 0.0
  try:
    y1 = exp( a/2*(a*s2+2*diff) ) * erfc( (a*s2+diff) / sqrt(2*s2) )
  except:
    y1 = 0.0
  y2 = 0.0
  try:
    y2 = exp( b/2*(b*s2-2*diff) ) * erfc( (b*s2-diff) / sqrt(2*s2) )
  except: 
    y2 = 0.0
  return I * ( y1 + y2 )

def compressEvents(*kargs,**kwargs):
  """Given an event workspace, compress the events"""
  sys.path.append('/home/jmborr/devel/mantidproject/mantid/Code/debug/bin')
  from mantid.simpleapi import Load,CompressEvents,SaveNexus,DeleteWorkspace
  filename=os.path.basename(kwargs['Filename'])
  outfilename=kwargs['Filename']+'compress'+str(kwargs['tolerance'])+'us.nxs'
  if 'FileOut' in kwargs.keys(): outfilename=kwargs['FileOut']
  try:
    for ibank in kwargs['banks'].split():
      """Load bank by bank, preferable when the file is bigger than 2GB"""
      Load(Filename=kwargs['Filename'],BankName='bank'+ibank,OutputWorkspace=filename)
      CompressEvents(InputWorkspace=filename,OutputWorkspace=filename,tolerance=float(kwargs['tolerance']))
      append=False
      if(os.path.exists(outfilename)): append=True
      SaveNexus(InputWorkspace=filename,Filename=outfilename,Append=append)
      DeleteWorkspace(filename)
  except:
    Load(Filename=kwargs['Filename'],OutputWorkspace=filename)
    CompressEvents(InputWorkspace=filename,OutputWorkspace=filename,tolerance=float(kwargs['tolerance']))
    SaveNexus(InputWorkspace=filename,Filename=outfilename)
  return outfilename

def conjoinWorkspaces(*kargs,**kwargs):
  """ Conjoin the group workspace after applying compressEvents """
  sys.path.append('/home/jmborr/devel/mantidproject/mantid/Code/debug/bin')
  from mantid.simpleapi import Load,ConjoinWorkspaces,RenameWorkspace,SaveNexus
  workspace=os.path.basename(kwargs['Filename']).split('.nxscompress')[0] # e.g VIS_1085_event.nxscompress0.1us.nxs
  Load(Filename=kwargs['Filename'],OutputWorkspace=workspace)
  ConjoinWorkspaces(InputWorkspace1=workspace+'_1',InputWorkspace2=workspace+'_2')
  ConjoinWorkspaces(InputWorkspace1=workspace+'_1',InputWorkspace2=workspace)
  RenameWorkspace(InputWorkspace=workspace+'_1',OutputWorkspace=workspace)
  SaveNexus(InputWorkspace=workspace,Filename=kwargs['outFilename'])
  return None

def misc(*kargs,**kwargs):
  jobList = ['fitted plots',
             'find peak heights',
             'find position of the maximum for the fit',
             'create HTML report',
             'emission time fit muParser',
             ]
  JOB = kwargs['job']  #what job to do
  if JOB == 'fitted plots':
    """generate a simple PNG for every pixel_???.dat file """
    template="""#From http://ringo.ams.sunysb.edu/index.php/Xmgrace
READ NXY "pixel__ID_.dat"
title "detector _ID_"
xaxis label "TOF (micro-seconds)"
yaxis label "Intensity"
AUTOSCALE
PRINT TO "pixel__ID_.png"
HARDCOPY DEVICE	"PNG"
PAGE SIZE 2560,	2048
DEVICE "PNG" FONT ANTIALIASING on
DEVICE "PNG" OP "compression:9"
PRINT
"""
    for idx in range(535,621):
      script = template.replace( '_ID_', str(idx) )
      open( '%s/junk.txt' % kwargs['inDir'], 'w' ).write( script )
      command = 'cd %s; xmgrace -batch junk.txt -nosafe -hardcopy' % kwargs['inDir']
      os.system( command )
  elif JOB == 'find peak heights':
    """Read file parameters.dat and find out the height of each back2backExponential"""
    inF = open( kwargs['parFile'] )
    ouF = open( kwargs['outFile'], 'w' )
    ouF.write( '#Elastic Adjacent\n' )
    for line in inF.readlines():
      if line[0] == '#': continue
      Z = [ float( x ) for x in line.split() ]
      elastic = { 'I':Z[1], 'X0':Z[2], 'A':Z[3], 'B':Z[4], 'S':Z[5] }
      yE = back2backExp( elastic['X0'],elastic )
      adjacent = { 'I':Z[6], 'X0':Z[7], 'A':Z[8], 'B':Z[9], 'S':Z[10] }
      yA = back2backExp( adjacent['X0'], adjacent )
      ouF.write( '%d %f %f\n' % (int(Z[0]), yE, yA) )
      ouF.close()
    inF.close()
  elif JOB == 'find position of the maximum for the fit':
    """Read file pixel_???.dat and find the position of the maximum for the fit"""
    from utilities.readingWritingFiles import read_to_cols
    outF = open( kwargs['outFile'], 'w' )
    outF.write('# detector maximum\n')
    for ID in range(535,621):
      fname = '/projects/development/VisionElasticFit/test/b2blb/pixel_%d.dat' % ID
      cols = read_to_cols( fname, outFmt='numpy.array',  xtype='float' )
      cols = cols.transpose()
      i = cols[2].argmax(axis=0)  #position of the maximum for the fit
      #trace()
      outF.write('%d %f\n'%( ID, cols[0][i] ))
    outF.close()
  elif JOB == 'create HTML report':
    hs = '<!DOCTYPE html>\n<html>\n<body>\n\n'
    hs += '<h1>Peak Center Positions</h1>\n'
    hs += '<p>The plot shows the evolution of the elastic and the adjacent line as a function of detector ID</p>\n'
    hs += ' <img src="./b2blb/peakCenter.png" alt="Peak center evolution with detector ID" width="1100" height="800">\n'
    hs += '<p>The fitting procedure has a hard time finding the position of the elastic line for the detector range 590-610, as can be seen from change in the slope of the elastic line. However, the slope for the position of the maximum for the fit of the whole spectrum is smooth</p>\n'
    hs += '<h1>Peak Intensity</h1>\n'
    hs += '<p>The plot shows the evolution of the peak intensity as a function of detector ID</p>\n'
    hs += '<img src="./b2blb/peakHeights.png" alt="Peak intensity evolution with detector ID" width="1100" height="800">\n'	    
    hs += '<h1>Fit Plots</h1>\n'
    hs += '<p>Below are simple plots of the fit in the neighborhood of the elastic line. The model is two back-to-back-exponential peaks and one linear background</p>\n'
    for i in range(535, 621):
      hs += '<img src="./b2blb/pixel_%d.png" alt="spectrum %d" width="1100" height="800">\n' % ( i, i )
    hs += '</body>\n</html>'
    open( kwargs['outFile'], 'w' ).write( hs )
  elif JOB=='emission time fit muParser':
    """Read in the selected (E_i, t_0) points to create a multi-segmented line of log(E_i) versus log(t_0)
    Each line will have the functional form t_0=a/E^alpha"""
    from math import log,exp
    pairs=[]; buf=''
    for l in open(kwargs['infile'],'r').readlines():
      if '#' not in l: pairs.append([float(x) for x in l.split()])
    for i in range(1,len(pairs)):
      E1=pairs[i-1][0]; e1=log(E1); T1=pairs[i-1][1]; t1=log(T1)
      E2=pairs[i][0]; e2=log(E2); T2=pairs[i][1]; t2=log(T2)
      alpha=-(t2-t1)/(e2-e1)
      a=(e1*t2-e2*t1)/(e1-e2); a=exp(a)
      buf+='(incidentEnergy &lt; %f) ? %f*incidentEnergy^(-%f) : '%(E2,a,alpha)
    buf += '0.0'
    print buf
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
