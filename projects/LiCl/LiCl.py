#!/usr/bin/python

import argparse, sys, os, numpy
from math import exp, erfc, sqrt
from pdb import set_trace as trace

sys.path.append('/home/jmborr/devel/mantidproject/mantid/Code/debug/bin')

projd = '/projects/research/LiCl'

def NormalizeSpectrum(ws, start_x=None, end_x=None):
  """Normalize spectrum to unity in the given range"""
  from mantid.simpleapi import ConvertToHistogram,NormaliseToUnity
  ConvertToHistogram(InputWorkspace=ws,OutputWorkspace=ws)
  kwargs={}
  if start_x: kwargs['RangeLower']=start_x
  if end_x: kwargs['RangeUpper']=end_x
  NormaliseToUnity(InputWorkspace=ws,OutputWorkspace=ws,**kwargs)

def LoadEugeneData( strT, dataDir ):
  """load files into Mantid workspaces"""
  from mantid.simpleapi import LoadAscii
  for Qval in ( '3', '5', '7', '9' ):
    fname = 'LiCl_0_%s_%sK_ForDAVE.TXT' % ( Qval, strT )
    LoadAscii( Filename = dataDir + '/' + fname, Unit = 'DeltaE',Separator='Tab', OutputWorkspace = 'data_0_%s_%s' % ( Qval, strT ) )
    fname = 'LiCl_0_%s_006K_ForDAVE.TXT' % Qval
    LoadAscii( Filename = dataDir + '/' + fname, Unit = 'DeltaE',Separator='Tab', OutputWorkspace = 'data_0_%s_006' % Qval )

def prepareSystem(*kargs,**kwargs):
  #jobList = ['format PDB for VMD',]
  job = kwargs['job']  #what job to do

  if job == 'format PDB for VMD':
    """Slightly change the format of the PDB that was output by
    ptraj so that VMD can read it """
    newpdb = ''
    for line in  open(kwargs['inpdb']).readlines():
      if 'Li+' in line:
        newpdb += line.replace('Li+ Li+  ','LIT LIT X').strip() + '      O1  LI\n'
      elif 'Cl-' in line:
        newpdb += line.replace('Cl- Cl-  ','CLA CLA X').strip() + '      O2  CL\n'
      elif ' O ' in line:
        newpdb += line.replace('O   WAT  ','OH2 TIP3X').strip() + '      W1   O\n'
      elif ' H1 ' in line:
        newpdb += line.replace('H1  WAT  ','H1  TIP3X').strip() + '      W1   H\n'
      elif ' H2 ' in line:
        newpdb += line.replace('H2  WAT  ','H2  TIP3X').strip() + '      W1   H\n'
    newpdb += 'END\n'
    open(kwargs['outpdb'],'w').write(newpdb)

  elif job=='output experimental nexus':
    """ output experimental spectra to a Nexus file for a given temperature"""
    from mantid.simpleapi import ScaleX,CloneWorkspace, AppendSpectra, SaveNexus
    strT=kwargs['T']
    start_x=None
    if 'StartX' in kwargs.keys(): start_x=float(kwargs['StartX'])
    end_x=None
    if 'EndX' in kwargs.keys(): end_x=float(kwargs['EndX'])
    LoadEugeneData( strT, kwargs['dataDir'] )
    Qvals=['3', '5', '7', '9']
    Qval=Qvals.pop(0)
    iws='data_0_%s_%s' % ( Qval, strT )
    wsname='data_0_%s' %strT
    ScaleX(InputWorkspace=iws,OutputWorkspace=iws,factor=0.001) # input energies are in micro-eV
    CloneWorkspace(InputWorkspace=iws, outputWorkspace=wsname)
    while Qvals:
      Qval=Qvals.pop(0)
      iws='data_0_%s_%s' % ( Qval, strT )
      ScaleX(InputWorkspace=iws,OutputWorkspace=iws,factor=0.001) # input energies are in micro-eV
      AppendSpectra(InputWorkspace1=wsname, InputWorkspace2=iws,OutputWorkspace=wsname)
    NormalizeSpectrum(wsname, start_x=start_x, end_x=end_x)
    SaveNexus(InputWorkspace=wsname, Filename=kwargs['expNexus'])
  elif job=='output simulated nexus':
    """ convert I(Q,t) to S(Q,E) and save to a nexus file"""
    from simulation.src.scattering.sassenatasks import genSQE
    hdfname=kwargs['sassenafile']
    nxsname=kwargs['outfile']
    startX=float(kwargs['startX'])
    endX=float(kwargs['endX'])
    temp=float(kwargs['T'])
    indexes=[ int(x) for x in kwargs['indexes'].split() ]
    genSQE(hdfname,nxsname,wsname='simulated',indexes=indexes,
           LoadSassena={'TimeUnit':1.0,}, #simulation output every 1ps
           SassenaFFT={'FFTonlyRealPart':1,'DetailedBalance':1,'Temp':temp},
           SaveNexus={'Title':'Simulated S(Q,E)',},
           NormaliseToUnity={'RangeLower':startX,'RangeUpper':endX}
           )
  elif job=='create normalized convolution':
    """Load nexus S(Q,E) and convolve with resolution nexus file,
      then use Mantid algorithm NormaliseToUnity to normalize, then save as
      Nexus file
    """
    from mantid.simpleapi import (LoadNexus, Rebin, NormaliseToUnity, SaveNexus)
    from numpy import convolve
    wss=LoadNexus(Filename=kwargs['simulated'], OutputWorkspace='simulated')
    v=wss.readX(0);  start_x=v[0];  end_x=v[-1]
    wsr=LoadNexus(Filename=kwargs['resolution'], OutputWorkspace='resolution')
    v=wsr.readX(0);  width=v[1]-v[0] # retrieve the (constant) width from the first bin
    wss=Rebin(InputWorkspace='simulated', Params=(start_x,width,end_x), OutputWorkspace='simulated') #make same bin width
    for i in range(wsr.getNumberHistograms()):
      v=wsr.readY(i)
      w=wss.readY(i)
      x=convolve(v,w,mode='same')
      wss.setY(i,x)
    v=wsr.readX(0)
    Rebin(InputWorkspace='simulated', Params=(v[0],width,v[-1]), OutputWorkspace='simulated')
    NormaliseToUnity(InputWorkspace='simulated', OutputWorkspace='convolved')
    SaveNexus(InputWorkspace='convolved', Filename=kwargs['convolved'])
  elif job=='create resolution file':
    """ Load nexus file at 6K, the elastic line, and produce a mirror image,
    which is the resolution function"""
    from mantid.simpleapi import (LoadNexus, ScaleX, SaveNexus)
    LoadNexus(Filename=kwargs['elasticLine'],OutputWorkspace='elastic')
    ScaleX(InputWorkspace='elastic', OutputWorkspace='resolution',factor=-1)
    SaveNexus(InputWorkspace='resolution', Filename=kwargs['resolution'])
  elif job=='create normalized resolution file':
    """Load nexus file at 6K, the elastic line, and follow these steps:
    1- produce a mirror image for each workspace (each Q-value)
    2- average all Q-slices and substitute each Q-slice with the average
    3- use Mantid algorithm NormaliseToUnity (recall this normalisation
       does not take into account bin width
    """
    from mantid.simpleapi import (LoadNexus, ScaleX, SumSpectra, NormaliseToUnity, CloneWorkspace, AppendSpectra, SaveNexus)
    wse=LoadNexus(Filename=kwargs['elasticLine'],OutputWorkspace='elastic')
    ScaleX(InputWorkspace='elastic', OutputWorkspace='resolution',factor=-1)
    SumSpectra(InputWorkspace='resolution',OutputWorkspace='resolution')
    CloneWorkspace(InputWorkspace='resolution',OutputWorkspace='temp')
    for _ in range(wse.getNumberHistograms()-1): 
      AppendSpectra(InputWorkspace1='resolution',InputWorkspace2='temp',OutputWorkspace='resolution')
    NormaliseToUnity(InputWorkspace='resolution',OutputWorkspace='resolution')
    SaveNexus(InputWorkspace='resolution', Filename=kwargs['resolution'])
  else:
    print 'Job not recognized. Nothing to do'

def runSystem( *kargs, **kwargs ):

  #joblist = ['extract initial conformations from annealing','extract extended info from PDB header']
  JOB = kwargs['job'] #what job to do
  if JOB=='extract initial conformations from annealing':
    """extract coordinates from the annealing run"""
    from tempfile import mkstemp
    #template VMD script to output coordinates
    pyscr="""molecule load psf _ANNDIR_/toppar/crd.md18_vmd_autopsf.psf dcd _ANNDIR_/annealing.dcd
foreach temp {10 190 200 210 220 230 240 250 260 270 280 290} {
  set frame [expr 1+290-$temp]
  animate write pdb _EQDIR_/T$temp/init.pdb beg $frame end $frame
}
"""
    pyscr = pyscr.replace( '_ANNDIR_', kwargs['anndir'] )
    pyscr = pyscr.replace( '_EQDIR_', kwargs['eqdir'] )
    pt,scriptFile = mkstemp( dir='/tmp' )
    open(scriptFile, 'w').write( pyscr )
    print '\n#########################################################'
    print '\n Manually run "vmd -dispdev none -e %s" to extract the conformations\n' % scriptFile
    print ' #########################################################\n'
  elif JOB=='extract extended info from PDB header':
    """When using VMD to extract a frame to a PDB file from a DCD containing unit cell information,
    the cell information is stored in the first line of the PDB file with a "CRYST1" header.
    We use this to create an .xsc file."""
    XYZ = open(kwargs['inpdb']).readline().split()[1:4] # X, Y, Z cell dimensions
    xscLines = open(kwargs['inxsc']).readlines() # read input .xsc file
    sXYZ = xscLines[2].split()
    sXYZ[1] = XYZ[0]; sXYZ[5] = XYZ[1]; sXYZ[9] = XYZ[2] # substitute cell dimensions
    xscLines[2] = ' '.join( sXYZ )
    open( kwargs['outxsc'], 'w' ).write( ''.join( xscLines )+'\n' ) # save output .xsc file
  else:
    print 'Job not recognized. Nothing to do'

def perturbSystem( *kargs, **kwargs ):
  #joblist = ['change Hydrogen charge',]
  JOB = kwargs['job'] #what job to do
  if JOB=='change Hydrogen charge':
    """Load the PSF file and change the hydrogen and oxygen charges to whatever value
    is passed for the charge of the hydrogen. Save to another psf file"""
    Hq = float( kwargs['Hq'] )
    psf = open( kwargs['inpsf'] ).read()
    psf = psf.replace( ' 0.417000', '%9.6f'%Hq ) # 0.417 is the charge of the hydrogen
    psf = psf.replace( '-0.834000', '%9.6f'%(-2.*Hq,) ) # change the oxygen charge
    open( kwargs['outpsf'], 'w' ).write( psf )
  else:
    print 'Job not recognized. Nothing to do'

def analyze( *kargs, **kwargs ):
  #joblist = ['parse out file','find average','fit spectra at constant T', 'instant volume']
  JOB = kwargs['job'] #what job to do

  if JOB=='parse out file':
    """Read the .out file from a NAMD simulation and extract the "ETITLE:' lines to a file"""
    os.system( 'grep "ETITLE:" %s|head -1 > ./junk.txt; sleep 1s' % kwargs['inOutFile'] )
    items = open('./junk.txt').read()[8:].split()
    data = '#' + open('./junk.txt').read()[8:]
    for i in range( len(items) ):
      data = data.replace( items[i]+'    ', items[i] + '(%02d)' % (1+i,) )
    os.system( 'grep "ENERGY:" %s > %s; sleep 1s' % ( kwargs['inOutFile'], kwargs['outDatFile'] ) )
    for line in open( kwargs['outDatFile'] ).readlines():
      data += line[7:]
    open( kwargs['outDatFile'], 'w' ).write( data )
    #os.system('/bin/rm junk.txt') # clean-up

  elif JOB=='find average':
    """read input files and print the average for the given column"""
    from utilities.readingWritingFiles import read_column
    files = [ x.strip() for x in os.popen( 'ls %s' % kwargs['inFiles'] ).readlines() ]
    column = int( kwargs['column'] )
    values = [ ]
    for xfile in files:
      values += read_column( xfile, column, isFloat=1 )
    print '%7.2f' % numpy.average( numpy.array( values ) )

  elif JOB=='instant volume':
    """calculate and store the volume at each step, as given by the *.xst file"""
    from utilities.readingWritingFiles import read_to_cols
    from math import sqrt
    x=read_to_cols(kwargs['xstfile'], comment='#', separator=' ', multiple=True, outFmt='list',
                   outComm=False, xtype='float', nskip=3)
    buf='# time(ps) volume (A^3)\n'
    avv=0
    av2=0
    N=len(x[0])
    for i in range(N):
      v=x[1][i]*x[5][i]*x[9][i]
      avv+=v
      av2+=v*v
      buf+='%4.0f %8.2f\n'%(x[0][i], v) # Lx,Ly,Lz correspond to indexes 1, 5, 9
    open(kwargs['volfile'],'w').write(buf)
    print '%f %f'%(avv/N,sqrt(av2/N-(avv/N)**2))

  elif JOB=='elastic line versus simulation time':
    """S(Q,E=0) is the integral of I(Q,t). We vary the limits of integration [-tmax, tmax]
    to ascertain how S(Q,E=0) changes with tmax. For sufficiently large tmax, S(Q,E=0)
    should plateau."""
    from mantid.simpleapi import LoadSassena,SortByQVectors,ConvertToHistogram,Integration,mtd
    timeUnit=float(kwargs['timeUnit'])
    LoadSassena(Filename=kwargs['fqtfile'],timeUnit=timeUnit,OutputWorkspace='temp')
    SortByQVectors(InputWorkspace='temp')
    ws=ConvertToHistogram(InputWorkspace='temp_fqt.Re',OutputWorkspace='temp_fqt.Re')
    N=len(ws.readX(0))/2 # number of positive time points
    delta=int(N/100)
    buf='# tmax S(Q_1,E=0) S(Q_2,E=0) S(Q_3,E=0) ..\n'
    for it in range(delta,N,delta):
      tmax=timeUnit*it
      buf += ' %f'%tmax
      ws=Integration(InputWorkspace='temp_fqt.Re',RangeLower=-tmax,RangeUpper=tmax,OutputWorkspace='junk')
      for iq in range(ws.getNumberHistograms()): buf+= ' %f'%ws.readY(iq)[0]
      buf+='\n'
    open(kwargs['outfile'],'w').write(buf)

  elif JOB=='fit spectra at constant T':
    """fit Eugene data to simulation
    """
    def OuputExperimentalNexusFile(strT, expNexus):
      """ output experimental and spectra to a Nexus files"""
      from mantid.simpleapi import CloneWorkspace, AppendSpectra, SaveNexus
      Qvals=['3', '5', '7', '9']
      #trace()
      Qval=Qvals.pop(0)
      wsname='data_0_%s' %strT
      CloneWorkspace(InputWorkspace='data_0_%s_%s' % ( Qval, strT ), outputWorkspace=wsname)
      while Qvals:
        Qval=Qvals.pop(0)
        AppendSpectra(InputWorkspace1=wsname, InputWorkspace2='data_0_%s_%s' % ( Qval, strT ),OutputWorkspace=wsname)
      SaveNexus(InputWorkspace=wsname, Filename=expNexus)
    def LoadSassenaOutput( strT ):
      """load Sassena hdf5 files, perform Fourier transform"""
      from mantid.simpleapi import LoadSassena,SassenaFFT
      from mysassena.version import addVersionStamp
      for Hq in ( '32', '34', '36', '38', '40', '42', '44', '46', '48', '50' ):
        fname = kwargs['simDir']+'/Q%s/T%s/production/fqt.hd5' % ( Hq, strT )
        ws_name='Hq%s_T%s' % ( Hq, strT )
        addVersionStamp(fname,'1.4.1')
        LoadSassena(Filename=fname, OutputWorkspace =ws_name, TimeUnit=1.0)
        SassenaFFT(InputWorkspace=ws_name, FFTonlyRealPart=True, DetailedBalance=True, Temp=float(strT))
    def interpolSimW( strT, Qval ):
      """create interpolated workspaces in Hydrogen charge for a given Q-value"""
      from mantid.simpleapi import ExtractSingleSpectrum,CloneWorkspace,Scale,ConjoinWorkspaces,SumSpectra
      # First extract spectra with given Q value
      for Hq in ( '32', '34', '36', '38', '40', '42', '44', '46', '48', '50' ):
        ExtractSingleSpectrum(InputWorkspace = 'Hq%s_T%s_sqw' % ( Hq, strT ), OutputWorkspace = 'Q%s_Hq%s0_T%s_sqw' % ( Qval, Hq, strT ), WorkSpaceIndex = int(Qval) - 1 )
      # Now interpolate, create 9 spectra in between each spectra pair
      Hqs = ( '320', '340', '360', '380', '400', '420', '440', '460', '480', '500' )
      return Hqs #Do not interpolate
      Hqx = []
      for i in range( len(Hqs) - 1):
        Hqx.append( Hqs[i] )
        for k in range( 1, 10 ):
          dk = 2*k
          CloneWorkspace( InputWorkspace = 'Q%s_Hq%s_T%s_sqw' % ( Qval, Hqs[i], strT ), OutputWorkspace = 'junk0' )
          Scale( InputWorkspace = 'junk0', OutputWorkspace = 'junk0', Factor = 1-0.1*k, Operation = 'Multiply' )
          CloneWorkspace( InputWorkspace = 'Q%s_Hq%s_T%s_sqw' % ( Qval, Hqs[i+1], strT ), OutputWorkspace = 'junk1' )
          Scale( InputWorkspace = 'junk1', OutputWorkspace = 'junk1', Factor = 0.1*k, Operation = 'Multiply' )
          ConjoinWorkspaces( InputWorkspace1 = 'junk0', InputWorkspace2 = 'junk1', CheckOverlapping='0' )
          SumSpectra( InputWorkspace = 'junk0', OutputWorkspace = 'Q%s_Hq%3d_T%s_sqw' % ( Qval, int(Hqs[i])+dk, strT ) )
          Hqx.append( str( int(Hqs[i]) + dk ) )
      Hqx.append( Hqs[-1] )
      return Hqx  #all values of Hq
    def fitSimToDat( simW, datW, resW, model='model1' ):
      """fit simulation spectra to data spectra. The model is the simulation plus
      a background convoluted with the experimental resolution
      The experimental resolution is taken as the mirror image of resW
      model1: resolution x ( simulated + FlatBkg)
      model2: elastic_line + FlatBkg + resolution x simulated
      """
      from mantid.simpleapi import ScaleX,Fit,mtd,NormaliseToUnity,DeleteWorkspace
      # Normalize both simW and datW in the range to be fit
      start_x=-50.0
      end_x=50.0
      NormalizeSpectrum(simW,start_x,end_x)
      NormalizeSpectrum(datW,start_x,end_x)
      ScaleX(InputWorkspace=resW, OutputWorkspace='resolution',factor=-1)
      funcStr = ''
      if model == 'model1':
        funcStr  += '(composite=Convolution'
        funcStr += ';name=TabulatedFunction,FileName="",Workspace=%s,Scaling=1.0,ties=(Scaling=1)' % 'resolution'
        funcStr += ';(name=TabulatedFunction,FileName="",Workspace=%s' % ( simW )
        funcStr += ';name=FlatBackground,A0=0.001,constraints=(0<A0)))'
      elif model == 'model2':
        funcStr  += '(composite=Convolution'
        funcStr += ';name=TabulatedFunction,FileName="",Workspace=%s,Scaling=1,ties=(Scaling=1)' % 'resolution'
        funcStr += ';name=TabulatedFunction,FileName="",Workspace=%s,Scaling=1e-10)' % simW  #resolution x simulated
        funcStr += ';name=TabulatedFunction,FileName="",Workspace=%s,Scaling=1' % resW  #elastic line
        funcStr += ';name=FlatBackground,A0=0.001,constraints=(0<A0)' #flat background
      f=Fit( funcStr, datW, StartX=start_x, EndX=end_x, CreateOutput='1' )
      print f.getParameterNames()
      parW = mtd[ datW + '_Parameters' ] #retrieve the cost of the fitting
      printFittedParms(parW)
      DeleteWorkspace('resolution')  #clean-up
      return parW.column(1)[-1]
    def printFittedParms(parW):
      names=parW.column(0)
      values=parW.column(1)
      for i in range(len(names)):
        sys.stdout.write('%s = %f '%(names[i],values[i]))
      print ''

    from mantid.simpleapi import Scale,DeleteWorkspace
    strT=kwargs['Temperature']
    model='model1' #fitting model
    if 'model' in  kwargs.keys(): model=kwargs['model']
    LoadEugeneData( strT, kwargs['dataDir'])
    if 'expNexus' in kwargs.keys():
      OuputExperimentalNexusFile(strT, kwargs['expNexus'])  # output experimental data to file expNexus (Nexus format)
    LoadSassenaOutput( strT )
    cumChi2 = None
    strOut = '#T=%s\n#Q-val Hyd-q Chi2\n' % strT
    for Qval in ('3', '5', '7', '9'):
      datW = 'data_0_%s_%s' % ( Qval, strT )  #data
      resW = 'data_0_%s_006' % Qval  #resolution
      Hqx = interpolSimW( strT, Qval )
      if not cumChi2: cumChi2 = [0.] * len(Hqx) #initialize for the first time
      iHq = 0
      for Hq in Hqx:
        simW = 'Q%s_Hq%s_T%s_sqw' % ( Qval, Hq, strT )
        sys.stdout.write('Q=%s Hq=%s '%(Qval,Hq))
        Chi2 = fitSimToDat( simW, datW, resW, model=model )
        DeleteWorkspace( simW )
        for suffix in ( 'NormalisedCovarianceMatrix', 'Parameters', 'Workspace'):
          DeleteWorkspace( datW+'_'+suffix)
        #strOut += '0.%s 0.%s %6.3f\n' % ( Qval, Hq, Chi2 )
        #strOut += '0.%s %6.3f\n' % ( Hq, Chi2 )
        strOut += '0.%s %f\n' % ( Hq, Chi2 )
        cumChi2[ iHq ] += Chi2  #Chi2 for all Q-values
        iHq += 1
      #open( kwargs['outFile'], 'w' ).write( strOut)
      strOut += '&\n'
    open( kwargs['outFile'], 'w' ).write( strOut)
    for iHq in range( len(Hqx) ):
      print '0.%s %f' % ( Hqx[iHq], cumChi2[iHq] )

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
        sys.stderr.write('Error: error in '+service+' or service '+service+' not found \n')
        exitCode=1
    sys.exit(exitCode)

