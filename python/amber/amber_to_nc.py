#!/usr/bin/python

# Conversion of a DCD trajectory (CHARMM/XPlor) to MMTK's trajectory format.
#
# Usage: dcd_to_nc pdb_file dcd_file nc_file
#

import getopt,os,string,sys,Numeric
import pdb
from MMTK import Units
from MMTK.ParticleProperties import ParticleVector
from MMTK.Universe import InfiniteUniverse,OrthorhombicPeriodicUniverse
from MMTK.PDB import PDBConfiguration
from MMTK.Trajectory import Trajectory
from amber.amber2MMTK import getBox,fixPDB,amberReader
from utilities.small_utilities import Bye

usage = """Usage: amber_to_nc prmpdb prmcrd nc_file [options]

prmpdb and prmcrd must describe the same system. nc_file is the
name of the netCDF trajectory that is created.

Options:

  --box=size
     specify a periodic system of a given size, in case prmtop and prmcrd
     lack the box info. The size can be a single number (cubic box) or
     three comma-separated values (orthorhombic box); in both cases
     the dimensions must be in Angstroms.

  --dt=dtime
     specify time between consecutive frames, since this parameter is
     missing in prmcrd files
     
  --vel=prmvel
     specify a velocity trayectory file

  --format=crd,cdf
     trajectory format (crd:ASCII amber (default); cdf: NetCDF)

  --block-size=number
     specifies the block structure of the netCDF trajectory. The
     default value of 1 optimizes the trajectory for step-by-step access
     to conformations. Larger values favour atom-by-atom access to
     one-particle trajectories for all times, which is required for the
     calculation of dynamic quantities. The highest sensible value is
     the number of steps in the trajectory.

  --pdbfixed=pdbfixed
    specify the name of the fixes to prmpdb in order to be
    successfully read by MMTK.
    Default: prmpdb+'.fixed'

  --pdbtest=pdbtest
    after MMTK Universe create, we output a PDB file. Maybe the
    numbering of atoms will be different than that of prmpdb.
    Default: prmpdb+'.mmtk'
"""


def check_options(format,opts):
  keys=opts.keys()
  buf=''
  if format=='crd':
    mesg='ERROR: %s needed if crd format\n'
    for req in ('dtime',):
      if not opts[req]: buf+=mesg%req
  if buf:
    sys.stderr.write(buf+'\n'+usage)
    sys.exit(1)
  
def check_files(prmpdb,prmcrd,nc_file):
  if not os.path.exists(prmpdb):
    sys.stderr.write('File %s not found\n' % prmpdb)
    raise SystemExit,1
  if not os.path.exists(prmcrd):
    sys.stderr.write('File %s not found\n' % prmcrd)
    raise SystemExit,1
  if os.path.exists(nc_file):
    sys.stderr.write('File %s already exists. ' % nc_file)
    while 1:
      answer = raw_input('Overwrite? [y/n] ')
      if answer.lower() in ('n','no'):
        raise SystemExit,0
      if answer.lower() in ('y','yes'):
        break

def convertAmberToMMTK(prmpdb,prmcrd,nc_file,
                       prmvel=None,format='crd',dtime=None,box=None,
                       pdbfixed=None,pdbtest=None,block_size=1):
  """box: list of box sizes, in Angstroms"""
  check_files(prmpdb,prmcrd,nc_file)

  #find number atoms in system, adapt the PDB file to MMTK requirem.
  pdbfixed,fixes,natoms=fixPDB(prmpdb,pdbfixed=pdbfixed)
  if fixes:
    print 'Adapted PDB file to MMTK convention with fixes:\n'+\
          '\n'.join([' * '+fix for fix in fixes])

  #find periodic boundary info
  if not box:
    box=getBox(prmcrd,natoms,format=format)
  try:
    box=Numeric.array(list(box),Numeric.Float)*Units.Ang
  except:
    sys.stderr.write('ERROR: could not find Box info. Exiting...\n')
    sys.exit(1)

  #Create the universe.
  options={'prmvel':prmvel,'format':format,'dtime':dtime,'box':box}
  check_options(format,options)
  universe = OrthorhombicPeriodicUniverse(box)

  # Create all objects from the PDB file. The PDB file *must* match the
  # AMBER trajectory file (same atom order), and you *must* create all
  # objects defined in the PDB file, otherwise interpretation of the
  # AMBER trajectory file is not possible.

  universe.addObject( PDBConfiguration(pdbfixed).createAll() )
  if prmvel:
    v0=Numeric.zeros((natoms,3),Numeric.Float)
    vels=ParticleVector(universe,data_array=v0)
    universe.setVelocities(vels)

  if not pdbtest: pdbtest=prmpdb+'.test'
  universe.writeToFile(pdbtest)
  
  #os.system('/bin/rm '+pdbfixed)

  # Create the trajectory object for output.
  files=[os.path.basename(file) for file in\
         (prmpdb,prmcrd,prmvel) if file ]
  convmsg='Converted from '+','.join(files)
  t=Trajectory(universe,nc_file,"w",convmsg,block_size=block_size)
  readAmber=amberReader(t,universe,natoms,prmcrd,**options)
  readAmber()
  universe.writeToFile(pdbtest.replace('.pdb','.finalconf.pdb'))
  t.close()


if __name__=='__main__':
  # Parse options
  try:
    options_keys=['box=','dt=','vel=','format=','block-size=']
    options,file_args=getopt.getopt(sys.argv[4:],'',options_keys)
  except getopt.GetoptError:
    sys.stderr.write(usage)
    raise SystemExit, 1
  box=None
  dtime=None
  prmvel=None
  format='crd' ; formats=('crd','cdf')
  block_size = 1
  pdbfixed=None
  pdbtest=None
  for option, value in options:
    if option == '--box':
      sizes = map(float, string.split(value, ','))
      if len(sizes) != 1 and len(sizes) != 3:
        sys.stderr.write(usage)
        raise SystemExit,1
      if len(sizes)==1: sizes=[sizes[0],]*3
      if option == '--dt':
        dtime=float(value)
      if option == '--vel':
        if not os.path.exists(value):
          sys.stderr.write('File %s not found\n' % prmvel)
          raise SystemExit,1
        prmvel=value
      if option == '--format':
        if value not in formats:
          sys.stderr.write('Specify a valid format ('+\
                           ','.join(formats)+')')
          raise SystemExit,1
        format=value
      if option == '--block-size':
        block_size = int(value)
        if block_size < 1:
          sys.stderr.write("Block size must be positive.")
          raise SystemExit,1
      if option=='--pdbfixed':
        if not os.path.exists(value):
          sys.stderr.write('File %s not found\n' % pdbfixed)
          raise SystemExit,1
        prmvel=value
      if option=='--pdbtest':
        if not os.path.exists(value):
          sys.stderr.write('File %s not found\n' % pdbtest)
          raise SystemExit,1
        prmvel=value
        
  # Check file arguments
  req_keys=['prmpdb','prmcrd','nc_file']
  if len(sys.argv)<1+len(req_keys):
    sys.stderr.write(usage)
    raise SystemExit, 1

  prmpdb, prmcrd, nc_file = sys.argv[1:4]
  convertAmberToMMTK(prmpdb,prmcrd,nc_file,
                     prmvel=prmvel,format='crd',dtime=dtime,box=box,
                     pdbfixed=pdbfixed,pdbtest=pdbtest,
                     block_size=block_size)

  sys.exit(0)
