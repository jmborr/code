#!/home/cs/calligari/bin/pyt
#
# This script converts generic AMBER9 netCDF output files in MMTK netCDF file
# for analysis with nMoldyn code.
#
# Corrected the bug on units for periodic box
# Written by: Paolo Calligari
#
# Last Revision: 27/02/2008
#

from Numeric import *
from MMTK import *
from MMTK.ParticleProperties import Configuration
from Scientific.IO.NetCDF import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import Protein
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
from MMTK.ForceFields import LennardJonesForceField
from MMTK.Universe import OrthorhombicPeriodicUniverse
from sys import argv, exit
from Scientific.Geometry import Vector
from string import *
 
pdbfile = argv[1]
ncin = argv[2]

#Open the netCDF input file
fil = NetCDFFile(ncin, 'r')

ncout = ncin[:-3] + '_filtered' + ncin[-3:]
parts = string.split(ncin, '_')
comment = 'From AMBER9 NetCDF file'

print comment

try:
    box= fil.variables['cell_lengths']
    box_array= box[:]*Units.Ang
    lbox_x=box_array[0][0]
    lbox_y=box_array[0][1]
    lbox_z=box_array[0][2]
    ang= fil.variables['cell_angles']
#    if ang[0][0] == 90.0 and ang[0][1] == 90.0:
#    universe = OrthorhombicPeriodicUniverse((lbox_x,lbox_y,lbox_z),
#                                        LennardJonesForceField(1.*Units.Ang))
except:
    universe=InfiniteUniverse()

# Total time for which trajectory undergoes the translation
Nsteps=len(fil.variables['time'])

# Get data about the universe.
conf=PDBConfiguration(pdbfile)
molecules=conf.createAll()
universe.addObject(molecules)

trajectory=Trajectory(universe, ncout, 'w', comment)
snapshot = SnapshotGenerator(universe, actions=[TrajectoryOutput(trajectory, ['all'], 0, None, 1)])

num=len(molecules.atomList())

for i in range(Nsteps):
    t=fil.variables['time'][i]
    conf_array= fil.variables['coordinates'][i][:num]*Units.Ang
    config=Configuration(universe, conf_array)
    conf = universe.contiguousObjectConfiguration(None,
                                                  config)
    universe.setConfiguration(conf)
    snapshot(data = {'time': t})

trajectory.close()
