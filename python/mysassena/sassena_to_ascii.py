#!/usr/bin/env python

import sys
import argparse
from pdb import set_trace as tr

sys.path.append('/opt/Mantid/bin')
from mantid.simpleapi import LoadSassena, mtd, SaveAscii

parser = argparse.ArgumentParser(description='''Translate a Sassena HDF5 file into a series of ASCII files. Example: python sassena_to_ascii.py fqt_coh.h5 fqt_coh 2.0''')
parser.add_argument('infile',help='input Sassena file')
parser.add_argument('basename',help='basename for the output ASCII files')
parser.add_argument('timeunit', type=float, help='time in between simlation frames, in picoseconds')
args = parser.parse_args()

LoadSassena(Filename=args.infile, TimeUnit=args.timeunit, SortByQVectors=1,OutputWorkspace='fqt')
for dataset in ('qvectors', 'fq', 'fq0', 'fq2', 'fqt.Re', 'fqt.Im'):
    wsname = 'fqt_' + dataset 
    if mtd.doesExist(wsname):
        outfile = args.basename + '_' + dataset + '.dat'
        SaveAscii(InputWorkspace=wsname, Filename=outfile, Separator='Space',
                  CommentIndicator='#', Version=1)

