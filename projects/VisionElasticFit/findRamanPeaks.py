'''
Created on Feb 21, 2013

@author: jmborr
'''

import sys
sys.path.append('/home/jmborr/devel/mantidproject/mantid/Code/debug/bin')
from mantid.simpleapi import *
import argparse
import numpy
from pdb import set_trace as trace

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='guess peaks for the Raman spectrum')
  parser.add_argument('--ramanTxt',help='raman spectrum as a tab-separated text file')
  parser.add_argument('--outFile',help='output file')
  args = parser.parse_args()

  ws=LoadAscii(Filename=args.ramanTxt,Separator='Tab',Unit='Dimensionless')
  
  buf='#centre height width'
  # Scan by fine variations of FWHM
  for fwhm in range (4,40,2):
    peaks = FindPeaks( ws, PeaksList=ws.getName()+'_peaks', FWHM=fwhm, PeakFunction='Gaussian',HighBackground=True,BackgroundType='Linear')
    for iRow in range(0, peaks.rowCount()):
      buf += '%f %f %f\n'%(peaks.row( iRow )['centre'], peaks.row( iRow )['height'], peaks.row( iRow )['width'])
    peaks = FindPeaks( ws, PeaksList=ws.getName()+'_peaks', FWHM=fwhm, PeakFunction='Gaussian',HighBackground=True,BackgroundType='Linear')
    for iRow in range(0, peaks.rowCount()):
      buf += '%f %f %f\n'%(peaks.row( iRow )['centre'], peaks.row( iRow )['height'], peaks.row( iRow )['width'])

  # Scan by coarse variations of FWHM and peak position
  tofRange=[0, 3700]
  tofList = [];  delta = 10;  tof=tofRange[0]+delta; tofList.append( tof )
  while tof < tofRange[1]:
    tofList.append(tof);  tof += delta

  for fwhm in range (4,40,2):
    peaks = FindPeaks( ws, PeaksList=ws.getName()+'_peaks', FWHM=fwhm, PeakFunction='Gaussian',HighBackground=True,BackgroundType='Linear',PeakPositions=tofList)
    for iRow in range(0, peaks.rowCount()):
      buf += '%f %f %f\n'%(peaks.row( iRow )['centre'], peaks.row( iRow )['height'], peaks.row( iRow )['width'])
    peaks = FindPeaks( ws, PeaksList=ws.getName()+'_peaks', FWHM=fwhm, PeakFunction='Gaussian',HighBackground=False,BackgroundType='Linear',PeakPositions=tofList)
    for iRow in range(0, peaks.rowCount()):
      buf += '%f %f %f\n'%(peaks.row( iRow )['centre'], peaks.row( iRow )['height'], peaks.row( iRow )['width'])

  open(args.outFile,'w').write(buf)