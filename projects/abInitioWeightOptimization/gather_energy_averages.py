#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp

root='/gpfs1/scratch/jose/abInitioWeightOptimization/100to199'
headers=chomp(open(root+'/success_hardNoGeomViol.list','r').readlines())
for header in headers:
    f=root+'/'+header[1]+'/'+header+'/energy_averages.dat'
    #join all lines
    print header+' '.join( chomp(os.popen('grep -v "#" '+f).readlines()) )
    
sys.exit(0)
