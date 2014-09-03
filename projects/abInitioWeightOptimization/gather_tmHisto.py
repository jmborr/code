#!/usr/bin/python
'''
The typical tmHisto.dat file looks like:
 TM    N(TM)
0.550  6122
0.650 10000
0.750 10000
0.850 10000
0.950 10000
'''
import sys,os,re
from utilities.small_utilities import chomp
from math import sqrt

tms=['0.550','0.650','0.750','0.850','0.950']
av=[0.0]*5  #average bin
dev=[0.0]*5 #standard deviation of the bin

root='/gpfs1/scratch/jose/abInitioWeightOptimization/100to199'
headers=chomp(open(root+'/success_hardNoGeomViol.list','r').readlines())
for header in headers:
    f=root+'/'+header[1]+'/'+header+'/tmHisto.dat'
    #join all lines
    l=open(f,'r').readlines()
    for i in range(1,len(l)):
        h=float( l[i].split()[1] )
        av[i-1]+=h
        dev[i-1]+=h*h

print ' TM    <N(TM)> dev(N(TM))'

for i in range(len(av)):
    av[i]=int(av[i]/len(headers))
    dev[i]=int(sqrt(dev[i]/len(headers)-av[i]*av[i]))
    print '%s %5d %5d'%(tms[i],av[i],dev[i])
    
sys.exit(0)
