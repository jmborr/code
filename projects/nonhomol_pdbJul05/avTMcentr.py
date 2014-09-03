#!/usr/bin/pyton

"""
for each sequence calculate average TM-score between centroids,
then average for all sequences

"""

import sys,os
from utilities.small_utilities import chomp
from spicker.spickerResultsManager import spickOut
from math import sqrt

headers=chomp(open('/gpfs1/active/jose/projects/nonhomol_pdbJul05/small_proteins/successAplusB.list','r').readlines())
#headers=chomp(open('/gpfs1/active/jose/projects/nonhomol_pdbJul05/small_proteins/toylist','r').readlines())
root='/gpfs1/active/jose/projects/nonhomol_pdbJul05/small_proteins'

avv=0
avv2=0

for header in headers:
    p=spickOut(dir=root+'/'+header[1]+'/'+header)
    av=0.0
    for id1 in p.ids:
        for id2 in p.ids:
          av+=p.TMidId(id1,id2)
    av=(av-p.nc)/(p.nc*(p.nc-1))
    avv+=av
    avv2+=av*av
    print header,av

avv=avv/len(headers)
dev=sqrt(avv2/len(headers) - avv*avv )

print 'average=',avv,'deviation=',dev
