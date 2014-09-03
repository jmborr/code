#!/usr/bin/pyton

""" for each sequence calculate average TM-score, coverage, and
rmsd_coverage, then average for all sequences that have best model
within some Z-rRMSD-to-native range

"""

import sys,os,re
from utilities.small_utilities import chomp
from spicker.spickerResultsManager import spickOut
from math import sqrt

pre=re.compile('tm=\s(\d\.\d+).+cov=\s(\d\.\d+)\srmsd_cov=\s+(\d+\.\d+)')

headers=chomp(open('/gpfs1/active/jose/projects/nonhomol_pdbJul05/small_proteins/successA.list','r').readlines())
#headers=chomp(open('/gpfs1/active/jose/projects/nonhomol_pdbJul05/small_proteins/toylist','r').readlines())
root='/gpfs1/active/jose/projects/nonhomol_pdbJul05/small_proteins'

avvTM=0.0  ;  avv2TM=0.0  ;  avvCov=0.0  ;  avv2Cov=0.0  ;  avvR=0.0  ; avv2R=0.0
nc=5
for header in headers:
    p=spickOut(dir=root+'/'+header[1]+'/'+header,nat='CA')
    avTM=0.0  ;  avCov=0.0  ;  avR=0.0
    ranks=p.rankIDsByDens()[0:nc] #rank centroids by density, select the first 'nc'
    zid=p.rankIDsByZrRMSDtoNat(list=ranks)[0]
    Zbest=p.zrrmsds[zid] #get best Zrmsd in 'ranks' list
    if Zbest>-1: #only average those significant predictions
        for id1 in ranks:  #only first five centroids
            for id2 in ranks:
                if id1==id2: continue
                l=p.TMverboseIdId(id1,id2,d0=8.0)
                [tm,cov,rmsd_cov]=pre.search(l).groups()
                tm=float(tm)  ;  cov=float(cov)  ;  rmsd_cov=float(rmsd_cov)
                avTM+=tm  ;  avCov+=cov  ;  avR+=rmsd_cov
        x=nc*(nc-1)
        avTM=avTM/x  ;  avCov=avCov/x  ;  avR=avR/x
        avvTM+=avTM  ;  avv2TM+=avTM*avTM
        avvCov+=avCov  ;  avv2Cov+=avCov*avCov
        avvR+=avR  ;  avv2R+=avR*avR
        print header,zid,Zbest,avTM,avCov,avR

avvTM=avvTM/len(headers)  ;  devTM=sqrt(avv2TM/len(headers) - avvTM*avvTM )
avvCov=avvCov/len(headers)  ;  devCov=sqrt(avv2Cov/len(headers) - avvCov*avvCov )
avvR=avvR/len(headers)  ;  devR=sqrt(avv2R/len(headers) - avvR*avvR )

print 'TM: ',avvTM,' +- ',devTM,' cov: ',cov,' +- ',devTM,'rmsdCov: ',avvR,' +- ',devR
