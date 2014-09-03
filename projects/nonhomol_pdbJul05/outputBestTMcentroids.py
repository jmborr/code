#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp
from spicker.spickerResultsManager import spickOut

ih=inpHand('Usage: outputBestTMcentroids.py [options]\nGiven a list of headers, and a xxxxx-globbed directory,\ngoes to each de-globbed directory and assumes it will find\nsummary.txt, centroids.pdb, closest.pdb and a native file.\nThen finds centroid with best TM score to native, and outputs it as file "bestTMcentroid.pdb" in the de-globbed directory\n',
           ' -a _RA_ind (xxxxx-globbed) directory where summary.txt is (def=curr dir)',
           ' -b _RA_list list of headers',
           ' -c __nat native filename, supposed to be in SAME DIRECTORY as each summary.txt (xxxxx globbing allowed, def=CA)',
           ' -d __topc number of top clusters (ranked by density) from which select the best TM (def=5)',
           ' -e __seqdg (xxxxx-globbed) seq.dat file, in case we want TM-score for only non-coil assigned residues'
           )
ih.parse(locals(),sys.argv)

#defaults
if not ind: int=os.getcwd()
if not nat: nat='CA'
if not topc: topc=5
else: topc=int(topc)

for header in chomp(open(list,'r').readlines()):
    dir=deglobb(ind,header)  #de-globb directory
    nat2=deglobb(nat,header) #de-globb native file name
    if seqdg: seqd=deglobb(seqdg,header)
    p=spickOut(dir=dir,nat=nat2,target=header) #instantiate a spickOut object
    if p.readError: #some error reading files, go to next header
        sys.stderr.write(header+' '+p.readError+'\n')
        continue
    rankdens=p.rankIDsByDens()[0:topc] #Cluster's ids by decreasing density, up to 'topc' clusters
    if seqdg:
        bestTMId=p.rankIDsByTMonlySecStrToNat(seqd,list=rankdens)[0]
        print header+' %5.3f'%(p.TMonlySecStrToNat(seqd,id=bestTMId))
    else:
        p.initTMtoNat()      #TMscore of each cluster centroid to native
        bestTMId=p.rankIDsByTMtoNat(list=rankdens)[0]            #biggest TM to native
        print header+' %5.3f'%(p.tms[bestTMId])
        p.printCentroid(dir+"/bestTMcentroid.pdb",id=bestTMId)
