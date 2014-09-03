#!/usr/bin/python

#EXAMPLE:
#/gpfs1/active/jose/code/projects/nonhomol_pdbJul05/genSummaryAll.py -a /gpfs1/active/jose/projects/createInputsTasser/x/xxxxx/abinitio  -b successAplusB.list -e summary.all

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp
from spicker.spickerYangResultsManager import spickOut

ih=inpHand('Usage: genSummaryAll.py [options]\nGiven a list of headers, and a xxxxx-globbed directory,\ngoes to each de-globbed directory and assumes it will find\nsummary.txt, centroids.pdb, closest.pdb and a native file.\nThen calculates rmsd,TM for biggest and densest\nclusters. Then, out of the "top" most dense clusters, finds\nthe one with the smallest rmsd, and the one with the biggest TM score\n',
           ' -a _RA_ind x-globbed directory where spicker results are (def=curr dir)',
           ' -b _RA_list list of headers',
           ' -c __nat native filename, (xxxxx globbing allowed, def=CA)',
           ' -d __topc number of top clusters (ranked by density) from which select the best rmsd, and best TM (def=5)',
           ' -e _A_outf output filename  (def=summary.all)',
           ' -f __slack permissivity between first and topc cluster(0<=slack<1, def=0.0)',
           )
ih.parse(locals(),sys.argv)

#defaults
if not ind: int=os.getcwd()
if not nat: nat='CA'
if not topc: topc=5
else: topc=int(topc)
if not outf: outf='./summary.all'
if not slack: slack=0.0
else: slack=float(slack)

pf=open(outf,'w')
pf.write('#Header Biggest(densRank,rmsd,TM) Densest(sizeRank,rmsd,TM)  BestRmsd(densRank,SizeRank,rmsd) BestTM(densRank,sizeRank,TM)\n')

for header in chomp(open(list,'r').readlines()):
    dir=deglobb(ind,header)  #de-globb directory
    nat2=deglobb(nat,header) #de-globb native file name
    p=spickOut(dir=dir,nat=nat2) #instantiate a spickOut object
    if p.readError: #some error reading files, go to next header
        pf.write(header+' '+p.readError+'\n')
        continue   
    p.initTMtoNat()      #TMscore of each cluster centroid to native
    p.initRmsdsToNat()   #rmsd    of each cluster centroid to native
    
    rankdens=p.rankIDsByDens()[0:topc] #Cluster's ids by decreasing density, up to 'topc' clusters

    bestRmsdId=p.rankIDsByRmsdtoNat(list=rankdens)[0] #id in rankdens with smallest rmsd to native
    bestTMId=p.rankIDsByTMtoNat(list=rankdens)[0]            #biggest TM to native

    #permissivity   
    if (1-slack)*p.rmsds[ p.densest ]   < p.rmsds[ bestRmsdId ]     : bestRmsdId  =p.densest
    if (1+slack)*p.tms[ p.densest ]     > p.tms[ bestTMId ]         : bestTMId    =p.densest
    
    #summarize our calculations
    a0=p.densRanks[p.biggest] #rank of most populated cluster if ranked by decreasing density
    a1=p.rmsds[p.biggest]    #most populated cluster
    a3=p.tms[p.biggest]

    b0=p.densest             #densest cluster
    b1=p.rmsds[p.densest]
    b3=p.tms[p.densest]

    c0=p.densRanks[bestRmsdId]
    c1=bestRmsdId        #cluster in the "topc" most dense clusters with smallest rmsd to native
    c2=p.rmsds[bestRmsdId]

    e0=p.densRanks[bestTMId]
    e1=bestTMId          #cluster in the "topc" most dense clusters with biggest TM to native
    e2=p.tms[bestTMId]

    #output line
    pf.write('%7s %2d %5.2f %4.2f %2d %5.2f %4.2f %2d %2d %5.2f %2d %2d %4.2f\n'%(header,a0,a1,a3, b0,b1,b3, c0,c1,c2, e0,e1,e2 ))
    
pf.close()
