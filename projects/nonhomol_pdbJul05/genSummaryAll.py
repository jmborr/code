#!/usr/bin/python

#EXAMPLE:
#/gpfs1/active/jose/code/projects/nonhomol_pdbJul05/genSummaryAll.py -a /gpfs1/scratch/jose/createInputsTasser/x/xxxxx/abinitio  -b /gpfs1/scratch/jose/nonhomol_pdbJul05/alpha_beta/successAplusBBelow200.list -e summary.all

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp
from spicker.spickerResultsManager import spickOut

ih=inpHand('Usage: genSummaryAll.py [options]\nGiven a list of headers, and a xxxxx-globbed directory,\ngoes to each de-globbed directory and assumes it will find\nsummary.txt, centroids.pdb, closest.pdb and a native file.\nThen calculates rmsd,Z-rRMSD,TM for biggest and densest\nclusters. Then, out of the "top" most dense clusters, finds\nthe one with the smallest rmsd, the one with the smallest\nZ-rRMSD and the one with the biggest TM score. Generate top.pdb and best.pdb for top and best clusters.\n',
           ' -a _RA_ind directory where headers are (xxxxx globbing, def=curr dir)',
           ' -b _RA_list list of headers',
           ' -c __nat native filename, supposed to be in SAME DIRECTORY as each summary.txt (xxxxx globbing allowed, def=CA)',
           ' -d __topc number of top clusters (ranked by density) from which select the best rmsd, best Z-rRMSD and best TM (def=5)',
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
pf.write('#Header Biggest(densRank,rmsd,Z-rRMSD,TM) Densest(sizeRank,rmsd,Z-rRMSD,TM)  BestRmsd(densRank,SizeRank,rmsd) BestZ-rRMSD(densRank,sizeRank,Z-rRMSD) BestTM(densRank,sizeRank,TM)\n')

for header in chomp(open(list,'r').readlines()):
    dir=deglobb(ind,header)  #de-globb directory
    nat2=deglobb(nat,header) #de-globb native file name
    p=spickOut(dir=dir,nat=nat2,target=header) #instantiate a spickOut object
    if p.readError: #some error reading files, go to next header
        pf.write(header+' '+p.readError+'\n')
        continue   
    p.initTMtoNat()      #TMscore of each cluster centroid to native
    p.initZrRMSDtoNat()  #Z-rRMSD of each cluster centroid to native
    p.initRmsdsToNat()   #rmsd    of each cluster centroid to native
    
    rankdens=p.rankIDsByDens()[0:topc] #Cluster's ids by decreasing density, up to 'topc' clusters

    bestRmsdId=p.rankIDsByRmsdtoNat(list=rankdens)[0] #id in rankdens with smallest rmsd to native
    bestZrRMSDId=p.rankIDsByZrRMSDtoNat(list=rankdens)[0] #smallest Z-rRMSD to native
    bestTMId=p.rankIDsByTMtoNat(list=rankdens)[0]            #biggest TM to native

    #Write best.pdb and top.pdb to the de-globbed directory
    p.printCentroid(dir+'/top.pdb')               #print top centroid
    p.printCentroid(dir+'/best.pdb',id=bestTMId)  #print best centroid (by TM)

    #permissivity   
    if (1-slack)*p.rmsds[ p.densest ]   < p.rmsds[ bestRmsdId ]     : bestRmsdId  =p.densest
    if (1+slack)*p.zrrmsds[ p.densest ] < p.zrrmsds[ bestZrRMSDId ] : bestZrRMSDId=p.densest
    if (1+slack)*p.tms[ p.densest ]     > p.tms[ bestTMId ]         : bestTMId    =p.densest
    
    #summarize our calculations
    a0=p.densRanks[p.biggest] #rank of most populated cluster if ranked by decreasing density
    a1=p.rmsds[p.biggest]    #most populated cluster
    a2=p.zrrmsds[p.biggest]
    a3=p.tms[p.biggest]

    b0=p.densest             #densest cluster
    b1=p.rmsds[p.densest]
    b2=p.zrrmsds[p.densest]
    b3=p.tms[p.densest]

    c0=p.densRanks[bestRmsdId]
    c1=bestRmsdId        #cluster in the "topc" most dense clusters with smallest rmsd to native
    c2=p.rmsds[bestRmsdId]

    d0=p.densRanks[bestZrRMSDId]
    d1=bestZrRMSDId      #cluster in the "topc" most dense clusters with smallest Z-rRMSD to native
    d2=p.zrrmsds[bestZrRMSDId]

    e0=p.densRanks[bestTMId]
    e1=bestTMId          #cluster in the "topc" most dense clusters with biggest TM to native
    e2=p.tms[bestTMId]

    #output line
    pf.write(header+ ' %2d %5.2f %5.1f %4.2f %2d %5.2f %5.1f %4.2f %2d %2d %5.2f %2d %2d %5.1f %2d %2d %4.2f\n'%(a0,a1,a2,a3, b0,b1,b2,b3, c0,c1,c2, d0,d1,d2, e0,e1,e2 ))
    
pf.close()
