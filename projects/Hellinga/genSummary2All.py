#!/usr/bin/python

#EXAMPLE:
#/gpfs1/active/jose/code/projects/nonhomol_pdbJul05/genSummaryAll.py -a /gpfs1/active/jose/projects/createInputsTasser/x/xxxxx/abinitio  -b successAplusB.list -e summary.all

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp,junkName,Bye
from spicker.spickerYangResultsManager import spickOut

ih=inpHand('Usage: genSummaryAll.py [options]\nGiven a list of headers, and a xxxxx-globbed directory,\ngoes to each de-globbed directory and assumes it will find\noutput from spicker and a native file.\nThen calculates rmsd,TM for biggest and densest\nclusters. Then, out of the "top" most dense clusters, finds\nthe one with the smallest rmsd, and the one with the biggest TM score\n',
           ' -a _RA_ind globbed directory where spicker results are (def=curr dir)',
           ' -b _RA_listf list of headers',
           ' -c _RA_nat globbed native file (def=./CA)',
           ' -d __topc number of top clusters (ranked by density) from which select the model, and best TM (def=5)',
           ' -e _A_outf output filename  (def=./summary.all)',
           ' -f __slack permissivity between topc and any more dense cluster(0<=slack<1, def=0.0)',
           ' -g __tarf0 globbed tar basename of the file containing all spiker results. If such a file exists, then it is assumed to be under the globbed directory where spicker results are , so pass only the basename (def:none)',
           ' -i __strtype do we use closc instead of combo? (def:yes)'
           )
ih.parse(locals(),sys.argv)

#defaults
if not ind: ind=os.getcwd()
if not nat: nat=os.path.join(ind,'CA')
if not topc: topc=5
else: topc=int(topc)
if not outf: outf='./summary.all'
if not slack: slack=0.0
else: slack=float(slack)
if tarf0: tarf0=os.path.basename(tarf0) #make sure we only pass the basename
if strtype and strtype in 'noNo': strtype='combo'
else: strtype='closc'

list=chomp(open(listf,'r').readlines())
nleft=len(list)
buf=''
for header in list:
    print header,'remaining=',nleft
    dir=deglobb(ind,header)  #de-globb directory
    dirT=dir #temporary directory...or may be not
    nat2=deglobb(nat,header) #; print 'nat2=',nat2 #de-globb native file name
    if tarf0:
        tarf=os.path.join(dir,deglobb(tarf0,header))
        dirT=os.path.join(dir,junkName()) #we will work in a temporary directory
        cmd='/bin/mkdir -p '+dirT+' && /bin/tar xf '+tarf+' -C '+dirT #; print cmd
        os.system(cmd)
    p=spickOut(dir=dirT,nat=nat2) #instantiate a spickOut object
    if p.readError: #some error reading files, go to next header
        pf.write(header+' '+p.readError+'\n')
        continue   
    p.initTMtoNat(str=strtype)  #TMscore of each cluster centroid to native
    p.initRmsdsToNat()          #rmsd    of each cluster centroid to native

    rankdens=p.rankIDsByDens()[0:topc] #Cluster's ids by decreasing density, up to 'topc' clusters

    bestTMId=p.rankIDsByTMtoNat(list=rankdens)[0]            #biggest TM to native
    
    #permissivity. Find densest cluster with TM-score above the lower
    #limit. We'll select this cluster as the model representative of
    #the sequence, rather than cluster with highest TM score
    modelId=bestTMId
    lowerLimit=(1-slack)*p.tms[modelId]
    rankmodelId=p.densRanks[modelId]
    for id in rankdens[0:rankmodelId]:
        if p.tms[id] > lowerLimit:
            modelId=id
            break

    #summarize our calculations
    a0=p.densRanks[p.biggest] #rank of most populated cluster if ranked by decreasing density
    a1=p.rmsds[p.biggest]    #most populated cluster
    a3=p.tms[p.biggest]
    a4=p.numberdens[p.biggest]
    a5=p.dens[p.biggest]

    b0=p.densest             #densest cluster
    b1=p.rmsds[p.densest]
    b3=p.tms[p.densest]
    b4=p.numberdens[p.densest]
    b5=p.dens[p.densest]
    
    c0=p.densRanks[modelId]
    c1=modelId        #cluster in the "topc" most dense clusters that we select as model
    c2=p.tms[modelId]
    c4=p.numberdens[modelId]
    c5=p.dens[modelId]
    
    e0=p.densRanks[bestTMId]
    e1=bestTMId          #cluster in the "topc" most dense clusters with biggest TM to native
    e2=p.tms[bestTMId]
    e4=p.numberdens[bestTMId]
    e5=p.dens[bestTMId]

    f0=p.mutual_avTM(str=strtype,list=rankdens)
    
    #store results line
    buf+='%7s %2d %4.2f %7.5lf %5.2f %4.2f  %2d %4.2f %7.5lf %5.2f %4.2f  %2d %2d %4.2f %7.5lf %5.2f  %2d %2d %4.2f %7.5lf %4.2f %4.2f\n'%(header,a0,a4,a5,a1,a3,  b0,b4,b5,b1,b3,  c0,c1,c4,c5,c2,  e0,e1,e4,e5,e2, f0 )
    if dirT!=dir: os.system('/bin/rm -r '+dirT)
    nleft-=1

#print results
open(outf,'w').write('#Header Biggest(densRank,numberdens,dens,rmsd,TM) Densest(sizeRank,numberdens,dens,rmsd,TM)     model(densRank,SizeRank,numberdens,dens,rmsd) BestTM(densRank,sizeRank,numberdens,dens,TM)  mutual_avTM\n'+buf)
