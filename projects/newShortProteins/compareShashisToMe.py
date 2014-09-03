#!/usr/bin/python
"""
For each of the Shashi's top five models ranked by size, find my model
(out of the first five top models ranked by density) that has the
highest TM score.
"""
import sys,os
from utilities.small_utilities import chomp,junkName
from spicker.spickerYangResultsManager import spickOut,TMscore
from random import randint
from inputArgs.inputArgs import addAbsPath


#compare to mouse genes


#compare to benchmark
listf='/gpfs1/scratch/jose/newShortProteins/benchmark/success.list' #non-identical list of targets
#listf='/gpfs1/scratch/jose/newShortProteins/benchmark/toy.list'
myroot='/gpfs1/scratch/jose/newShortProteins/benchmark'
hisroot='/gpfs1/scratch/shashi/bench-100/Tasserout'

#junk directory where I will untar shashi's spicker results
junkd=junkName() ; os.system('/bin/mkdir -p '+junkd)

print '#header, shashi-size-rank, my-size-rank, my-dens-rank, tm'

headers=chomp(open(listf,'r').readlines())
for header in headers:
    tarred=hisroot+'/'+header+'.tasout.bz'
    if not os.path.exists(tarred): continue   #no shashi ouput for this target
    os.system('tar jxf '+tarred+' -C '+junkd) #untar to junk directory
    t=spickOut(dir=junkd)
    if t.readError: continue #error reading shashi output
    hisIDlist=t.ids[0:5] #first five clusters ranked by size (the default ranking in output spiker)
    mydir=myroot+'/'+header[1]+'/'+header #; print mydir
    s=spickOut(dir=mydir) 
    if s.readError: continue #error reading my output
    densranked=s.rankIDsByDens()[0:5] #top five clusters, ranked by density
    #compare each of my top five clusters to Shashi's top five clusters
    for id in hisIDlist:
        myBestDensID=0
        relatedSizeID=0
        tmMax=0.0
        for j in range(0,len(densranked)):
            mySizeID=densranked[j]
            tm=float(TMscore(s.combo[mySizeID],t.combo[id]))
            if tm>tmMax:
                tmMax=tm
                myBestSizeID=mySizeID
                relatedDensID=j+1
        print '%s %2d %2d %2d %5.3f'%(header,id,myBestSizeID,relatedDensID,tmMax)
        
    os.system('/bin/rm '+junkd+'/*') #empty the contents of junk directory
