#!/usr/bin/python
"""
compare my first five top models ranked by density, against Shashi's top five models ranked by size
"""
import sys,os
from utilities.small_utilities import chomp,junkName
from spicker.spickerYangResultsManager import spickOut,TMscore
from random import randint
from inputArgs.inputArgs import addAbsPath


#compare to mouse genes
#listf='/gpfs1/scratch/jose/newShortProteins/0002.nonredlist' #non-identical list of targets
#myroot='/gpfs1/scratch/jose/newShortProteins/tasser'
#hisroot='/gpfs1/scratch/shashi/genomes/0002'

#compare to benchmark
listf='/gpfs1/scratch/jose/newShortProteins/benchmark/success.list' #non-identical list of targets
#listf='/gpfs1/scratch/jose/newShortProteins/benchmark/toy.list'
#listf='/gpfs1/scratch/jose/newShortProteins/benchmark/filter.lowtm'
myroot='/gpfs1/scratch/jose/newShortProteins/benchmark'
hisroot='/gpfs1/scratch/shashi/bench-100/Tasserout'

#junk directory where I will untar shashi's spicker results
junkd=junkName() ; os.system('/bin/mkdir -p '+junkd)

print '#header, my-size-rank, my-dens-rank, shashi-size-rank, tm'

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
    for i in range(0,len(densranked)):
        id=densranked[i]
        histBestID=0
        tmMax=0.0
        for hisID in hisIDlist:
            tm=float(TMscore(t.combo[hisID],s.combo[id]))
            if tm>tmMax:
                tmMax=tm
                histBestID=hisID
                myBestID=id
        print '%s %2d %2d %2d %5.3f'%(header,myBestID,i+1,histBestID,tmMax)
        
    os.system('/bin/rm '+junkd+'/*') #empty the contents of junk directory
