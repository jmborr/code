#!/usr/bin/python
"""

Out of the Shashi's top five clusters (ranked by size) and out of my
top five clusters (ranked by density), find the pair (S*,J*) having
the highest TM score.

Out of my top five clusters, find the cluster (J@) having highest TM
score to Shashi's top cluster (Stop)

Do TM score of S* to native, and Stop to native. The point is to make
a scatter plot of TM(S*,nat) versus TM(J*,S*) to detect any
correlation, and also the scatter plot of TM(Stop,nat) versus
TM(J@,Stop), also to detect any putative correlation

Calculate avtm, average TM score among my five top clusters

Obtain Jtop and Jbest, top and best cluster to native
"""
import sys,os
from utilities.small_utilities import chomp,junkName
from spicker.spickerYangResultsManager import spickOut,TMscore
from random import randint
from inputArgs.inputArgs import addAbsPath


#compare to mouse genes


#compare to benchmark
#listf='/gpfs1/scratch/jose/newShortProteins/benchmark/success.list' #non-identical list of targets
#listf='/gpfs1/scratch/jose/newShortProteins/benchmark/toy.list'

#this list is subset of hard targets for which prospector can't find a template with TM score to
#native above 0.35
listf='/gpfs1/scratch/jose/newShortProteins/benchmark/filter.lowtm'

myroot='/gpfs1/scratch/jose/newShortProteins/benchmark'
hisroot='/gpfs1/scratch/shashi/bench-100/Tasserout'
pdbdir='/gpfs1/active/adrian/short_prot/pdb'
#outd='/gpfs1/scratch/jose/newShortProteins/benchmark/very_hard'
#junk directory where I will untar shashi's spicker results
junkd=junkName() ; os.system('/bin/mkdir -p '+junkd)

print '# ./compareShashiToNat.py > compareShashiToNat.benchmark.dat'

print '#Out of the Shashi\'s top five clusters (ranked by size) and out of my top five clusters (ranked by density), find the pair (S*,J*) having the highest TM score'

print '#Out of my top five clusters, find the cluster (J@) having highest TM score to Shashi\'s top cluster (Stop)'

print '#avtm: average TM score among my top five clusters, ranked by density'

print '#Sbest: Shashi\'s best cluster to native, Jbest: my best cluster to native'

print '#header, Stop, J@-ranked-by-dens, J@-ranked-by-size, TM(Stop,J@), TM(Stop,nat), S*, J*-ranked-by-dens, J*-ranked-by-size, TM(S*,J*), TM(S*,nat), avTM, avTMtoTop, jtopBySize, jtopByDens, TM(Jtop,nat), density(Jtop),jbestBySize, jbestByDens, TM(Jbest,nat), Sbest, TM(Sbest,nat)'

headers=chomp(open(listf,'r').readlines())
for header in headers:
    nat=pdbdir+'/'+header+'.pdb' #native file
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
    myBestSizeIDtoStop=0  # J@ ranked by size
    relatedDensIDtoStop=0 # J@ ranked by density
    tmMaxToStop=0.0       # TM(Stop,J@)
    tmToNatToStop=0.0     # TM(Stop,nat)
        
    histBestSizeID=0     # S*
    myBestSizeID=0       # J* ranked by size
    relatedDensID=0      # J* ranked by density
    tmMax=0.0            # TM(S*,J*)        
    tmToNat=0.0          # TM(S*,nat)       

    sbest=0              # Sbest
    tmSbestToNat=0.0     # TM(Sbest,nat)
    
    jtopBySize=1         # Jtop ranked by size
    jtopByDens=0         # Jtop ranked by density
    tmJtopToNat=0.0      # TM(Jtop,nat)
    densJtop=0.0         # density(Jtop)
    jbestBySize=0        # Jbest ranked by size
    jbestByDens=0        # Jbest ranked by density
    tmJbestToNat=0.0     # TM(Jbest,nat)
    
    for id in hisIDlist: #go cluster by cluster in Shashi's top five clusters
        
        for j in range(0,len(densranked)): #go cluster by cluster in my first top clusters
            mySizeID=densranked[j]
            tm=float(TMscore(t.combo[id],s.combo[mySizeID]))
            if tm>tmMax:
                histBestSizeID=id
                myBestSizeID=mySizeID
                relatedDensID=j+1
                tmMax=tm

        if id==1:
            myBestSizeIDtoStop=myBestSizeID
            relatedDensIDtoStop=relatedDensID
            tmMaxToStop=tmMax
            tmToNatToStop=float(TMscore(t.combo[id],nat))

        tm=float(TMscore(t.combo[id],nat))
        if tm>tmSbestToNat:
            sbest=id
            tmSbestToNat=tm
            
    tmToNat=float(TMscore(t.combo[histBestSizeID],nat))

    jtopByDens=densranked[0]
    densJtop=s.dens[s.densest]
    tmJtopToNat=float(TMscore(s.combo[1],nat))

    #Obtain Jtop and Jbest
    for i in range(0,len(densranked)):
        id=densranked[i]
        tm=float(TMscore(s.combo[id],nat))
        if tm>tmJbestToNat:
            jbestBySize=id
            jbestByDens=i+1
            tmJbestToNat=tm
    
    #calculate average TM score among my top five clusters
    avtm=0.0
    n=0
    for i in range(0,len(densranked)-1):
        id1=densranked[i]
        for j in range(i+1,len(densranked)):
            id2=densranked[j]
            avtm+=float(TMscore(s.combo[id1],s.combo[id2]))
            n+=1
    avtm=avtm/n

    #calculate average TM score to the top cluster
    id1=densranked[0]
    avTMtoTop=0.0
    for i in range(1,len(densranked)):
        id2=densranked[i]
        avTMtoTop+=float(TMscore(s.combo[id1],s.combo[id2]))
    avTMtoTop=avTMtoTop/(len(densranked)-1)

    #os.system('/bin/cp '+s.combo[jbestBySize]+' '+outd+'/'+header+'/best_combo_abinitio.pdb')
    #os.system('/bin/cp '+t.combo[sbest]+' '+outd+'/'+header+'/best_combo_normal.pdb')

    print '%s %2d %2d %2d %5.3f %5.3f %2d %2d %2d %5.3f %5.3f %5.3f %5.3f %2d %2d %5.3f %6.4f %2d %2d %5.3f %2d %5.3f'%(header,1,relatedDensIDtoStop,myBestSizeIDtoStop,tmMaxToStop,tmToNatToStop,histBestSizeID,relatedDensID,myBestSizeID,tmMax,tmToNat,avtm,avTMtoTop,jtopBySize,jtopByDens,tmJtopToNat,densJtop,jbestBySize,jbestByDens,tmJbestToNat,sbest,tmSbestToNat)

    os.system('/bin/rm '+junkd+'/*') #empty the contents of junk directory
