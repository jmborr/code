#!/usr/bin/python

import sys,os,re
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand,deglobb
from random import randint
from math import sqrt
from spicker.spickerResultsManager import spickOut,TMscore

description='we pass a native pdb file and a model pdb file. Both files should contain "side group" (SG) atoms. Native file should also contain CA atoms.  For the native file and for the model, we calculate the contact map, with a contact defined between two SG atoms, where the contact distance is given by matrix DISFIT. We calculate for native and model how many residues in the termini do not have contacts, and compare. We create another native file without the dangling residues. Then look at the top most dense centroids and calculate best RMSD, Z-rRMSD and TM-score with respect to this cropped native file.'

headerf='/gpfs1/active/jose/projects/nonhomol_pdbJul05/analysis/success200.list'
natSG0='CA_SG'
modelSG0='bestTMcentroid.pdb_SG'
dir0='/gpfs1/active/jose/projects/createInputsTasser/x/xxxxx/abinitio'

#header for the output
print '#ALL true-dang false-dang false-dock, N-TER (same 3), C-TER (same 3), bestRMSD bestZrRMSD bestTM of cropped models to cropped native'
    
#read DISFIT creating a two-key dictionary dist2[aa1][aa2]
dist2={}
p=open('/gpfs1/active/jose/code/f77/tertiary/DISTFIT','r')
#list containing three-letter codes for amino acids
#     GLY   ALA   SER   CYS   VAL   THR   ILE   PRO...
aas=chomp(p.readline()).split()
N=len(aas) #; print N #number of codes
#read remaining lines. They contain one three-letter codes and then distances
#GLY  450   450   494   482   531   537   576...
for line in chomp(p.readlines()):
    list=line.split()
    aa1=list[0]
    vals=list[1:]
    dist2[aa1]={}
    for i in range(0,N):
        aa2=aas[i]
        x=float(vals[i])/100.0 #remember DISFIT distance units are 0.01Angstroms
        dist2[aa1][aa2]=x*x #we store the square of the distances
p.close()

#set up some global variables
junk='junk'+`randint(0,999)` #junk file
topc=5 #look in the first five clusters
INF=100000.0 # a very big number
headers=chomp(open(headerf,'r').readlines())

#cicle through all headers
for header in headers:
    
    dir=deglobb(dir0,header) #particular directory
    natSG=dir+'/'+natSG0     #particular native file name
    modelSG=dir+'/'+modelSG0 #particular model file name
    #print dir,natSG,modelSG ; sys.exit(1)

    #read " SG " atoms of native file native file and store in "junk" file
    os.system('grep -e " SG " '+natSG+' > '+junk)

    #read junk file and store SG amino acid codes and coordinates
    aas=[] ; x=[] ; y=[] ; z=[]
    N=0
    for line in open(junk,'r').readlines():
        aas.append(line[17:20])
        x.append(float(line[30:38]))
        y.append(float(line[38:46]))
        z.append(float(line[46:54]))
        N+=1

    #find number of dangling amino acids in the N terminal of the native state
    contact=0 ; n=0
    for i in range(0,N-4):
        aai=aas[i] ; xi=x[i] ; yi=y[i] ; zi=z[i]
        for j in range(i+4,N):
            dx=xi-x[j] ; dy=yi-y[j] ; dz=zi-z[j]
            #stop if i makes a contact with j
            if dx*dx+dy*dy+dz*dz < dist2[aai][ aas[j] ]:
                #print i+1,aai,j+1,aas[j],sqrt(dx*dx+dy*dy+dz*dz), sqrt(dist2[aai][ aas[j] ])
                contact=1
                break
        if contact==1:
            nterdanglingnat=n #;print 'nterdanglingnat=',nterdanglingnat
            break
        n+=1
    #find number of dangling amino acids in the C terminal of the native state
    contact=0 ; n=0
    span=range(4,N) ; span.reverse()
    for i in span:
        aai=aas[i] ; xi=x[i] ; yi=y[i] ; zi=z[i]
        for j in range(0,i-4):
            dx=xi-x[j] ; dy=yi-y[j] ; dz=zi-z[j]
            #stop if i makes a contact with j
            if dx*dx+dy*dy+dz*dz < dist2[aai][ aas[j] ]:
                contact=1
                break
        if contact==1:
            cterdanglingnat=n #;print 'cterdanglingnat=',cterdanglingnat
            break
        n+=1
    
    #read " SG " atoms of model file and store in "junk" file
    os.system('grep -e " SG " '+modelSG+' > '+junk)

    #read junk file and store codes and coordinates
    aas=[] ; x=[] ; y=[] ; z=[]
    N=0
    for line in open(junk,'r').readlines():
        aas.append(line[17:20])
        x.append(float(line[30:38]))
        y.append(float(line[38:46]))
        z.append(float(line[46:54]))
        N+=1

    #find number of dangling amino acids in the N terminal of the model
    contact=0 ; n=0
    for i in range(0,N-4):
        aai=aas[i] ; xi=x[i] ; yi=y[i] ; zi=z[i]
        for j in range(i+4,N):
            dx=xi-x[j] ; dy=yi-y[j] ; dz=zi-z[j]
            #stop if i makes a contact with j
            if dx*dx+dy*dy+dz*dz < dist2[aai][ aas[j] ]:
                contact=1
                break
        if contact==1:
            nterdanglingmodel=n #;print 'nterdanglingmodel=',nterdanglingmodel
            break
        n+=1
    #find number of dangling amino acids in the C terminal of the model
    contact=0 ; n=0
    span=range(4,N) ; span.reverse()
    for i in span:
        aai=aas[i] ; xi=x[i] ; yi=y[i] ; zi=z[i]
        for j in range(0,i-4):
            dx=xi-x[j] ; dy=yi-y[j] ; dz=zi-z[j]
            #stop if i makes a contact with j
            if dx*dx+dy*dy+dz*dz < dist2[aai][ aas[j] ]:
                contact=1
                break
        if contact==1:
            cterdanglingmodel=n #;print 'cterdanglingmodel=',cterdanglingmodel
            break
    n+=1

    #calculate true-dangling, false-dangling, and false-docked for N
    #terminal, C terminal and the sum of the two
    d=min([nterdanglingnat,nterdanglingmodel])
    e=max([0,nterdanglingmodel-nterdanglingnat])
    f=max([0,nterdanglingnat-nterdanglingmodel])

    g=min([cterdanglingnat,cterdanglingmodel])
    h=max([0,cterdanglingmodel-cterdanglingnat])
    i=max([0,cterdanglingnat-cterdanglingmodel])

    #sum results for N-terminal and C-terminal
    a=d+g
    b=e+h
    c=f+i

    #print/store these results
    results='%2d %2d %2d %2d %2d %2d %2d %2d %2d'%(a,b,c,d,e,f,g,h,i)

    #we store now nondangling residues of native files as a cropped native
    #file, only CA atoms read " CA " atoms of native file and store in
    #memory
    cas=os.popen('grep -e " CA " '+natSG).readlines()
    dirname=chomp(os.popen('dirname '+natSG).readline())
    cacrop=dirname+'/CAcrop'
    caall=dirname+'/CA'
    p=open(cacrop,'w')
    for i in range(0+nterdanglingnat,len(cas)-cterdanglingnat): p.write(cas[i])
    p.close()

    #find best rmsd, Z-rRMSD, and TM-score with the cropped CA file
    s=spickOut(dir=dir,target=header) #load spicker's results
    rankdens=s.rankIDsByDens()[0:topc] #Cluster's ids by decreasing density, up to 'topc' clusters
    bestTM=0 ; bestRMSD=INF ; bestZrRMSD=INF #initialize
    segment=[1+nterdanglingnat,len(cas)-cterdanglingnat]
    for id in rankdens:
        s.printCentroid(junk,id=id,segment=segment)                    #print cropped centroid
        RMSD=float(chomp(os.popen('rmsd.x '+cacrop+' '+junk).readline())) #rmsd of model to native
        ZrRMSD=float(chomp(os.popen('zrrmsd.x '+cacrop+' '+junk).readline()))
        TM=float(TMscore(cacrop,junk))
        if RMSD < bestRMSD: bestRMSD=RMSD
        if ZrRMSD < bestZrRMSD: bestZrRMSD=ZrRMSD
        if TM > bestTM: bestTM=TM

    #find best rmsd, Z-rRMSD, and TM-score with the intact CA file
    bestTMall=0 ; bestRMSDall=INF ; bestZrRMSDall=INF #initialize
    segment=[1,len(cas)]
    for id in rankdens:
        s.printCentroid(junk,id=id,segment=segment)                    #print cropped centroid
        RMSD=float(chomp(os.popen('rmsd.x '+caall+' '+junk).readline())) #rmsd of model to native
        ZrRMSD=float(chomp(os.popen('zrrmsd.x '+caall+' '+junk).readline()))
        TM=float(TMscore(caall,junk))
        if RMSD < bestRMSDall: bestRMSDall=RMSD
        if ZrRMSD < bestZrRMSDall: bestZrRMSDall=ZrRMSD
        if TM > bestTMall: bestTMall=TM

    #print results
    print header+'  '+results+' %5.2f %6.2f %5.3f %5.2f %6.2f %5.3f'%(bestRMSD,bestZrRMSD,bestTM,bestRMSDall,bestZrRMSDall,bestTMall)
    
#cleaning
os.system('/bin/rm '+junk)

#exit with success
sys.exit(0)
