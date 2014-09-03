#!/usr/bin/python

import os,sys,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp
from tmScore.tm_score_cpp import getd0,copyD2,gettm,allocD2,print_array

N=2000 #maximum number of amino acids

#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: findTMScores.py [options]\nFor each template in chain.dat, do TM-Score to native file. Use d0 the length of the whole native file, in order to compare tm-scores of templates having different coverages',
               '  -a _AR_natf native pdb file',
               '  -b _AR_chainf chain.dat file')
    ih.parse(locals(),sys.argv)

    #store CA coordinates of native file
    natl=chomp(os.popen('cat '+natf+' |grep " CA "|grep "ATOM "').readlines())
    #store templates of chain.dat in a set of dictionaries
    ptf=open(chainf,'r') ; line=chomp(ptf.readline())
    templs=[]
    while line:
        items=line.split()
        x=''
        while len(items)==4:
            if not x:
                templs.append({}) #add one empty template entry
                templ=templs[-1]  #reference ot the created entry
                x='1'
            templ[int(items[0])]='%8.3f%8.3f%8.3f'% \
            (float(items[1]),float(items[2]),float(items[3]))
            line=chomp(ptf.readline())
            items=line.split()
        line=chomp(ptf.readline())
    ptf.close()
    if not templs[-1]: del templs[-1] #remove last template entry if empty

    mres=len(natl)
    d0=getd0(mres)     #find d0 for the full native length
    
    #for each template, retrieve corresponding atoms from native file, then do TM-score
    print 'naligned coverage tm'
    x=allocD2(N,3) #allocate double ** arrays
    y=allocD2(N,3)
    for templ in templs:
        natr=[] #list containing native coordinates
        templr=[]
        resnn=templ.keys()  ; resnn.sort() #residue numbers for aligned residues
        nres=len(resnn) #number of aligned residues in the template
        for resn in resnn:
            resi=resn-1
            for r in natl[resi][30:54].split(): natr.append(float(r)) # resi is index of natl list
            for r in templ[resn].split(): templr.append(float(r)) # resn is key of templ dict
        copyD2(nres,3,x,natr) #copy list into double** array
        copyD2(nres,3,y,templr)
        tm=gettm(x,y,nres); #calculate tm score
        print '%3d %4.2lf %4.2lf'%(nres,float(nres)/mres,tm)

#        extern double gettm(double **x, double **y, int nn, double d0=0., double *tx=NULL, double *ty=NULL, double **U=NULL, int Lmin=18, int n_chunks=10, double tmco=0.85);

sys.exit(0)
