#!/usr/bin/python

import sys,os
from utilities.small_utilities import chomp
from random import randint
from inputArgs.inputArgs import inpHand,deglobb
from TMalign.TMalignManageResults import TMalignOut
from spicker.spickerResultsManager import TMscore
from math import sqrt
from utilities.codedir import CAdir

def printHeader():
    return  '#modelHeader, pdbBest, L(pdbBest), tmHit, Z_TM, gd, 1-gd, tmHit*(1-gd), (tmnat, tmHit-tmNat)'
#=================================================================
#this function goes through a list of PDB files and finds the best
#match to a model structure in PDB format. We can pass a header for
#the model (modelHeader), if there's a chance that the model is
#actually a file of the PDB library
def findBestPDBhitByTMalign(modelf,modelHeader='',pdbent='',pdblistf='',nat=''):
    if not modelHeader: modelHeader=os.path.basename(modelf)
    if not pdbent: pdbent=CAdir+'/xxxxx.pdb' #generic pdb file
    #create a list of pdb headers
    if not pdblistf: pdblist=chomp(os.popen('ls -1 '+CAdir+' | cut -d"." -f 1 ').readlines())
    else:  pdblist=chomp(open(pdblistf,'r').readlines())
    tmHit=0 #tm score of the best hit
    avtm=0  #average tm score
    avtm2=0 #average of the square of the tm score
    n=0     #number of TMalign's
    junk='junk'+`randint(0,999)` #some junk file to be over-written many times
    for pdbHeader in pdblist: 
        #print 'pdbHeader=',pdbHeader
        if modelHeader==pdbHeader: continue #avoid coincidence
        pdbf=deglobb(pdbent,pdbHeader)
        if not os.path.exists(pdbf): continue
        #print 'modelf=',modelf,'pdbf=',pdbf ;sys.exit(1)
        os.system('TMalign '+pdbf+' '+modelf+' > '+junk)
        s=TMalignOut(junk) #import TMalign results into object 's'
        tm=s.tm  #;print 'tm=',tm #; sys.exit(1) #tm-score
        n+=1 ;print n
        avtm+=tm
        avtm2+=tm*tm
        if tm>tmHit:
            tmHit=tm
            gd=s.gapDensity() #;print 'gd=',gd ;sys.exit(1)#gap density in the first structure
            pdbBest=pdbHeader
            L=int(chomp(os.popen('grep -e " CA " '+pdbf+'|wc -l').readline()))            
            print modelHeader,pdbBest,L,tmHit,gd,avtm,avtm2,n #; sys.exit(1)
    avtm=avtm/n #average tm score
    sigma=sqrt(avtm2/n-avtm*avtm) #standard deviation of the tm scores
    Z_TM=(tmHit-avtm)/sigma #Z-score of the best hit
    res='%s %s %4d %5.3f %6.1f %5.3f %5.3f %5.3f'%(modelHeader,pdbBest,L,tmHit,Z_TM,gd,1.0-gd,tmHit*(1.-gd))
    if nat:
        natf=deglobb(nat,modelHeader)      
        tmnat=float(TMscore(modelf,natf))
        res+=' %5.3f %6.3f'%(tmnat,tmHit-tmNat)
    #clean up
    os.system('/bin/rm '+junk)
    return res
#=================================================================


if __name__=='__main__':   
    inpHand('Usage: findBestPDBhit.py\nFind the pdb structure with maximum TMalign to the model\nstructure. We can pass a model file, one x-globbed generic\nmodel file plus a header, or  one x-globbed generic model\nfile plus a list of headers',
            ' -a _AR_model (possibly x-globbed) model in pdb format',
            ' -b __modelHeader model header (in place of -c)',
            ' -c _A_modelList list of model headers (in place of -b)',
            ' -d _A_pdbent x-globbed pdb structure (def: CAdir/xxxxx.pdb)',
            ' -e _A_pdblistf list of pdb headers (def:a list generated from CAdir)',
            ' -f _A_nat (x-globbed) native file'
            ).parse(locals(),sys.argv)

    if not modelHeader and not modelList: modelHeader=os.path.basename(model)
    if modelHeader: modelList=[modelHeader,] #create list of headers out of one header
    print printHeader()

    for modelHeader in modelList:
        modelf=deglobb(model,modelHeader)
        #we can pass empty pdbent and pdblistf, since they will be initialized.
        print findBestPDBhitByTMalign(modelf,modelHeader=modelHeader,pdbent=pdbent,pdblistf=pdblistf)
        
    #exit with success
    sys.exit(0)
