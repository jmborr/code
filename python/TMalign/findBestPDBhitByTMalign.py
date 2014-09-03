#!/usr/bin/python

import sys,os
from utilities.small_utilities import chomp,junkName,Bye
from random import randint
from inputArgs.inputArgs import inpHand,deglobb
from TMalign.TMalignManageResults import TMalignOut
from spicker.spickerResultsManager import TMscore
from math import sqrt
from utilities.codedir import CAdir

def printHeader():
    return  '#modelHeader(1), L(2) pdb(3), L(4), AlignedLength(5), TM-score(tm)(6), rmsd_aligned(7), Z_TM(8), gap-density(gd)(9), 1-gd(10), tm*(1-gd)(11), (tmnat(12), tm-tmNat(13))'
#=================================================================
#this function goes through a list of PDB files and finds the best
#match to a model structure in PDB format. We can pass a header for
#the model (modelHeader), if there's a chance that the model is
#actually a file of the PDB library
def findBestPDBhitByTMalign(modelf,modelHeader='',pdbent='',pdblistf='',pdbdir='',nat='',nouts=1):
    if not os.path.exists('./TMalign'):
        Bye('ERROR from findBestPDBhitByTMalign: please put TMalign executable in working dir')
    if not pdbdir: pdbdir=CAdir
    if not modelHeader: modelHeader=os.path.basename(modelf)
    if not pdbent: pdbent=pdbdir+'/xxxxx.pdb' #generic pdb file
    #create a list of pdb headers
    if not pdblistf: pdblist=chomp(os.popen('ls -1 '+pdbdir+' | cut -d"." -f 1 ').readlines())
    else:  pdblist=chomp(open(pdblistf,'r').readlines())

    L0=int(chomp(os.popen('grep -e " CA " '+modelf+'|wc -l').readline()))
    #calculate all TM-scores, and a couple of averages
    avtm=0  #average tm score
    avtm2=0 #average of the square of the tm score
    n=0     #number of TMalign's
    junk=junkName() #some junk file to be over-written many times
    tms={} #store all tm-scores
    other={}
    for pdbHeader in pdblist: 
        if modelHeader==pdbHeader: continue #avoid coincidence
        pdbf=deglobb(pdbent,pdbHeader)
        if not os.path.exists(pdbf): continue
        #print 'modelf=',modelf,'pdbf=',pdbf ;sys.exit(1)
        os.system('./TMalign '+pdbf+' '+modelf+' > '+junk)
        s=TMalignOut(junk) #import TMalign results into object 's'
        tm=s.tm
        tms[pdbHeader]=tm   #;print 'tm=',tm #; sys.exit(1) #tm-score
        rmsd=s.rmsd         #rmsd of aligned chunk
        AlignedLength=s.AlignedLength #aligned chunk
        gd=s.gapDensity()
        other[pdbHeader]={'rmsd':rmsd,'AlignedLength':AlignedLength,'gd':gd}
        n+=1
        avtm+=tm
        avtm2+=tm*tm
        
    avtm=avtm/n #average tm score
    sigma=sqrt(avtm2/n-avtm*avtm) #standard deviation of the tm scores

    #sort TM-scores from biggest to smallest
    items=[(tm,pdbHeader) for pdbHeader,tm in tms.items()]
    items.sort()    #from smallest to biggest
    items.reverse() #from biggest to smallest

    #write output for the first nouts hits
    if nouts=='all': nouts=n
    res=''
    for (tm,pdbHeader) in items[0:nouts]:
        pdbf=deglobb(pdbent,pdbHeader)
        L=int(chomp(os.popen('grep -e " CA " '+pdbf+'|wc -l').readline()))        
        Z_TM=(tm-avtm)/sigma #Z-score of the best hit
        rmsd=other[pdbHeader]['rmsd']
        AlignedLength=other[pdbHeader]['AlignedLength']
        gd=other[pdbHeader]['gd']
        res+='%s %3d %s %3d %3d %5.3lf %5.2lf %6.1lf %5.3lf %5.3lf %5.3lf\n'%(modelHeader,L0,pdbHeader,L,AlignedLength,tm,rmsd,Z_TM,gd,1.0-gd,tm*(1.-gd))
        if nat:
            natf=deglobb(nat,modelHeader)      
            tmnat=float(TMscore(modelf,natf))
            res+=' %5.3lf %6.3lf\n'%(tmnat,tm-tmNat)
            
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
            ' -f _A_nat (x-globbed) native file',
            ' -g _A_pdbdir directory where to find the pdb file library (def: utilities::codedir::CAdir)',
            ' -i __nouts number of hits to output (def:1, the best hit. NOTE: We can pass keyword "all" in order to store all of them )'
            ).parse(locals(),sys.argv)

    if not modelHeader and not modelList: modelHeader=os.path.basename(model)
    if modelHeader: modelList=[modelHeader,] #create list of headers out of one header
    if not pdbdir: pdbdir=CAdir
    if not nouts: nouts=1
    elif  nouts!='all': nouts=int(nouts)
    print printHeader()

    for modelHeader in modelList:
        modelf=deglobb(model,modelHeader)
        #we can pass empty pdbent and pdblistf, since they will be initialized.
        sys.stdout.write( findBestPDBhitByTMalign(modelf,modelHeader=modelHeader,pdbent=pdbent,pdblistf=pdblistf,pdbdir=pdbdir,nouts=nouts) )
        
    #exit with success
    sys.exit(0)
