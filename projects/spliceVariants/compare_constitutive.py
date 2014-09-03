#!/usr/bin/python

"""Compare tm-score of the constitutive regions for 'easy' predictions."""

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp,junkName,Bye
from utilities.codedir import scratchdir,nodecodedir,codedir
from seq.alignYangManageResults import alignOut

#global variables
rootd_fasta=os.path.join(scratchdir,'spliceVariants/input') #directory containing fasta files
algf=junkName()                              #temporary file to store sequence alignments
alignx=os.path.join(codedir,'f77/seq/align') #align program
tmScorex=os.path.join(nodecodedir,'python/tmScore/tm_score.py')   #TM-score program

if __name__=='__main__':
    inpHand('Usage: compare_constitutive.py [options]',
            '  -a _AR_summf summary.all file',
            '  -b _A_rootd_out rood dir for output (def:scratchdir/out3)',
            '  -c _A_outf output file (def:STDOUT)',
            ).parse(locals(),sys.argv)

    #optional input
    if not rootd_out: rootd_out=os.path.join(scratchdir,'spliceVariants/out3')
    if not outf: outp=sys.stdout
    else: outp=open(outf,'w')
    
    #Gather all easy predictions
    roots={}
    topcombos={}
    for line in open(summf,'r').readlines():
        if line[0]=='#':continue
        items=line.split()
        root=items[0][0:6]
        var=items[0]
        type=items[1]
        coverage=float(items[3])
        topcombo=items[4] #combo ID for top model according to density
        if type=='easy':
            if root not in roots.keys(): roots[root]=[var,]
            else: roots[root].append(var)
            topcombos[var]=topcombo
    #Bye(roots)

    #do TM-score on the consitutive regions
    for root in roots.keys():
        isos=roots[root]
        ref=root+'00'
        if ref not in isos: continue
        reff=os.path.join(rootd_fasta,ref)
        tpc=int(topcombos[ref]) #combo ID for top model according to density
        refcombo=os.path.join(rootd_out,ref[0:6],ref[6:8],'combo%02d.pdb.rebuilt'%(tpc,))
        if not os.path.exists(refcombo): continue
        #Bye(reff+' '+comboref)
        if len(isos)>1:
            isos.remove(ref) #extract reference sequence, only isoforms remain
            for iso in isos:
                isof=os.path.join(rootd_fasta,iso)
                cmd=alignx+' '+reff+' '+isof+' > '+algf #;Bye(cmd)
                os.system(cmd) #align to find constitutive region
                #Bye(algf)
                tpc=int(topcombos[iso]) #combo ID for top model according to density        
                isocombo=os.path.join(rootd_out,iso[0:6],iso[6:8],'combo%02d.pdb.rebuilt'%(tpc))
                if not os.path.exists(isocombo): continue
                cmd=tmScorex+' -a '+refcombo+' -b '+isocombo+' -i '+algf+' -j identical -d y'
                tm=os.popen(cmd).readline().strip()
                #tm=float(tm) #;Bye(tm)
                #outp.write('%s %s %4.2lf\n'%(ref,iso,tm))
                outp.write('%s %s %s\n'%(ref,iso,tm))
            
os.system('/bin/rm '+algf)
sys.exit(0)
