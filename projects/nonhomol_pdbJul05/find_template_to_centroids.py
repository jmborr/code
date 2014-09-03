#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from random import randint
from seq.get_seq import get_seq
from utilities.codedir import yangpdb,projectsdir
from utilities.small_utilities import chomp
from seq.align_yang_f77 import seq_id
from tm_align.tmalign_Yang_f77 import tmalign_yang
from jobs.job import pastry

def read_snapshot(ipf,tmpf):
    opf=open(tmpf,'w')
    p=re.compile('TER|END')
    while 1:
        line = ipf.readline();
        if not line: return None #return None if end of file reached
        opf.write(line)
        if p.match(line): #we found TER or END, close output file and return
            opf.close()
            return 1

## function attractors
#  centroids: centroids.pdb file
#  0<seqidcutoff<=1
#  attf: attractors.pdb file, store the alignment of each centroid to
#  its attractor (template maximizing the TM-score
def attractors(centroids,seqidcutoff=0.30,attf='',templDir=yangpdb,
               templList=projectsdir+'/db/yangpdb_list'):
    tmpf=os.getcwd()+'/'+`randint(0,1000)`+'.pdb'      #;print tmpf
    if attf:
        tmpf2=os.getcwd()+'/'+`randint(0,1000)`+'.txt'
        os.system('/bin/rm '+attf+' 2>/dev/null') #erase old attractors.pdb file, if existing
    ipf=open(centroids,'r')
    #read first centroid
    read_snapshot(ipf,tmpf)
    #get sequence in one letter code
    targetseq=get_seq(pdbfile=tmpf)      #;print 'targetseq=',targetseq
    #read all template headers, replace '\n' with '.pdb'
    templs=chomp(open(templList,'r').readlines(),replace='.pdb')
    #rewind centroids file and do loop search for all centroids
    ipf.seek(0)
    while 1:
        if not read_snapshot(ipf,tmpf): break #store centroid pdb in tmpf
        tmscoreMax=0.0 #initialize best tm-score between centroid and templates
        besttempls=[]  #initialize list of best templates
        #cicle over all templates
        for templ in templs:
            templf=templDir+'/'+templ #;print 'templf='+templf #pdb file of template
            if not os.path.exists(templf): continue #non-existing file
            #get template sequence in one letter code
            templseq=get_seq(pdbfile=templf) #;print 'templseq=',templseq
            #if target and template have low sequence identity, then..
            if seq_id(targetseq,templseq)<seqidcutoff:
                tmscore=tmalign_yang(templf,tmpf) #important template before target
                if tmscore>tmscoreMax:                    
                    tmscoreMax=tmscore  ;print tmscoreMax
                    besttempl=templ

        #store best template found, remove '.pdb'
        besttempl=besttempl.replace('.pdb','')
        besttempls.append(besttempl)
        
        #find aligned portion of the best template to the target, if required
        if attf:
            tm='%5.3f'%(tmscoreMax)                       
            templf=templDir+'/'+besttempl+'.pdb'
            si='%4.2f'%(seq_id(targetseq,get_seq(pdbfile=templf)))
            pastry('echo "HEADER '+templf+' seqid='+si+'" >> '+attf)
            #store aligned portions in second temporary file
            tmalign_yang(templf,tmpf,alignf=tmpf2)
            #append aligned portions to the attractors file
            pastry('cat '+tmpf2+'>> '+attf)
            
    #remove temporary file(s)
    pastry('/bin/rm '+tmpf)  
    if tmpf2: pastry('/bin/rm '+tmpf2) 
    return besttempls
    
#execute as standalone script
if __name__=='__main__':
    ih=inpHand('Usage: find_template_to_centroids.py [options]',
               '  -a _AR_centroids centroids.pdb file',
               '  -b _A_outf output file name (def: ./attractors.pdb)',
               '  -c __seqidcutoff sequence identity cutoff (def:0.30)',
               '  -d _A_templDir directory where templates are (def:yangpdb)',
               '  -e _A_templList list of templates in templDir (def:projectsdir/db/yangpdb_list)')    
    ih.parse(locals(),sys.argv)
    
    if not seqidcutoff: seqidcutoff=0.30
    else: seqidcutoff=float(seqidcutoff)

    if not outf: outf=os.getcwd()+'/attractors.pdb'

    if not templDir: templDir=yangpdb
    if not templList:templList=projectsdir+'/db/yangpdb_list'

    #find best templates to each centroids, store the aligned portions in outf
    attractors(centroids,seqidcutoff=seqidcutoff,attf=outf,
               templDir=templDir,templList=templList)
        
