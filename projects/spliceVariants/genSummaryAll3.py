#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp,junkName,Bye
from spicker.spickerYangResultsManager import spickOut,TMscore

"""Predict top cluster (with different ranking methods). Then compare
to the \"native\" as given by the best annotated PDB chunk """

#global variables

rex={'tm':re.compile('TM-score    =\s+(\d\.\d+)'),
     'rmsd':re.compile('RMSD of  the common residues=\s+(\d+\.\d+)')
     }

def updateInputFiles(header,root,root_prosp=''):
    """retrieve neccessary files to calculate predictor-properties

    (template file from prospector, chain.dat, seq.dat, combo files,
    rst.dat)"""

    if not root_prosp:root_prosp=root
    os.system('/bin/rm *pdb *.dat 2>/dev/null') #remove previous files inside directory
    inpdir=root+'/'+header[:6]+'/'+header[6:]
    inpdir_prosp=root_prosp+'/'+header[:6]+'/'+header[6:]
    if os.system('tar xf '+inpdir_prosp+'/'+header+'.out.prospector.tar pdbbpdborienrev/'+header+'rap3orienrev5s.pdb -C .') or os.system( 'tar xf '+inpdir+'/'+header+'.in.tasser.tar chain.dat  -C .') or os.system('tar xf '+inpdir+'/'+header+'.out.spicker.tar -C .'):
        return None #some file is missing
    type=open('./chain.dat','r').readline().split()[1] #easy,medm,hard
    os.system('/bin/mv ./pdbbpdborienrev/'+header+'rap3orienrev5s.pdb .')
    return header+'rap3orienrev5s.pdb',type,inpdir

#Execution begins here
ih=inpHand('Usage: genSummaryAll.py [options]\nGiven a list of headers, and a xxxxx-globbed directory,\ngoes to each de-globbed directory and retrieves all neccessary info to generate the set of properties that Ying\'s TM score stimator needs\n',
           ' -a _RA_list list of headers',
           ' -b _A_outf output filename  (def=./summary.all)',
           ' -c _RA_root_scratch root directory for tasser and spicker',
           ' -d _A_root_prosp root directory for prospector output (def: value of root variable)',
           )
ih.parse(locals(),sys.argv)

#defaults
currd=os.getcwd()
if not outf: outf=currd+'/summary.all'
if not root_prosp: root_prosp=root_scratch

pf=open(outf,'w')
pf.write('#(1)header (2)easy/medm/hard (3)L (4)coverage-of-best-PDB-template (5)topID-by-Ying-rank (6)TM-score to best_PDB_templ (7)rmsd to best_PDB_templ  topID-by-Density-rank (6)topID-by-Cscore-rank\n')
pf.write('#   1      2   3   4    5  6\n')
 
wd=junkName()  ;  os.system('/bin/mkdir '+wd)  ;  os.chdir(wd) #working directory
for header in chomp(open(list,'r').readlines()):
    #put all files of interest in the working directory
    templatefile,type,inpdir=updateInputFiles(header,root_scratch,root_prosp=root_prosp)
    if not templatefile: continue
    p=spickOut(dir='.') #instantiate a spickOut object
    if p.readError: #some error reading files, go to next header
        #pf.write(header+' ERROR: '+p.readError+'\n')
        continue
    print header
    data={}
    data['L']=p.length
    #get the top cluster ranked through three different approaches
    topIDs={}
    topIDs['Ying']=p.rankIDsByYingPredictedTM(templatef=templatefile,type=type)[0]
    topIDs['Density']=p.densest
    topIDs['Cscore']=p.rankIDsByCscores(templatef=templatefile)[0]
    best_pdb_templ=os.path.join(inpdir,'best_PDB_templ.pdb')
    #coverage of best_pdb_templ
    if os.path.exists(best_pdb_templ):
        cov=float(os.popen('wc '+ best_pdb_templ).readline().split()[0])/data['L']
    else:
        cov=0.0
    buf='%8s %4s %3d %4.2lf'%(header,type,data['L'],cov)
    #copy related combo and make appropriate links
    for method in ('Ying','Density','Cscore'):
    #for method in ('Density',):
        combo='combo'+`topIDs[method]`+'.pdb'
        os.system('/bin/cp '+combo+' '+inpdir)
        link=inpdir+'/top_by_'+method+'.pdb'
        os.system('/bin/rm '+link+' 2>/dev/null && /bin/ln -s '+inpdir+'/'+combo+' '+link)
        if os.path.exists(best_pdb_templ):
            tmoutput=''.join(os.popen('TMscore '+combo+' '+best_pdb_templ).readlines())
            tm=float(rex['tm'].search(tmoutput).group(1))
            rmsd=float(rex['rmsd'].search(tmoutput).group(1))
        else:
            tm=0.0
            rmsd=99.9
        buf+=' %2d %4.2lf %4.1lf'%(topIDs[method],tm,rmsd)
        
    pf.write(buf+'\n')

pf.close()

os.chdir(currd)  ;  os.system('/bin/rm -rf '+wd) #clean-up
sys.exit(1)
