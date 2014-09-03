#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp,junkName,Bye
from spicker.spickerYangResultsManager import spickOut,TMscore

"""

"""
root='/gpfs1/scratch/jose/spliceVariants/out2'

def updateInputFiles(header):
    """retrieve neccessary files to calculate predictor-properties

    (template file from prospector, chain.dat, seq.dat, combo files,
    rst.dat)"""

    os.system('/bin/rm *pdb *.dat 2>/dev/null') #remove previous files inside directory
    inpdir=root+'/'+header[:6]+'/'+header[6:]
    if os.system('tar xf '+inpdir+'/'+header+'.out.prospector.tar pdbbpdborienrev/'+header+'rap3orienrev5s.pdb -C .') or os.system( 'tar xf '+inpdir+'/'+header+'.in.tasser.tar chain.dat  -C .') or os.system('tar xf '+inpdir+'/'+header+'.out.spicker.tar -C .'):
        return None #some file is missing
    type=open('./chain.dat','r').readline().split()[1] #easy,medm,hard
    os.system('/bin/mv ./pdbbpdborienrev/'+header+'rap3orienrev5s.pdb .')
    return header+'rap3orienrev5s.pdb',type,inpdir

#Execution begins here
ih=inpHand('Usage: genSummaryAll.py [options]\nGiven a list of headers, and a xxxxx-globbed directory,\ngoes to each de-globbed directory and retrieves all neccessary info to generate the set of properties that Ying\'s TM score stimator needs\n',
           ' -a _RA_list list of headers',
           ' -b _A_outf output filename  (def=./summary.all)',
           )
ih.parse(locals(),sys.argv)

#defaults
currd=os.getcwd()
if not outf: outf=currd+'/summary.all'

pf=open(outf,'w')
pf.write('#(1)header (2)easy/medm/hard (3)L (4)topID-by-Density-rank (5)topID-by-Cscore-rank (6)topID-by-Ying-rank\n')
pf.write('#   1      2   3   4  5  6\n')
 
wd=junkName()  ;  os.system('/bin/mkdir '+wd)  ;  os.chdir(wd) #working directory
for header in chomp(open(list,'r').readlines()):
    #put all files of interest in the working directory
    templatefile,type,inpdir=updateInputFiles(header)
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
    #copy related combo and make appropriate links
    for method in ('Ying','Density','Cscore'):
        combo='combo'+`topIDs[method]`+'.pdb'
        os.system('/bin/cp '+combo+' '+inpdir)
        link=inpdir+'/top_by_'+method+'.pdb'
        os.system('/bin/rm '+link+' 2>/dev/null && /bin/ln -s '+inpdir+'/'+combo+' '+link)

    #output line
    pf.write('%8s %4s %3d %2d %2d %2d\n'%
             (header,type,data['L'],topIDs['Density'],topIDs['Cscore'],topIDs['Ying'])             
            )
    #pf.write( '%8s %4s %3d %2d\n'%(header,type,data['L'],topIDs['Density']) )

pf.close()

os.chdir(currd)  ;  os.system('/bin/rm -rf '+wd) #clean-up
sys.exit(1)
