#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp,junkName,Bye
from spicker.spickerYangResultsManager import spickOut,TMscore

#retrieve neccessary files to calculate predictor-properties (template file from prospector, chain.dat, seq.dat, combo files, rst.dat)

#root='/gpfs1/scratch/jose/G8/out'
root='/gpfs1/scratch/jose/spliceVariants/wrong'

def updateInputFiles_G8(header):
    os.system('/bin/rm *pdb *.dat 2>/dev/null') #remove previous files inside directory
    if os.system('tar jxf '+root+'/'+header+'/'+header+'_tasserinp.tbz2 -C .'):return None
    if os.system('tar jxf '+root+'/'+header+'/'+header+'_out.tasser.tbz2 -C .'):return None
    if os.system('/bin/cp '+root+'/'+header+'/'+header+'rap3orienrev4s.pdb .'):return None
    return header+'rap3orienrev4s.pdb'

#retrieve neccessary files to calculate predictor-properties (template file from prospector, chain.dat, seq.dat, combo files, rst.dat)

def updateInputFiles_spliceVariants_wrong(header):
    os.system('/bin/rm *pdb *.dat 2>/dev/null') #remove previous files inside directory
    if os.system('/bin/cp '+root+'/prospector/'+header[1]+'/'+header+'/pdbbpdborienrev/'+header+'rap3orienrev5s.pdb .') or os.system('/bin/cp '+root+'/tasser/'+header[1]+'/'+header+'/seq.dat    .') or  os.system('/bin/cp '+root+'/tasser/'+header[1]+'/'+header+'/chain.dat  .') or os.system('/bin/cp '+root+'/tasser/'+header[1]+'/'+header+'/rst.dat    .') or os.system('/bin/cp '+root+'/tasser/'+header[1]+'/'+header+'/combo*.pdb .') or os.system('/bin/cp '+root+'/tasser/'+header[1]+'/'+header+'/closc*.pdb .'):
        return None #some file is missing
    return header+'rap3orienrev5s.pdb'

#return secondary structure content
def sec_str_content():
    content={'H':0.0,'E':0.0,'C':0.0} #number of residues assigned helix(H), strand(E) or coil(C)
    for line in open('seq.dat','r').readlines():
        state=line.split()[2]
        if state=='1': content['C']+=1.0
        elif state=='2': content['H']+=1.0
        else: content['E']+=1.0
    total=content['C']+content['H']+content['E']
    for key in content: content[key]=content[key]/total
    return content

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
pf.write('#Results are for most dense cluster only\n')
pf.write('#(1)header (2)easy/medm/hard (3)L (4)cluster-number-of-top-cluster (5)Z (6)closest template match by TM  (7)TMscore(closest-template-match,top-cluster) (8)<rmsd> (9)number-density (10)density (11)H (12)E (13)C (14)templ with highest Z (15)Z of templ with highest Z\n')
pf.write('#   1      2   3   4    5      6     7      8    9    10     11   12   13      14    15\n')
 
junkd=junkName()  ;  os.system('/bin/mkdir '+junkd)  ;  os.chdir(junkd) #working directory
for header in chomp(open(list,'r').readlines()):
    #put all files of interest in the working directory
#    templatefile=updateInputFiles_G8(header)
    templatefile=updateInputFiles_spliceVariants_wrong(header)
    if not templatefile: continue
    p=spickOut(dir='.') #instantiate a spickOut object
    if p.readError: #some error reading files, go to next header
        #pf.write(header+' ERROR: '+p.readError+'\n')
        continue
    print header
    data={}
    data['L']=p.length
    data['Z']=p.calcZscoreOfCluster(templatefile) #Z score
    data['template header']=p.templateOfCluster[p.densest]['libPDB']

    junk=junkName()
    open(junk,'w').writelines(p.templateOfCluster[p.densest]['templ'])
    data['TMscore(closest-template-match,top-cluster)']=float( TMscore(junk,p.combo[p.densest]) )
    os.system('/bin/rm '+junk)
    
    data['<rmsd>']=p.Rc_in[p.densest]
    data['number density']=p.numberdens[p.densest]
    data['density']=p.dens[p.densest]
    data['sec str']=sec_str_content() #dictionary with secondary stucture content
    data['easy medm hard']=chomp( open('chain.dat','r').readline() ).split()[1]
    data['templ with highest Z']=p.templates[0]['libPDB']
    data['Z of templ with highest Z']=p.templates[0]['Z']
    #output line
    pf.write('%8s %4s %3d %2d %5.2lf %8s %4.2lf %5.2lf %4.2lf %7.5lf %4.2lf %4.2lf %4.2lf %8s %5.2lf\n'%(header,data['easy medm hard'],data['L'],p.densest,data['Z'],data['template header'],data['TMscore(closest-template-match,top-cluster)'],data['<rmsd>'],data['number density'],data['density'],data['sec str']['H'],data['sec str']['E'],data['sec str']['C'],data['templ with highest Z'],data['Z of templ with highest Z']))
    
pf.close()

os.chdir(currd)  ;  os.system('/bin/rm -rf '+junkd) #clean-up
sys.exit(1)
