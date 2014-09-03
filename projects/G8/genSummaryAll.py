#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand,deglobb
from utilities.small_utilities import chomp,junkName,Bye
from spicker.spickerYangResultsManager import spickOut,TMscore
from tasser.chainDatReader import chainDatReader

#retrieve neccessary files to calculate predictor-properties (template file from prospector, chain.dat, seq.dat, combo files, rst.dat)

root='/gpfs1/scratch/jose/G8/out'                #for updateInputFiles_G8
#root='/gpfs1/scratch/jose/spliceVariants/wrong' #for updateInputFiles_spliceVariants_wrong
#root='/gpfs1/active/shashi'                     #for updateInputFiles_shashi_benchmark

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

def updateInputFiles_shashi_benchmark(header):
    dir=root+'/Tasserinputs/'+header
    if os.system('/bin/cp '+dir+'/rmsinp .'):return None

    #check sequence length above 100 residues
    if int(chomp(open('rmsinp','r').readline()).split()[1]) <= 100:return None
    
    if os.system('/bin/cp '+dir+'/chain.dat .') or os.system('/bin/cp '+dir+'/seq.dat .'):return None    
    os.system('/bin/rm *pdb *.dat CA 2>/dev/null') #remove previous files inside directory
    dir=root+'/Tasser_3035/'+header
    if os.system('tar jxf '+dir+'/'+header+'_spA.bz -C .'):return None
    if os.system('/bin/cp '+dir+'/CA .'):return None
    dir=root+'/NEWWORK/Prospector_30_35/pdbfpdborienrev'
    if os.system('/bin/cp '+dir+'/'+header+'rap3orienrev4s.pdb .'):return None
    return header+'rap3orienrev4s.pdb'
    
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
pf.write('#(1)header (2)easy/medm/hard (3)L (4)cluster-number-of-top-cluster (5)Z (6)closest template match by TM  (7)TMscore(closest-template-match,top-cluster) (8)<rmsd> (9)number-density (10)density (11)H (12)E (13)C (14)templ with highest Z (15)Z of templ with highest Z (16)rank of top cluster in the list of templates in chain.dat (17)template header in chain.dat for which has biggest TMscore to top cluster (18)Max( TM(chain.dat template,top cluster) ) (19)Max( TM(chain.dat template,top cluster) )- TM(first chain.dat template,top cluster)\n')
pf.write('#   1      2   3   4    5      6     7      8    9    10     11   12   13      14    15   16     17    18   19\n')

junkd=junkName()  ;  os.system('/bin/mkdir '+junkd)  ;  os.chdir(junkd) #working directory
for header in chomp(open(list,'r').readlines()):
    #put all files of interest in the working directory
    templatefile=updateInputFiles_G8(header)
#    templatefile=updateInputFiles_spliceVariants_wrong(header)
#    templatefile=updateInputFiles_shashi_benchmark
    if not templatefile: continue

    p=spickOut(dir='.') #instantiate a spickOut object
    if p.readError: #some error reading files, go to next header
        #pf.write(header+' ERROR: '+p.readError+'\n')
        continue
    print header

    #calculate the rank of the densest cluster with respect to the list of templates in chain.dat
    ch=chainDatReader(natf=p.combo[p.densest])  #read the templates of chain.dat
    tmbest=0 ; tm1=0 ;  n=0
    while ch.readTempl():
        n+=1  ;  tm=ch.TMscore()
        if n==1: tm1=tm
        if tm>tmbest:
            tmbest=tm  ;  rank=n  ;  templheader=ch.template_header
            
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

    data['rank_chain']=rank
    data['templheader']=templheader
    data['tmbest']=tmbest
    data['dTM']=tmbest-tm1 #diference in TM of the best template
    
    #output line
    pf.write('%8s %4s %3d %2d %5.2lf %8s %4.2lf %5.2lf %4.2lf %7.5lf %4.2lf %4.2lf %4.2lf %8s %5.2lf %2d %8s %4.2lf %4.2lf\n'%(header,data['easy medm hard'],data['L'],p.densest,data['Z'],data['template header'],data['TMscore(closest-template-match,top-cluster)'],data['<rmsd>'],data['number density'],data['density'],data['sec str']['H'],data['sec str']['E'],data['sec str']['C'],data['templ with highest Z'],data['Z of templ with highest Z'],data['rank_chain'],data['templheader'],data['tmbest'],data['dTM']))
    
pf.close()

os.chdir(currd)  ;  os.system('/bin/rm -rf '+junkd) #clean-up
sys.exit(1)
