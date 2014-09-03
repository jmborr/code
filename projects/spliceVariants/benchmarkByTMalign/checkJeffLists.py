#!/usr/bin/python
'''
the file format of  101m_.hist file is as follows
         20           25                
101m_     1.000    1.000  154   0.000 #fist line is the target itself*/
1ut0A     0.299    0.910  151   1.381 #homolog seqIDbyTMalign/L(101m_) TM length-align rmsd-alignm
...

'''
import os,sys
from utilities.small_utilities import chomp,junkName,Bye
from TMalign.TMalignManageResults import TMalignOut
from seq.alignYangManageResults import gappedSeq

templd='/gpfs1/archive/adrian/dat/pdb_aug06.apr06' #template directory
pdbd=templd+'/CA'
seqd=templd+'/input'
histd='/gpfs1/active/skolnick/struhistaug06'
scd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
#listf=scd+'/toy.list'
listf=scd+'/list'
junk1=junkName() ; os.system('/bin/touch '+junk1) #temporary file
junk2=junkName() ; os.system('/bin/touch '+junk2)

outd='/gpfs1/scratch/jose/spliceVariants/benchmarkByTMalign'
outfdat=outd+'/alignmentsByTM.dat' ; outdat=open(outfdat,'w') #store alignments in sequence
outfidx=outd+'/alignmentsByTM.idx' ; outidx=open(outfidx,'w') #index file to outfdat
#for the moment we don't output superpositions, since it uses space 4GB in hard drive
#supfdat=outd+'/superposByTM.dat'   ; supdat=open(supfdat,'w') #store alignments in coordinates
#supfidx=outd+'/superposByTM.idx'   ; supidx=open(supfidx,'w') #index file to supfdat


#First gather length of chains
length={}
headers=chomp(open(listf,'r').readlines())
for header in headers:
    pdb1=pdbd+'/'+header+'.pdb'
    length[header]=int(os.popen('grep "ATOM " '+pdb1+' | grep " CA " | wc -l').readline().strip())

#Now get the alignments.
n=0
for header in headers:
    pdb1=header+'.pdb'; n+=1
    #We copy to local directory because this way output from TMalign
    #will show the header string as Chain 1    
    os.system('/bin/cp '+pdbd+'/'+pdb1+' .')
    l1=length[header]
    if l1>300: continue #keep only alignments with sequences below 300 residues
    print pdb1,n #will tell us of progress
    histf=histd+'/'+header+'.hist' #file with list of structural homologs
    pin=open(histf,'r') ; pin.readline() ; pin.readline()
    while(1):
        line=pin.readline()
        if not line: break #reached end of file
        list=line.split()
        homolog=list[0] #;print homolog
        if length[homolog]>l1: continue #keep only alignments for which header is the longest chain
        pdb2=homolog+'.pdb' #;print pdb2
        os.system('/bin/cp '+pdbd+'/'+pdb2+' .')
        os.system('TMalign '+pdb2+' '+pdb1+' -o '+junk1+' > '+junk2)
        os.system('/bin/rm '+pdb2)
        posalg=`outdat.tell()`
        info=homolog+' '+header+' - '+header+' '+homolog+'\n'
        outidx.write(posalg+':'+info)
        outdat.write('# '+info)
        outdat.writelines(open(junk2,'r').readlines()) #store the alignment of the sequence
        outdat.write('\n')
        #possup=`supdat.tell()`            
        #supidx.write(possup+':'+header+' '+homolog+'\n'+possup+':'+homolog+' '+header+'\n')
        #supdat.writelines('# '+header+' '+homolog+' - '+homolog+' '+header+'\n')
        #supdat.writelines(open(junk1,'r').readlines()) #store alignment of the coordinates
    os.system('/bin/rm '+pdb1)
            
os.system('/bin/rm '+junk1+' '+junk2) #cleaning of temporary files
outdat.close()
outidx.close()        
#supdat.close()
#supidx.close()
sys.exit(0)
