from utilities.small_utilities import chomp
#read DISFIT creating a two-key dictionary disfit2[aa1][aa2]
disfit={}
disfit2={}
p=open('/gpfs1/active/jose/code/f77/tertiary/DISTFIT','r')
#list containing three-letter codes for amino acids
#     GLY   ALA   SER   CYS   VAL   THR   ILE   PRO...
aas=chomp(p.readline()).split()
N=len(aas) #; print N #number of codes
#read remaining lines. They contain one three-letter codes and then distances
#GLY  450   450   494   482   531   537   576...
for line in chomp(p.readlines()):
    list=line.split()
    aa1=list[0]
    vals=list[1:] # contact cut-off distances
    disfit[aa1]={}
    disfit2[aa1]={}    
    for i in range(0,N):
        aa2=aas[i]
        x=float(vals[i])/100.0 #remember DISFIT distance units are 0.01Angstroms
        disfit[aa1][aa2]=x
        disfit2[aa1][aa2]=x*x #we store the square of the distances
p.close()
