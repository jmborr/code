#!/usr/bin/python

'''
'''
import os,sys
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import deglobb

begin='<html><head><title>Implementation of ab initio TASSER on X representative PDB proteins below 200 residues</title></head>\n<body>\n<center>\n<h1><font size="5">\nImplementation of TASSER on X representative PDB proteins below 200 residues\n</font></h1>\n</center>\nIn this webpage, we present the result of protein structure predictions by ab initio TASSER (Threading ASSEmbly Refinement) algorithm on X representative PDB benchmark proteins.\nBecause of the huge size of the trajectory files (~80G), which includes coordinates of C-alpha and side-chain center of mass of snapshot structures of Monte Carlo simulations, they are currently not uploaded on this server. However, those files are available when requested via Email: <a href="mailto:jmborr@gatech.edu">jmborr@gatech.edu</a>.\n<br>\n<br>\n\n<h2>Summary of folding results:</h2>\n<ul>\n\t<li> ID --- PDB code of the proteins, with a link to the native CA coordinates</li>\n\t<li> Lenth</li>\n\t<li> Top Cluster, according to ranking of the clusters by sctrutural density after application of the <a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=pubmed&amp;dopt=Abstract&amp;list_uids=15011258">SPICKER\n</a> algorithm: TM score and RMSD (in Angstrom) with a link to the model CA coordinates</li>\n\t<li> Best Cluster out of the five top clusters: TM score, RMSD (in Angstrom) and rank. Includes a link to the model CA coordinates</li>\n\n<table border="1">\n<tbody>\n\t<tr>\n\t\t<td align="center">ID</td>\n\t\t<td align="center">length</td>\n\t\t<td align="center">Top: TM RMSD</td>\n\t\t<td align="center">Best: TM RMSD rank</td>\n\t</tr>\n\n'

end='\t</tbody>\n</table>\n</body></html>'

root='/gpfs1/scratch/jose/createInputsTasser/x/xxxxx/abinitio'
headers=[] ; infos={}
for l in open('/gpfs1/scratch/jose/nonhomol_pdbJul05/analysis/summaryMoreBelow200.all','r').readlines():
    if l[0]=='#': continue #the line is a comment
    header=l[0:5]
    headers.append(header)
    infos[header]= l[5:] #l[0:4] is the header, the rest is more info

print begin
table={}
for header in headers:
    dir=deglobb(root,header)
    if header in infos.keys(): info=infos[header].split()
    else: continue
    L=info[0] ; tmT=info[9] ; rmsdT=info[7] ; tmB=info[18] ; rmsdB=info[12] ; rank=info[17]
    table[header]='<tr>\n\t<td align="center"><a href="abinitio/CA/'+header+'.CA.bz2">'+header+'</a></td>\n\t\t<td align="center">'+L+'</td>\n\t\t<td align="center"><a href="abinitio/model/'+header+'.top.bz2">'+tmT+' - '+rmsdT+'</a></td>\n\t\t<td align="center"><a href="model/'+header+'.best.bz2">'+tmB+' - '+rmsdB+' - '+rank+'</a></td>\n\t</tr>\n\n'

headers=table.keys() ; headers.sort()
for header in headers: print table[header]
print end


