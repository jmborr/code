#for each variant set, create all possible pairwise alignments
#see /gpfs1/scratch/jose/spliceVariants/refsAboveC1/readme
import os,sys
from utilities.small_utilities import chomp
inputdir='/gpfs1/scratch/jose/spliceVariants/input'
list=chomp( open('variants.list','r').readlines() )
for i in range(0,len(list)-1):
    root=list[i][0:-1]
    si=list[i][-1]
    seqi=inputdir+'/'+list[i]
    for j in range(i+1,len(list)):
        sj=list[j][-1]
        seqj=inputdir+'/'+list[j]
        filename=root+'_'+si+'_'+sj+'.alignment.dat'
        os.system('$HOME/bin/align '+seqi+' '+seqj+' 0 > '+filename)
sys.exit(0)
