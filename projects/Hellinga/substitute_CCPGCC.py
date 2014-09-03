#!/usr/bin/python

"""Use CCPGCC sequence to substitute several sequences in designated
positions (turns neat at one terminus or in the interior of the
sequence). The selected position will not substitute the three HIS
that the sequences contain"""

import os,sys
from utilities.codedir import scratchdir
from utilities.small_utilities import chomp,Bye
from seq.fastaManager import importFastaEntry,writeFasta
from seq.letters import one2three

#global variables
targets=('a0216','a0512','a0653','a0857','a0984')
XXX=[one2three[X] for X in 'CCPGCC'] #;Bye(XXX)
subst={'close_to_terminus':{'a0216':16,'a0512':145,'a0653':29,'a0857':145,'a0984':145},
       'internal_turn':{'a0216':108,'a0512':41,'a0653':55,'a0857':50,'a0984':95}
       } #marks where in the sequence to begin substitut. First res has index "1", not "0"
root=os.path.join(scratchdir,'Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/Cys_rich_turn')

#main program
for dir in ('close_to_terminus','internal_turn'):
    for target in targets:
        listf=os.path.join(root,dir,target,'seqsAbove0.4.list') #;Bye(listf)#list of headers
        print dir,target
        for header in chomp(open(listf,'r').readlines()):
            x=subst[dir][target]-1 #marks position where substitute CCPGCC
            #create new fasta file            
            fastaf=os.path.join(root,dir,target,header,header) #;Bye(fastaf)
            h,seq=importFastaEntry( open(fastaf,'r') )
            seq=seq[:x]+'CCPGCC'+seq[x+7:] #;Bye(seq)
            writeFasta(header,seq,filename=fastaf,n=50) #rewrite fasta file with substituted seq.
            #create new seq.dat file
            seqdatf=os.path.join(root,dir,target,header,'seq.dat') #;Bye(seqdatf)
            ls=open(seqdatf,'r').readlines()
            for i in range(6): ls[x+i]=ls[x+i][:6]+XXX[i]+ls[x+i][9:]
            open(seqdatf,'w').write( ''.join(ls) )
                        
sys.exit(0)
