#!/usr/bin/python
import sys,os,re
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import addAbsPath,inpHand

'''
Calculate average sequence identities from the howSimilarHongPredToPsiPred.dat files

'''

rootd='/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/hongpred'
gs='a0651  a1172  a1213  c0130  c1037  c1154'
comparisons=['hongpred vs  psipred=',' psipred vs original=','hongpred vs original=']
for g in gs.split():
    compf=os.path.join(rootd,g,'howSimilarHongPredToPsiPred.dat')
    for comparison in comparisons:
        avs={'A':0.0, 'H':0.0, 'B':0.0, 'L':0.0}
        lines=os.popen('grep -e "'+comparison+'" '+compf).readlines()
        for line in lines:
            avs['A']+=float(line[21:25])
            avs['H']+=float(line[26:30])
            avs['B']+=float(line[31:35])
            avs['L']+=float(line[36:40])
        N=len(lines)
        sys.stdout.write(g+' '+comparison)
        for x in ('A','H','B','L'):
            sys.stdout.write( ' %4.2lf'%(avs[x]/N) )
        sys.stdout.write('\n')

sys.exit(0)
