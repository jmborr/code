#!/usr/bin/python

import sys,os
from spicker.spickerTMresultsManager import spickOut
from utilities.small_utilities import chomp

#listf=sys.argv[1] #first argument
#root=sys.argv[2]
listf='/gpfs1/active/casp7/output/jose/assess/pdbs.list'
root='/gpfs1/active/casp7/output/jose/assess'

print '#header L   densest  bestTM  avtm'
for header in chomp(open(listf,'r').readlines()):
    T=chomp(os.popen('grep -e "'+header+' " '+root+'/'+'targetsNranges.dat |tr -s \' \'| cut -d\' \' -f 2').readline())
    ranks={}
    dir=root+'/'+T+'/'+header
    str=dir+'/'+header+'.pdb'  #;print 'dir=',dir,'str=',str;sys.exit(1)
    s=spickOut(dir=dir)
    if s.readError:  print s.readError
    else:
        ranked=s.rankIDsByDens()  #cluster ID's ranked by decreasing density
        for i in range(0,5): ranks[ranked[i]]=1+i #returns rank for a given cluster ID
        idd=ranked[0]             #cluster ID of densest cluster
        s.filltmAls(str)
        avtm=0 ; n=0
        for i in range(0,4):
            for j in range(i+1,5):
                avtm+=float( s.TMidId(ranked[i],ranked[j]) )
                n+=1
        avtm=avtm/n
        idt=s.rankIDsBytmAlstoStrt(list=ranked[0:5])[0] #id of best cluster by TMalign
        rt=ranks[idt] #rank of cluster with best TM score
        header=header+' '*(6-len(header)) #back-padd with blank spaces
        line='%5s %3d  1 %5.3lf %2d %5.3lf %5.3lf'%(header,s.length,s.tmAls[idd],rt,s.tmAls[idt],avtm)
        print line
