#!/usr/bin/python

'''
'''

import os,sys
from seq.letters import valid_A
from utilities.small_utilities import chomp,junkName,Bye
from utilities.codedir import codedir
from jobs.job import pastry

#input template from prospector, output chain.dat directly
def formatTempl(prospTempl,chainTempl='chain.dat',type='easy'):
    ppt=open(prospTempl,'r')
    cpt=open(chainTempl,'w')
    #number of templates equals number of TER lines
    ntempl=int(os.popen('grep TER '+prospTempl+'|wc -l').readline().strip())
    cpt.write('%6d %s\n'%(ntempl,type))
    n=0 #template number currently being read
    while n<ntempl:
        #read one template
        ncov=0 #number of residues covered by the template
        buf='' #buffer to hold lines
        while 1:
            line=ppt.readline()
            if line.find('TER') >= 0 or line=='':
                break #stop if end of template or end of file
            if line.find('ATOM')>=0 and line.find(' CA ')>=0:
                ncov+=1
                l=int(line[22:26])
                x=float(line[30:38]) ; y=float(line[38:46]) ; z=float(line[46:54])
                buf+='%6d%11.3lf%11.3lf%11.3lf\n'%(l,x,y,z)
        cpt.write('%6d\n'%(ncov)+buf) #write number of covered residues and coordinates
        n+=1
    ppt.close()
    cpt.close()
