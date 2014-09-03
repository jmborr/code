#!/usr/bin/python

'''
File to read a TASSER trayectory (like rep1.tra) and create a set of structures, one per snapshot

'''

import os,sys,re
from utilities.small_utilities import chomp,junkName,Bye
from inputArgs.inputArgs import addAbsPath,inpHand
from tasser.repTraReader import repTraReader

ih=inpHand('Usage: genStrucFromDecoys.py [options]',
           ' -a _RA_traf trayectory file',
           ' -b _RA_natf native file name (default="./CA")',
           ' -c _RA_strd directory where structures will be stored. Will create if it does not exists',
           ' -d _R_prefix single letter prefix for structure names',
           )
ih.parse(locals(),sys.argv)

if not os.path.exists('strd'): os.system('/bin/mkdir -p '+strd) #may create directory
prefix=prefix[0] #keep only first letter
if not re.compile('[a-z,A-Z]').match(prefix): Bye('Pass a letter for prefix, please')
reader=repTraReader(traf) #initialize reader object
reader.readNextSnapShot() #for some reason, first and second snapshots are same in filtered trajectory
i=0
while reader.readNextSnapShot():
    outf=os.path.join( strd,prefix+'%05d'%(i) )
    lastshot=reader.lastshot
    lastshot.writepdb(outf,natf,extraHeader='REMARK '+lastshot.header+'\n') #write snapshot to file
    i+=1
