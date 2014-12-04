#!/usr/bin/python

import os
import sys

root='/media/samsung'
pairs = {'/projects/development/' : root+'/projects/development/',
         '/projects/research/'    : root+'/projects/research/',
         '/projects2/development/': root+'/projects/development/',
         '/projects2/research/'   : root+'/projects/research/',
         '/home/jbq/'             : root+'/home/jbq/'
         }

if not os.path.exists(root):
    os.system('popwindow.py "Please mount samsung external hard drive in order to proceed with the syncing"')
    sys.exit(1)

mssg="""
********************************************************
*                                                      *
*            SYNCING DESKTOP FILES                     *
*                                                      *
*  DO NOT USE THIS CONSOLE WHILE OPERATION IN PROGRESS *
*                                                      *
********************************************************
"""
print mssg

for (origin,destination) in pairs.items():
    cmd = 'rsync -aur --exclude "*~" --exclude ".*~" --exclude "#*" '+\
        '--exclude ".#*" --exclude "junk*" --exclude "*core" '+\
        origin+' '+destination
    #print cmd
    os.system( cmd )

mssg="""
****************************
*                          *
* FINISHED SYNCING...BYE ! *
*                          *
****************************
"""
print mssg
os.system('popwindow.py "Syncing finished. Please unmount external drive"')
