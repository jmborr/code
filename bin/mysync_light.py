#!/usr/bin/python

import os
import argparse
from pdb import set_trace as tr

include_extensions=['cif', 'gif', 'in', 'jpeg', 'ipynb', 'docx', 'pdb', 'png', 'pse', 'tleap', 'conf']

parser = argparse.ArgumentParser(description='sync files of extensions %s'%', '.join(include_extensions))
parser.add_argument('origin',help='source directory')
parser.add_argument('target', help='target directory')
parser.add_argument('--maxsize', default='100MB', help='do not sync files larger than max-size. Default: 100MB')
parser.add_argument('--optional_extensions', help='pass extra extensions to sync, separated by a white space')
args=parser.parse_args()

if args.optional_extensions:
    include_extensions+=args.optional_extensions.split()

mssg="""
********************************************************
*                                                      *
*            SYNCING DESKTOP FILES                     *
*                                                      *
*  DO NOT USE THIS CONSOLE WHILE OPERATION IN PROGRESS *
*                                                      *
********************************************************
"""
os.system('clear')
print mssg
os.system('sleep 5s')

for extension in include_extensions:
    cmd  = 'cd "%s"; '%args.origin
    cmd += 'find . -name "*.%s" -print0 | '%extension

    rsync_string= ' rsync -aur --delete --files-from=- --from0'
    if args.maxsize: rsync_string+=' --max-size=%s'%args.maxsize
    rsync_string += ' ./ "%s"'%args.target

    cmd += rsync_string
    print cmd
    os.system( cmd )

#delete files in target not existing in origin. The --delete option in rsync does not work
buffer='#!/bin/bash\n\n'
for extension in include_extensions:
    cmd='''cd "%s"
IFS=$(echo -en "\\n\\b")   #tell find command that separator is newline, not whitespace
for file in `find . -name "*.%s"`;do
  if [ ! -f "%s/$file" ];then
    /bin/rm "$file"
  fi
done'''%(args.target, extension, args.origin)
    buffer += cmd + '\n\n'
open('/tmp/junk_mysync_light.sh','w').write(buffer)
os.system('bash /tmp/junk_mysync_light.sh; /bin/rm /tmp/junk_mysync_light.sh')


mssg="""
****************************
*                          *
* FINISHED SYNCING...BYE ! *
*                          *
****************************
"""

print mssg
os.system('popwindow.py "Syncing finished. Please unmount external drive"')

