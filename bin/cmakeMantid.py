#!/usr/bin/env python

import os
import argparse
import sys
from pdb import set_trace as tr

parser = argparse.ArgumentParser(description='Configure and make Mantid')
parser.add_argument('--builddir', default='debug', help='"debug", "release", or other. Default is "debug"')
parser.add_argument('--cleanbuild', default=False, help='Remove all contents in target build directory prior to build? Default: no')
parser.add_argument('--ncore', default=7, help='Default: 7')
parser.add_argument('--target', default='all', help='either of "all", "doc", "doctest". Default: "all"')
parser.add_argument('--doctest', default=None, help='doctest a particular algorithm or fitfuncion')
parser.add_argument('--popmsg', default='True', help='Popup a finish job at the end? Default=True')
args=parser.parse_args()

rootd=os.environ['HOME']+'/repositories/mantidproject/build'

script='''#!/bin/bash

mkdir -p _ROOTD_/__BUILDDIR__
__DOCLEANBUILD__
cd _ROOTD_/__BUILDDIR__

_MAKETARGET_

__POPMSG__

############ B U I L D I N G   S C R I P T    F I N I S H E D ###########
'''

script=script.replace('_ROOTD_', rootd)
script=script.replace('__BUILDDIR__', args.builddir)
buildtype='Debug'
if args.builddir=='release': buildtype='Release'
script=script.replace('__BUILDTYPE__', buildtype)
script=script.replace('_TARGET_', args.target)

if str(args.cleanbuild) in ('yes', 'y',  'Y', 'Yes', 'YES', 'True', 'true', '1'):
    xxx='/bin/rm -rf {0}/{1}/*'.format(rootd,args.builddir)
    script=script.replace('__DOCLEANBUILD__',xxx)
else:
    script=script.replace('__DOCLEANBUILD__','')

if args.target=='all':
    xxx='''cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=__BUILDTYPE__ -DCMAKE_ECLIPSE_GENERATE_SOURCE_PROJECT=TRUE -DCMAKE_ECLIPSE_MAKE_ARGUMENTS=-j__NCORE__ -DCMAKE_ECLIPSE_VERSION=4.3  -DUSE_PRECOMPILED_HEADERS=TRUE -DDOCS_HTML=TRUE ../../mantid/Code/Mantid

echo " Building using __NCORE__ cores..."
sleep 3s
make -j__NCORE__ all'''
    script=script.replace('_MAKETARGET_', xxx)

elif args.target=='doc':
    xxx='''echo " Building docs-qthelp..."
sleep 9s
make -j__NCORE__ docs-qthelp

echo -e " Building docs-html..."
sleep 9s
make -j__NCORE__ docs-html'''
    script=script.replace('_MAKETARGET_',xxx)

elif args.target=='doctest':
    if not args.doctest:
        xxx='''echo -e "\n Building docs-test, *** NOTE THIS --> Storing results in file make_docs-test.log...\n"
sleep 9s
make -j__NCORE__ docs-test &> make_docs-test.log
'''
    else:
        xxx='python docs/runsphinx_doctest.py -m {0}/{1}/bin -R {2}'.format(rootd,args.builddir,args.doctest)
    script=script.replace('_MAKETARGET_',xxx)
else:
    print 'target not recognized. Nothing to do'
    sys.exit(1)

script=script.replace('__NCORE__', str(args.ncore))


popmsg='popwindow.py 0s "finished building _TARGET_ in Mantid"'
if args.popmsg.lower() in ('false', 'no', 'n', '0'):
    popmsg=''
script=script.replace('__POPMSG__', popmsg)

############ S T A R T I N G    B U I L D I N G   S C R I P T  ###########
print script
#os.system(script)
