#!/usr/bin/python
"""Output several quantities related to secondary structure

output sequence length, fractional content of helixes and strands,
as well as number of helixes and strands
"""

import os,sys,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye,chomp,junkName
from secondary.ssManager import ssec

ih=inpHand('Usage: sec_analysis.py required args [optional args]\n'+locals()['__doc__'],
           ' -a _RA_seqdat seq.dat file of the original secondary assignments',
           ).parse(locals(),sys.argv)

x=ssec(seqdat=seqdat) #create secondary structure object from seqdat
if x.error: Bye(x.error)
x.initFractionalContent() #fractional content of every type
x.initChunks() #number of chunks of a given secondary type and minimal size

print ' L   %H   %E %(H+E) nH nE'
print '%3d %4.2f %4.2f %4.2f %2d %2d' % ( x.L,x.fracontent['helix'],x.fracontent['strand'],
                                    x.fracontent['helix']+x.fracontent['strand'],
                                   len(x.chunks['helix']),len(x.chunks['strand']) )
sys.exit(0)

