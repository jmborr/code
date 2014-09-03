#!/usr/bin/python

import sys
from random import shuffle

l=open('/gpfs1/scratch/jose/template_free_loops/list','r').readlines()
shuffle(l)
print ''.join(l)
sys.exit(0)
