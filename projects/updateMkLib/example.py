#!/usr/bin/python

from align_yang_f77 import align

seq1='AGVKLIKKKGDWLVY'
seq2='WKYSTAGLKLIILIGGSWTYNQ'

(nid,nali)=align(seq1,seq2)
print nid
