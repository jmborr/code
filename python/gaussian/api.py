#!/bin/usr/python

from g98 import gaussian
from aCLIMAX.aCLIMAXv5 import atom,freq,aCLIMAX

def g98freq_to_aCLIMAX(inf,outf):
    """translate g98 frequency ouput job to aCLIMAX
    inf: g98 output for frequency job
    outf: aCLIMAX input frequencies"""
    h_g=gaussian(inf)
    
    ac.outputFreqFile(outf)
    
