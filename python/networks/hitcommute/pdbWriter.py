"""
    Program that writes a custom PDB file for subsequent processing
"""

import os, string, sys

class PDBWriter:
    def __init__(self):
        return self
    
    def write(self, fname):
        f = open(fname, 'w')
        