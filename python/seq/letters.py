#!/usr/bin/python

import re

#PDB average amino acid composition according to Robinson and Robinson, PNAS 88 8880-8884 (1991)
background_probs={'A':0.07805,'C':0.01924,'D':0.05364,'E':0.06295,'F':0.03856,'G':0.07377,'H':0.02199,'I':0.05142,'K':0.05744,'L':0.09020,'M':0.02243,'N':0.04487,'P':0.05202,'Q':0.04264,'R':0.05130,'S':0.07120,'T':0.05841,'V':0.06441,'W':0.01330,'Y':0.03216}

valid_A='ACDEFGHIKLMNPQRSTVWY'
valid_AAA=('ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLN', 'GLU', 'GLX',
           'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
           'THR', 'TRP', 'TYR', 'UNK', 'VAL')

three2one={'ALA':'A','VAL':'V','PHE':'F','ILE':'I','LEU':'L','PRO':'P',
           'MET':'M','ASP':'D','GLU':'E','LYS':'K','ARG':'R','SER':'S',
           'THR':'T','TYR':'Y','HIS':'H','CYS':'C','CYM':'C','ASN':'N','GLN':'Q',
           'TRP':'W','GLY':'G','ASX':'B','GLX':'Z','UNK':'X'}
one2three={'A':'ALA','V':'VAL','F':'PHE','I':'ILE','L':'LEU','P':'PRO',
           'M':'MET','D':'ASP','E':'GLU','K':'LYS','R':'ARG','S':'SER',
           'T':'THR','Y':'TYR','H':'HIS','C':'CYS','N':'ASN','Q':'GLN',
           'W':'TRP','G':'GLY','B':'ASX','Z':'GLX','X':'UNK'}

def is_prot_seq(seq):
    if re.compile('[^A-Z]').search(seq):
        return False
    return True

def isMethyl(resname,name):
    """check if an atom belongs to a methyl group"""
    methyls={'ALA':' CB  HB1  HB2  HB3 ',
             'VAL':' CG1 1HG1 2HG1 3HG1 CG2 1HG2 2HG2 3HG2',
             'LEU':' CD1 1HD1 2HD1 3HD1 CD2 1HD2 2HD2 3HD2',
             'ILE':' CG2 1HG2 2HG2 3HG2 CD1 1HD1 2HD1 3HD1',
             'THR':' CG2 1HG2 2HG2 3HG2',
             'MET':' CE HE1 HE2 HE3 '
             }
    res = resname
    if len(res)==1: res = one2tree( res )
    if resname in methyls.keys():
        if name in methyls[resname]:
            return True
    return False
