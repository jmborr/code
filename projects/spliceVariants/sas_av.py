#!/usr/bin/python

"""Find average amino acid composition of domain surface with
$SC/db/scop/scop.1.7/sas.dat, which lease for every entry of scop list
of domains (release 1.71) filtered at 40% sequence identity the
solvent accesible surface of each amino acid

Results are:
A 0.043
C 0.005
D 0.073
E 0.113
F 0.029
G 0.041
H 0.029
I 0.030
K 0.119
L 0.056
M 0.017
N 0.054
P 0.050
Q 0.057
R 0.101
S 0.050
T 0.048
V 0.036
W 0.013
Y 0.036
"""

import os,sys
from seq.letters import valid_A,three2one

comp={}
for X in valid_A: comp[X]=0.0
total_area=0.0
for l in os.popen('grep -v -e ">" /gpfs1/scratch/jose/db/scop/scop.1.71/sas.dat',
                  'r').readlines():
    s=float(l[13:])
    comp[ three2one[l[0:3]] ]+=s
    total_area+=s

for X in valid_A:
    print X,'%5.3lf'%(comp[X]/total_area)



sys.exit(0)

