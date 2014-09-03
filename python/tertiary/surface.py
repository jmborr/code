#!/usr/bin/python

#Accessible surface area
#(calculated for the residue X in the tripeptide G-X-G)
#S.Miller et al, J. Mol. Biol., 196:641-656 (1987)
area1={
'A':   113.,
'R':   241.,
'D':   151.,
'N':   158.,
'C':   140.,
'E':   183.,
'Q':   189.,
'G':    85.,
'H':   194.,
'I':   182.,
'L':   180.,
'K':   211.,
'M':   204.,
'F':   218.,
'P':   143.,
'S':   122.,
'T':   146.,
'W':   259.,
'Y':   229.,
'V':   160. 
}

area3={
'ALA':  113.,
'ARG':	241.,
'ASP':	151.,
'ASN':	158.,
'CYS':	140.,
'GLU':	183.,
'GLN':	189.,
'GLY':	 85.,
'HIS':	194.,
'ILE':	182.,
'LEU':	180.,
'LYS':	211.,
'MET':	204.,
'PHE':	218.,
'PRO':	143.,
'SER':	122.,
'THR':	146.,
'TRP':	259.,
'TYR':	229.,
'VAL':  160. 	
}

area1sc={ #same as area1, but including only the side chain
'A':    67.,
'R':   196.,
'D':   106.,
'N':   113.,
'C':   104.,
'E':   138.,
'Q':   144.,
'G':     0.,
'H':   151.,
'I':   140.,
'L':   137.,
'K':   167.,
'M':   160.,
'F':   175.,
'P':   105.,
'S':    80.,
'T':   102.,
'W':   217.,
'Y':   187.,
'V':   117. 
}

area3sc={
'ALA':   67.,
'ARG':	196.,
'ASP':	106.,
'ASN':	113.,
'CYS':	104.,
'GLU':	138.,
'GLN':	144.,
'GLY':	  0.,
'HIS':	151.,
'ILE':	140.,
'LEU':	137.,
'LYS':	167.,
'MET':	160.,
'PHE':	175.,
'PRO':	105.,
'SER':	 80.,
'THR':	102.,
'TRP':	217.,
'TYR':	187.,
'VAL':  117. 	
}
