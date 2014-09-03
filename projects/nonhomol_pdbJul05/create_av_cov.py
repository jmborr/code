#!/usr/bin/python

from utilities.readingWritingFiles import read_column,write_lists
from utilities.lists import multiply
from stats.statistics import running_average


l=read_column('/gpfs1/active/jose/projects/nonhomol_pdbJul05/alpha_beta/coverages.all',3,isFloat=1)
x=read_column('/gpfs1/active/jose/projects/nonhomol_pdbJul05/alpha_beta/coverages.all',15,isFloat=1)

#v=multiply(l,x)
[avxl,avyl,stdl,stdml]=running_average(l,x,dx=50)
#[avxl,avyl,stdl,stdml]=running_average(x,l,dx=0.1)
write_lists('/gpfs1/active/jose/projects/nonhomol_pdbJul05/alpha_beta/junk',[avxl,avyl,stdl,stdml])
