#!/usr/bin/python
from jobs.job import job

j=job(name='list of headers, tasser + spicker, but no prospector',
      exe='preditc_strutc_list.py',
      exed='/gpfs1/active/jose/code/python/combo_jobs',
      args='-u rsub -l /gpfs1/active/jose/projects/science_set/run/list -t "codedir/python/tasser/tasser1.0/tasser1.0_header.py -o . -a CA -i" -z "codedir/c/spicker/spickerZZ.x -l tra.in -s seq.dat -a CA -o . -b NOTAR -m" -m "CA *pdb *dat in.dd out.d pair* rep* tra.in rmsinp *.fasta summary.txt" -o /gpfs1/active/jose/projects/science_set/run',
      story='submit a list of five-letter headers, then run tasser with no\ninfo from prospector. After this run spickerZZ.'
      )
j.store()

