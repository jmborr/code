#!/usr/bin/python

import os,sys

from utilities.codedir import codedir
from utilities.small_utilities import Bye,chomp
from inputArgs.inputArgs import inpHand,deglobb
from jobs.job import job,pastry,genJobList

root_code='/gpfs1/active/jose/code/projects/spliceVariants'
root_scratch='/gpfs1/scratch/jose/spliceVariants/out3'

joblist=genJobList()


for line in open(os.path.join(root_scratch,'summary.all'),'r').readlines():    
    if line[0]=='#': continue #is a comment line
    header,type,L,cov,id=line.split()[0:5] ; id=int(id)
    outdir=os.path.join(root_scratch,header[:6],header[6:])
    model=os.path.join(outdir,'combo%02d.pdb.rebuilt'%(id)) #;Bye(model)
    job(exe='runmin2.sh',exed=root_code,args=model+' '+outdir).qsub('min'+header,outdir,wallt='0.082',mem_limit='150',joblist=joblist,extraflags='-Wx=qos:personal')
    
sys.exit(0)
