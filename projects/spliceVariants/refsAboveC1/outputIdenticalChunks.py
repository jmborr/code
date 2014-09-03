#!/usr/bin/python

import sys,os,re
from seq.alignYangManageResults import gappedSeq,alignment,alignOut
from spicker.spickerYangResultsManager import spickOut
from utilities.small_utilities import chomp
from secondary.pdb2ss import outputSec

"""

For all pairwise alignments between the variants, find the aligned
chunks bigger than 21 residues that are identical and output their TM
score. The point is to find identical aligned sequences with low TM
score, indicating same sequences folding to very different folds. The
TM score is calculated by pulling the top models for each sequence,
from previous TASSER simulations. In addition, we use these models to
estimate the secondary structure using pdb2ss.py, and output the
\"alignment\" of the secondary structure.

For example, for file
/gpfs1/active/jose/code/projects/spliceVariants/refsAboveC1/O00584/O00584_1_2.alignment.dat,
it will create
/gpfs1/active/jose/code/projects/spliceVariants/refsAboveC1/O00584/O00584_1_2.alignment.chunks.dat
and several
/gpfs1/active/jose/code/projects/spliceVariants/refsAboveC1/O00584/O00584_1_X1_X2_2_Y1_Y2.sup,
where we structurally superimpose by TM score amino acis [X1,X2] of
O005841 with [Y1,Y2] of O005842

"""

dir='/gpfs1/scratch/jose/spliceVariants'
listf=dir+'/refsAboveC1/analysis/all_pairwise_alignments.list'
#listf=dir+'/refsAboveC1/analysis/toy.list'

remarks0='REMARK from the TM score structural alignment\n'

res={} ; res[1]=re.compile('_(\d+)_') ; res[2]=re.compile('_(\d+)\.')
for file in chomp(open(listf,'r').readlines()):
    outfile=dir+'/refsAboveC1/'+file[0:-4]+'.chunks.dat' #;print outfile;sys.exit(1)
    outsuproot=dir+'/refsAboveC1/'+file[2:15] #prefix for the superposition PDB files
    ids={}
    #an instance of variable "file" will be like "./O00584/O00584_1_2.alignment.dat"
    ids[1]=res[1].search(file).group(1)
    ids[2]=res[2].search(file).group(1)
    out=open(outfile,'w')
    var_root=file[2:8] ;
    #headers for the splice variants
    vars={} ; vars[1]=var_root+file[16] ; vars[2]=var_root+file[18]
    #find top model for each variant, and estimate secondary structure
    models={} ; sss={}
    for id in vars.keys():
        var=vars[id]
        s=spickOut(dir=dir+'/tasser/'+var[1]+'/'+var,rebuilt='.refined.rebuilt')
        models[id]=s.combo[s.densest] #name of the top combo file
        sss[id]=outputSec(models[id])
    #obtain an alignment object
    s=alignOut(dir+'/refsAboveC1/'+file).alg
    #extract identical chunks from the alignment with size equal or above 21 residues
    chunks=s.extractIdenticalChunks(minSize=21)
    #find the corresponding sequence segments
    for chunk in chunks:
        spans=chunk.parentUngappedSpans()
        segments={}
        segments[1]=`1+spans[1][0]`+'-'+`1+spans[1][1]`
        segments[2]=`1+spans[2][0]`+'-'+`1+spans[2][1]`
        out.write('segment '+segments[1]+' of '+vars[1]+' is aligned to segment ')
        out.write(segments[2]+' of '+vars[2]+'\n')
        #reformat spans into a single line. Remember the shift in the index convention
        sg=segments[1]+','+`1+spans[2][0]`+'-'+`1+spans[2][1]`
        #find TM score of the corresponding fragments in the variants
        outsup=outsuproot+'_'+ids[1]+'_'+segments[1]+'_'+ids[2]+'_'+segments[2]+'.pdb'
        #calculate TM score of the corresponding sequence fragments, and output the superposition
        cmd='tm_score.py -a '+models[1]+' -b '+models[2]+' -g '+sg+' -d yes -c '+outsup
        remarks=remarks0+''.join( os.popen(cmd).readlines() )
        #superimpose the secondary structures as an alignment object
        a=gappedSeq(sss[1][spans[1][0]:1+spans[1][1]])
        b=gappedSeq(sss[2][spans[2][0]:1+spans[2][1]])
        ssalignment=alignment(a,b)
        #print ssalignment.gs[1].gs;print ssalignment.gs[2].gs;sys.exit(1)
        out.write( 'PRIMARY\n'+chunk.info()+'\n'+remarks+'SECONDARY\n'+ssalignment.info()+'\n')
        #print 'PRIMARY\n'+chunk.info()+'\n'+remarks+'SECONDARY\n'+ssalignment.info()
    out.close()
        
