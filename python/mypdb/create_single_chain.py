#!/usr/bin/python

import os,sys,re
from utilities.codedir import codedir,pdbdir,pdbnonhom,CAdir
from inputArgs.inputArgs import inpHand
from jobs.job import pastry

def listCAs(pin,connect=False):
    '''
    return a list containing the lines of the pdb with CA atoms only
    '''
    cas=[]
    re1=re.compile('TER|END')
    line=pin.readline() #;print line
    while line:
        if re1.match(line): break #reached end of the chain
        if line[0:4]=='ATOM' and line[12:16]==' CA ': cas.append(line)
        line=pin.readline()
    if connect:
        con=['TER\n']
        for i in range(1,len(cas)):
            prev = int( cas[i-1][6:11] )
            curr = int( cas[ i ][6:11] )
            con.append( 'CONECT%5d%5d\n' % ( prev, curr ) )
        cas += con
    return cas


def create_single_chain(header): #return list of CA atoms
    cafile=CAdir+'/'+header+'.pdb'  #;print cafile;sys.exit(1)
    if os.path.exists(cafile):
        pt=open(cafile,'r')
        calines=pt.readlines()
    else:
        pastry(codedir+'/cpp/pdb_utils/bin/output_single_chain.x '+pdbdir+'/pdb'+header[0:4]+'.ent '+header[-1]+' junkall.pdb')
        pastry(codedir+'/cpp/pdb_utils/bin/createCAchain.x junkall.pdb junkca.pdb')
        calines=open('./junkca.pdb','r').readlines()
        os.system('/bin/rm -f junkall.pdb junkca.pdb')
    return map( (lambda x:x[:-1]),calines) #remove last new-line character

if __name__=='__main__':
    ih=inpHand('Usage: create_single_chain.py -i -o -f',
           ' -i _R_header (five-letter code)',
           ' -o _A_outd output directory (default=current dir)',
           ' -f __cafile name of output directory (default=\'CA\')',
           ' -h __help this script returns a file with only CA atoms with sole input of a header file'
           )

    ih.parse(locals(),sys.argv)
    if not outd: outd=os.getcwd()
    if not cafile: cafile='CA'
    calines=map( (lambda x:x+'\n'),create_single_chain(header))
    if not calines: sys.exit(1)
    outca= outd+'/'+cafile
    open(outd+'/'+cafile,'w').writelines(calines)
    sys.exit(0)
