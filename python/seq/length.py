#!/usr/bin/python
import sys,os
from inputArgs.inputArgs import inpHand,addAbsPath
from utilities.small_utilities import three2one,chomp #three2one dict for aa codes

#function to calculate length of a protein related to a four (or five) pdb code
def length_pdb(pdbf):
    #go first to appropriate line in pdb file
    pdbhandle=open(pdbf,'r')
    line=pdbhandle.readline()
    while line[0:5]!='ATOM ':
        line=pdbhandle.readline()
    #count number of residues
    n=0
    while line[0:3]!='TER ' and line[0:3]!='END ' and line:
        if line[0:5]=='ATOM ' and line[12:16]==' CA ':
            n=n+1  #;print n
        line=pdbhandle.readline()
    pdbhandle.close()
    return n

        
def length_header(header):
    pdbf='/library/pdb/pdb'+header[0:4]+'.ent'
    #print pdbf ; sys.exit(0) ;
    if len(header)<5: chId='_'
    else: chId=header[-1].upper(); #print chId ; sys.exit(0)
    if not os.path.exists(pdbf): ih.abort(' !!!Non existent PDB file '+pdbf)
    if chId=='_':
        for line in open(pdbf,'r').readlines():
            if line[0:5]=='ATOM ':
                chId=line[21]
                break
    #go to appropriate line in pdb file
    pdbhandle=open(pdbf,'r')
    line=pdbhandle.readline()
    while line[0:5]!='ATOM ' or line[21]!=chId:
        line=pdbhandle.readline()
    #print line ; sys.exit(0)
    #count number of residues
    n=0
    while line[0:4]!='TER ' and line[0:4]!='END ':
        #print line
        if line[0:5]=='ATOM ' and line[12:16]==' CA ':
            n=n+1
        line=pdbhandle.readline()
    pdbhandle.close()
    return n

#execute as standalone program
if __name__=='__main__':

    ih=inpHand('Usage: length.py  -i -l -o',
               '  -i __header five (or four) letter pdb header',
               '  -l _A_listh file with list of headers. Use in place of -i flag',
               '  -o _A_outf output file, otherwise output directed to STDOUT'
               )
    ih.parse(locals(),sys.argv)

    #did we pass one or more headers?
    if(header):
        list=(header,) #create a one-element list
    else:
        list=open(listh,'r').readlines()
        chomp(list) #remove last line character if character is '\n'

    #direct output to file or to STDOUT
    if(outf): out=open(outf,'w')
    else: out=sys.stdout
    
    for h in list:
        line="%s %4i\n" % (h,length_header(h)) #format the output
        out.write(line)
