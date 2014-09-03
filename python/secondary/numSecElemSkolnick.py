#!/usr/bin/python
import sys,os
from secondary.numSecElemSkolnick_f77 import nsec2
from utilities.codedir import CAdir
from utilities.small_utilities import chomp
from inputArgs.inputArgs import inpHand

def num_elem_sec(pdbfile='',header=''):
    if header:
        pdbfile=CAdir+'/'+header+'.pdb'
        if not os.path.exists(pdbfile):
            sys.stderr.write('ERROR (from numSecElemSkolnick.py): '+pdbfile+' does not exists!\n')
            sys.exit(1)
    return nsec2(pdbfile)
    
#execute as standalone program
if __name__=='__main__':

    ih=inpHand('Usage: numSecElemSkolnick.py [arguments]',
               '  -a __header five (or four) letter pdb header',
               '  -b _A_listh file with list of headers. Use in place of -a flag',
               '  -c _A_pdbf single pdb file (in place of -a and -b flag)',
               '  -d _A_outf output file, otherwise output directed to STDOUT'
               )
    ih.parse(locals(),sys.argv)

    #direct output to file or to STDOUT
    if(outf): out=open(outf,'w')
    else: out=sys.stdout

    list=()
    if(header):
        list=(header,) #create a one-element list
    elif(listh):
        list=open(listh,'r').readlines()
        chomp(list) #remove last line character if character is '\n'
    
    if list:
        for h in list:
            line="%s %2i\n" % (h,num_elem_sec(header=h)) #format the output
            out.write(line)
    else:
        line="%2i\n" % (num_elem_sec(pdbfile=pdbf)) #format the output
        out.write(line)
