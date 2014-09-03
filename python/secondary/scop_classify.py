#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.lists import read_file_lines
from utilities.codedir import scratchdir,scopd

scopf=os.path.join(scopd,'dir.cla.scop.txt_1.71')
pattern=re.compile('\s([a-k])\.')
verbosedict={'a':'All alpha proteins',
             'b':'All beta proteins',
             'c':'Alpha and beta proteins (a/b)',
             'd':'Alpha and beta proteins (a+b)',
             'e':'Multi-domain proteins (alpha and beta)',
             'f':'Membrane and cell surface proteins and peptides',
             'g':'Small proteins',
             'h':'Coiled coil proteins',
             'i':'Low resolution protein structures ',
             'j':'Peptides',
             'k':'Designed proteins',
             }

def scop_classify(headers):
    '''given a list of 5-letter headers, return a list with SCOP
    classifications for each header'''    
    scopls=read_file_lines(scopf)
    classifications=[]
    for header in headers:
        matchl=''
        #match only PDB ID
        root=header[0:4]
        rootmatches=[]
        for scopl in scopls:
            if root in scopl: rootmatches.append(scopl)
        #also match chain ID
        header2=header.lower() #scop uses lowercase for chain ID
        for rootmatch in rootmatches:
            if header2 in rootmatch:
                matchl=rootmatch
                break
        if not matchl and rootmatches: matchl=rootmatches[0] #if we could not find chain match
        if not matchl:
            sys.stderr.write('scop_classify could not find a match in SCOP for '+header+'\n')
            classifications.append( [header,'',''] )
        else:
            scopclass=pattern.search(matchl).group(1)
            verbose=verbosedict[scopclass]
            classifications.append( [header,scopclass,verbose] )            
    return classifications


if __name__=='__main__':
    inp=inpHand('Usage: scop_classify.py [options]',
                ' -a __header 5-letter header',
                ' -b _A_listf  list of 5-letter headers (in place of -a flag)',
                ' -c _A_outf output file (def: STDOUT)',
                )
    inp.parse(locals(),sys.argv)

    if not header and not listf: inp.abort('Either supply a header or a list of headers')
    if header: headers=[header,]
    if listf: headers=read_file_lines(listf)
    if not outf: outp=sys.stdout
    else: outp=open(outf,'w')

    buf=''
    for item in scop_classify(headers): buf+=' '.join(item)+'\n'        
    outp.write( ''.join(buf) )
    
    sys.exit(0)
