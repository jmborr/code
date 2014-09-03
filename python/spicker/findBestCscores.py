#!/usr/bin/python

import sys,os
from utilities.small_utilities import chomp
from utilities.returnVal import returnVal
from spicker.spickerYangResultsManager import spickOut
from inputArgs.inputArgs import inpHand,deglobb

rv=returnVal() #object to store results from findBestCscore, including error status and messages

def printHeader():
    return '#header, bestCscore, related Z, id, Highest Z'

def findBestCscore(dir,templatef):
    rv.clearErr()
    if not os.path.exists(dir):
        rv.setErr(msg='no '+dir+' directory found')        
        return rv
    if not os.path.exists(templatef):
        rv.setErr(msg='no '+templatef+' file found')
        return rv
    s=spickOut(dir=dir)
    if s.readError:
        rv.setErr(msg='error reading spicker output in '+dir)
        return rv
    #rank cluster ID's according to decreasing Cscore
    ranked=s.rankIDsByCscores(templatef=templatef)
    bestID=ranked[0]
    rv.val='%6.2f %6.2f %2d %6.2f'%(s.Cscores[bestID],s.ZscoreOfCluster[bestID],bestID,s.Zscores[0])
    return rv

if __name__=='__main__':

    ih=inpHand('Usage: findBestCscores.py\n',
               ' -a __header header',
               ' -b __headerList list of headers (in place of -a)',
               ' -c _RA_dir (possibly x-globbed) directory with spicker output',
               ' -d _RA_templf (possibly x-globbed) template.dat file)',
               ' -e _A_outf output file (def:STDOUT)'
               )
    ih.parse(locals(),sys.argv)
    if not header and not headerList: ih.abort('ERROR:supply either -a or -b flag!')
    if not outf: out=sys.stdout
    else: out=open(outf,'w')
    
    if header: headers=[header,]
    else: headers=chomp(open(headerList,'r').readlines())
    out.write(printHeader()+'\n')
    for header in headers:
        dir2=deglobb(dir,header) ; templatef=deglobb(templf,header) #;print dir2+'\n'+templatef
#        out.write(header+' '+findBestCscore(dir2,templatef)+'\n')
        findBestCscore(dir2,templatef)
        if rv.hasErr(): out.write(header+' '+rv.errMsg+'\n')
        else: out.write(header+' '+rv.val+'\n')
    sys.exit(0)
