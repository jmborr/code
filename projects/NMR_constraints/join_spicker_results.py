#!/usr/bin/python
import pdb
import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import chomp
from jobs.job import gen_lf,pastry

ih=inpHand('Usage: join_spicker_results.py',
           ' -a _A_ind root directory ("summary*.txt" files inside, def=curr dir)',
           ' -c __outf output filename  (def=summary.all)',
           ' -d _A_outd output directory (def=curr dir)',
           ' -e __slack permissivity between first and top cluster(0<=slack<1, def=0.0)',
           ' -f _A_listf filter by list of headers',
           ' -g __logname a logfile containing the passed word must be in the same directory as the file summary*.txt',
           ' -i __rco rco-cutof,f, to filter out those summary.txt files with rco-cutoff'
           )
ih.parse(locals(),sys.argv)

#set defaults
if not ind: ind=os.getcwd()             #;print ind
if not outd: outd=os.getcwd()
if not outf: outf=outd+'/summary.all'
else       : outf=outd+'/'+outf
if not slack:slack=1.0
else        : slack=1-float(slack)
listheaders=''
if listf: listheaders=chomp(open(listf,'r').readlines())
if logname:
    p=re.compile(logname)
    #print logname

pt_out=open(outf,'w')
pt_out.write('#header most-populated-centroid best-rmsd-centroid most-dense-centroid\n');
#print 'find '+ind+' -name "summary*.txt"'; sys.exit(1)
for file in chomp(os.popen('find '+ind+' -name "summary*.txt"').readlines()):
    #pdb.set_trace()
    if logname:
        basedir=os.path.dirname(file) #directory where "summary*.txt" resides 
        logfile=os.popen('ls '+basedir+'/*.log').readline()
        if not p.search(logfile): continue #don't process if required logname not present
    pt=open(file,'r')   #;print 'file=',file
    line=pt.readline()  #;print 'line=',line
    #set header
    if re.compile('Target').search(line):
        header=line[8:13] #first line contains the header
        print header,file
    else: header='xxxxx'
    #skip rest if target not in passed list of headers
    if listheaders:
        if header not in listheaders:
            pt.close()
            continue
    #print file
    #print logfile
    p=re.compile('rank size percent density')
    while not p.search(line):
        line=pt.readline()
        #print line

    #find most dense cluster
    line=pt.readline()
    (denserank,size,percent,topdensity,avrmsd)=line.split()
    topdensity=float(topdensity)
    p=re.compile('\s+\d+\s+\d+\s+\d+')
    line=pt.readline()
    while p.match(line):
        (rank,size,percent,density,avrmsd)=line.split()
        density=float(density)  #;print 'density=',density
        if density>topdensity:
            #print density,'>',topdensity
            denserank=rank
            topdensity=density
        line=pt.readline()
    #print 'denserank=',denserank,'topdensity=',topdensity
    p=re.compile('#cluster')
    #find fist line of centroids info
    while not p.search(line): line=pt.readline()
    #find top cluster
    line=pt.readline()
    (firstrank,firstrmsd)=line.split()
    firstrmsd=float(firstrmsd)
    if firstrank==denserank: densermsd=firstrmsd
    toprank=firstrank    
    toprmsd=firstrmsd
    p=re.compile('\s+\d+\s+\d+')    
    line=chomp(pt.readline())
    while p.match(line):
        (rank,rmsd)=line.split()
        rmsd=float(rmsd)
        if rank==denserank: densermsd=rmsd
        if rmsd<firstrmsd*slack and rmsd<toprmsd:
            toprank=rank
            toprmsd=rmsd
        line=chomp(pt.readline())
    pt.close()
    #write output to file    
    line=header+'  1 %5.2f %3s %5.2f %3s %5.2f\n'%(firstrmsd,toprank,toprmsd,denserank,densermsd)
    pt_out.write(line)
sys.exit(0)
