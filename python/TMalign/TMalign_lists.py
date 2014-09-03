#!/usr/bin/python

import sys,os
from utilities.small_utilities import Bye
from time import sleep
from jobs.job import pastry

"""Given a list of conformations (PDB file format), we want to
construct the matrix of TM-scores among them. We use fortran
executable TMalign_lists, but first, split the list into chunks and
send all pair combinations of the chunks as jobs to the cluster.
Example: TMalign_lists.py -a /tmp/jose/struct:/gpfs1/scratch/foo/struct.tbz2 -b list -c outd"""

def TMalign_lists(lib,listf,outd,prefix='',exe_type='yang02',sizechunk=250,wallt=0.259):

    """TMalign of a members of a library of structures"""
    
    from jobs.job import job,genJobList

    #settings for exe_type='yang02'
    exe='TMalign_lists'
    exed='/gpfs1/active/jose/code/f77/TM-align/TMalign2.0'
    if exe_type=='shashi':
        exe='TMalign_lists'
        exed='/gpfs1/active/jose/code/f77/TM-align/Shashi'
        if not sizechunk: sizechunk=100
        if not wallt: wallt=0.999
    elif exe_type=='fragment':
        exe='TMalign_lists'
        exed='/gpfs1/active/jose/code/f77/TM-align/Shashi/fragment'
        if not sizechunk: sizechunk=100
        if not wallt: wallt=0.99
    if not prefix: prefix='tm'
    structDir=lib.split(':')[0]
    joblist=genJobList()

    #create chunks and list files
    listL=open(listf,'r').readlines()
    basechunks=[]
    chunks=[]
    nchunk=0
    while listL:
        chunk=outd+'/list%03d'%(nchunk)
        chunks.append(chunk)
        basechunks.append('list%03d'%(nchunk))
        open(chunk,'w').write(''.join(listL[0:sizechunk]))
        nchunk+=1
        listL=listL[sizechunk:]

    #dispatch jobs
    ibatch=0
    for i in range(0,nchunk-1):
        header='%03d_%03d'%(i,i) ; jobname=prefix+header ; outfile=jobname
        outcmd='/bin/mv '+outfile+' '+outd+'/'+jobname
        args=args=' -sdir '+structDir+' -ss '+basechunks[i]+' -outf '+outfile+' -simpl'
        Job=job(exe=exe,exed=exed,args=args,shared=lib,outcmd=outcmd,inputs=[chunks[i],])
        Job.qsub(jobname,outd,wallt=wallt,mem_limit='150',ddisk='500',joblist=joblist)
        for j in range(i+1,nchunk):
            header='%03d_%03d'%(i,j) ; jobname=prefix+header ; outfile=jobname            
            outcmd='/bin/mv '+outfile+' '+outd+'/'+jobname
            args=' -sdir '+structDir+' -ss '+basechunks[i]+' -tt '+basechunks[j]+\
                  ' -outf '+outfile+' -tdir '+structDir+' -simpl'
            Job=job(exe=exe,exed=exed,args=args,shared=lib,\
                    outcmd=outcmd,inputs=[chunks[i],chunks[j]])
            Job.qsub(jobname,outd,wallt=wallt,mem_limit='150',ddisk='500',joblist=joblist)
            ibatch=Job.dormant(ibatch)
    #deal with very last job
    header='%03d_%03d'%(nchunk-1,nchunk-1) ; jobname=prefix+header ; outfile=jobname
    outcmd='/bin/mv '+outfile+' '+outd+'/'+jobname
    Job=job(exe=exe,exed=exed,args=' -sdir '+structDir+' -ss '+basechunks[-1]+\
            ' -outf '+outfile+' -simpl',shared=lib,outcmd=outcmd,inputs=[chunks[-1],])
    Job.qsub(jobname,outd,wallt=wallt,mem_limit='150',ddisk='500',joblist=joblist)


def TMalign_listsII(lib,listf,lib2,listf2,outd,prefix='',exe_type='yang02',sizechunk=100,\
                    wallt=None):

    """TMalign of a library against other library"""
    
    from jobs.job import job,genJobList

    #settings for exe_type='yang02'
    exe='TMalign_listsII'
    exed='/gpfs1/active/jose/code/f77/TM-align/TMalign2.0'
    if not sizechunk: sizechunk=250
    if not wallt: wallt=0.249
    if exe_type=='shashi':
        exe='TMalign_lists'
        exed='/gpfs1/active/jose/code/f77/TM-align/Shashi'
        if not sizechunk: sizechunk=100
        if not wallt: wallt=0.999
    elif exe_type=='fragment':
        exe='TMalign_lists'
        exed='/gpfs1/active/jose/code/f77/TM-align/Shashi/fragment'
        if not sizechunk: sizechunk=138
        if not wallt: wallt=0.999

    joblist=genJobList()
    #create chunks and list files
    if not prefix: prefix='tm'
    structDir=lib.split(':')[0]
    listL=open(listf,'r').readlines()
    basechunks=[] ; chunks=[] ; nchunk=0
    while listL:
        chunk=outd+'/list%03d'%(nchunk)
        chunks.append(chunk)
        basechunks.append('list%03d'%(nchunk))
        open(chunk,'w').write(''.join(listL[0:sizechunk]))
        nchunk+=1
        listL=listL[sizechunk:]

    structDir2=lib2.split(':')[0]
    listL2=open(listf2,'r').readlines()
    basechunks2=[] ; chunks2=[] ; nchunk2=0
    while listL2:
        chunk=outd+'/list%03d.2'%(nchunk2)
        chunks2.append(chunk)
        basechunks2.append('list%03d.2'%(nchunk2))
        open(chunk,'w').write(''.join(listL2[0:sizechunk]))
        nchunk2+=1
        listL2=listL2[sizechunk:]

    #dispatch jobs
    ibatch=0
    for i in range(0,nchunk):
        for j in range(0,nchunk2):
            header='%03d_%03d'%(i,j) ; jobname=prefix+header ; outfile=jobname            
            outcmd='/bin/mv '+outfile+' '+outd+'/'+jobname
            args=' -sdir '+structDir+' -tdir '+structDir2+' -ss '+basechunks[i]+\
                  ' -tt '+basechunks2[j]+' -outf '+outfile+' -simpl'
            Job=job(exe=exe,exed=exed,args=args,shared=lib+','+lib2,\
                    outcmd=outcmd,inputs=[chunks[i],chunks2[j]])
            Job.qsub(jobname,outd,wallt=wallt,mem_limit='150',ddisk='500',joblist=joblist)
            ibatch=Job.dormant(ibatch)
            

def gathertmAll2Alldat(tmAll2Alldir,tmAll2Allfile,listofstruct,tmtmregex='tm??_??'):

    """create file tmAll2All.dat file, remove tmAll2All directory"""

    cmd='cat '+tmAll2Alldir+'/'+tmtmregex+' > '+tmAll2Allfile
    os.system(cmd)
    nlines=int( os.popen('wc -l '+tmAll2Allfile).readline().split()[0] )
    N=int( os.popen('wc -l '+listofstruct).readline().split()[0] )
    trueNlines=N*(N-1)/2
    if nlines!=trueNlines: sys.stderr.write('ERROR: '+`nlines`+' versus correct '+`trueNlines`)
    #os.system('/bin/rm -r '+tmAll2Alld) #remove directory, we don't need it anymore
        

def gentmX2Allfiles(tmAll2Allfile,tmX2alldir,headers,tmco=0.4):

    """generate histogram files

    Keep a buffer for every histogram file, which we flush if it exceeds certain size"""

    from utilities.small_utilities import chomp,Bye,junkName
    
    #init buffer to hold hits
    buffer={}
    lengths={}
    for header in headers:
        buffer[header]=[]
        lengths[header]=0

    #local directory to avoid gpfs traffic
    localoutd='/tmp/jose/'+junkName()
    if os.path.exists(localoutd):
        sys.stderr.write('ERROR: local directory '+localoutd+' exists\n')
        sys.exit(1)
    os.system('/bin/mkdir -p '+localoutd)

    pin=open(tmAll2Allfile,'r')
    line=pin.readline()
    n=0
    while line:
        n+=1
        if len(line)<2:
            line=pin.readline()
            continue
        items=line.split()    
        if len(items)!=4:
            line=pin.readline()
            continue
        header1,header2,tm1,tm2=items #tm_i normalized by length of header_i
        if header1 in headers and header2 in headers:
            tmf=float(tm1)
            if tmf>tmco:
                buffer[header1].append(header2+' '+tm1)
                lengths[header1]+=1
            tmf=float(tm2) 
            if tmf>tmco:
                buffer[header2].append(header1+' '+tm2)
                lengths[header2]+=1
            #flush buffer component
            for headeri in (header1,header2):
                listx=buffer[headeri]
                if len(listx) > 200: #store up to 200 lines
                    open(localoutd+'/'+headeri+'.dat','a').write( '\n'.join(listx)+'\n')
                    buffer[headeri]=[]
        line=pin.readline()
    
    for header in headers:
        f=localoutd+'/'+header+'.dat'
        #flush buffer
        listx=buffer[header]
        open(f,'a').write( '\n'.join(listx)+'\n')
        buffer[header]=[]
        #prepend number of hits
        n=lengths[header]
        buff=''.join( open(f,'r').readlines() )
        open(f,'w').write(`n`+'\n'+buff)

    if not os.path.exists(tmX2alldir): os.system('/bin/rm -r '+tmX2alldir)
    #move all files. Can't do mv localoutd/* because "argument list is too long". Can't
    #send each file one by one because of heavy gpfs1 traffic. Solution is to "tar" first
    for header in headers:
        f=localoutd+'/'+header+'.dat'
        pastry('tar uf junk.tar '+f)
    pastry('tar xf junk.tar -C '+tmX2alldir)
    pastry('/bin/mv '+localoutd+' '+tmX2alldir) #move local stuff


#execute as standalone script
if __name__=='__main__':
    
    from inputArgs.inputArgs import inpHand

    ih=inpHand('Usage: TMalign_lists.py [options]',
               '  -a _R_lib library line with compressed structure directory ',
               '  -b _AR_listf file list of structures',
               '  -c _AR_outd output directory',
               '  -d __sizechunk split list in chunks of specified size (def:1000,max:1000)',
               '  -e __prefix prefix for jobname (def: tm)'
               )
    ih.parse(locals(),sys.argv)

    #init variables
    if not sizechunk: sizechunk=1000
    else: sizechunk=int(sizechunk)
    if not prefix: prefix='tm'
    if not re.compile('[a-z,A-Z]').search(prefix[0]): prefix='a'+prefix[1:]

    TMalign_list(lib,listf,outd,sizechunk=sizechunk,prefix=prefix)
    
    sys.exit(0)

    
