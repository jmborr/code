#!/usr/bin/python

import sys,os,re
from inputArgs.inputArgs import inpHand
from utilities.small_utilities import Bye
from amber.amber9 import ptraj

codedir=os.environ['CODED']
isthere=os.path.exists
joink=os.path.join
basename=os.path.basename

ih=inpHand('Usage: qsubpmemd.py',
           ' -a _RA_prodd absolute path to production directory',
           ' -b __service type of service (def:probe)',
           ' -c _R_TYP protein ID',
           ' -d _A_infile0 input file',
           ' -e _A_topfile0 topology file',
           ' -f _A_initrstfile0 restart file',
           ' -g _A_reffile0 reference file',
           ' -i __maxframe integer maxframe',
           ' -j __stripline line to strip atoms (def: "WAT,Na+")',
           ' -l __maxtime maxtime in nanoseconds, in place of maxframe'
           )
ih.parse(locals(),sys.argv)

#initialize default variables
if not service: service='probe'
if not stripline: stripline="WAT,Na+"
striplist=stripline.split(",")
#assumed certain file names
currd=joink(prodd,'curr')
prevd=joink(prodd,'prev')
debugd=joink(prodd,'debug')
initd=joink(prodd,'init')

arrangefile0=joink(codedir,'software/amber9/prodprotoc/Prod/arrange')
arrangefile=joink(prodd,'arrange')
crdfile=joink(currd,TYP+'.crd')
infile=joink(initd,TYP+'.in')
initrstfile=joink(initd,TYP+'.rst')
maxfrfile=joink(prodd,'maxframe.dat')
outfile=joink(currd,TYP+'.out')
submitfile0=joink(codedir,'software/amber9/prodprotoc/Prod/submit')
submitfile=joink(prodd,'submit')
reffile=joink(initd,TYP+'.ref')
refunsf=joink(initd,TYP+'.uns.ref')
rstfile=joink(currd,TYP+'.rst')
rstlistf=joink(prodd,'rst.list')
topfile=joink(initd,TYP+'.top')
topunsf=joink(initd,TYP+'.uns.top')
unstopfile=joink(initd,TYP+'.uns.top')
velfile=joink(currd,TYP+'.vel')

def check_rms(crd,rms):
    """has rms gone too far from unsolvated reference?"""
    rc=4.0 #rmsd cut-off
    ptraj(topunsf,crd).reference(refunsf)['self'].rsm(first='reference',outrms=rms).go()
    rm=0.0
    for line in open(rms,'r').readlines():
        x,r=line.split()  ;  r=float(r)
        if r>rm: rm=r
    if rm>rc:
        subject='RMSD='+str(rm)+' in '+crd
        cmd='echo "" | mail -s '+subject+' borreguero@gmail.com'
        os.system(cmd)
    
def framesnanosec(infile):
    """find number of frames to produce a nanosecond"""
    buf=''.join(open(infile,'r').readlines())
    n=int(re.compile('nstlim=(\d+)').search(buf).group(1))
    dt=float(re.compile(' dt=(\d+\.\d+)').search(buf).group(1)) #ps/step
    print `n*dt`+'picoseconds per simulation'
    
def getMaxframe(infile,maxtime):
    """find total number of frames for a particular maximum time"""
    buf=''.join(open(infile,'r').readlines())
    dt=float(re.compile(' dt=(\d+\.\d+)').search(buf).group(1)) #ps/step
    ntpr=int(re.compile('ntpr=(\d+)').search(buf).group(1)) #steps/frame
    return str(int(float(maxtime)*1000/(dt*ntpr)))
    
def conf_submit(submitfile):
    """configure the submitfile"""
    buf=''.join( open(submitfile0,'r').readlines() )
    buf=buf.replace('_JOBNAME_',TYP)
    buf=buf.replace('_PROD_',prodd)
    buf=buf.replace('_TYP_',TYP)
    buf=buf.replace('_SIZE_','64')
    open(submitfile,'w').write(buf)
    os.system('chmod u+x '+submitfile)

def conf_arrange(arrangefile):
    """configure the arrangefile"""
    buf=''.join( open(arrangefile0,'r').readlines() )
    buf=buf.replace('_JOBNAME_',TYP)
    buf=buf.replace('_PROD_',prodd)
    buf=buf.replace('_TYP_',TYP)
    open(arrangefile,'w').write(buf)
    os.system('chmod u+x '+arrangefile)

if service=='probe':
    """find out what action to take"""
    if not isthere(initd):
        sys.stdout('There is no init directory. Should I initialize?(y/n)')
        doweinit=sys.stdin.readline()
        if doweinit.lower()=='y':
            service='initialize' #start the simuation series
    elif isthere(currd):
        if isthere(crdfile):
            service='arrange'    #reap simulation result
        else:
            os.system('qsub '+submitfile) #send next simulation
    else:
        sys.stdout.write('There is init directory. Should I initialize?(y/n)')
        doweinit=sys.stdin.readline()
        if doweinit.lower()[0]=='y':
            service='initialize' #start the simuation series

if service=='initialize': #start the simulation series
    os.system('/bin/mkdir -p '+initd+' '+debugd)
    if not isthere(infile0) and not isthere(infile): Bye('no '+infile0)
    if infile0 and infile!=infile0:
        buf='/bin/cp '+infile0+' '+infile ; os.system(buf)
    print 'found/created infile '+basename(infile)
    if not isthere(topfile0) and not isthere(topfile): Bye('no '+topfile0)
    if topfile0 and topfile!=topfile0:
        os.system('/bin/cp '+topfile0+' '+topfile)
    print 'found/created topfile '+basename(topfile)
    if not isthere(initrstfile0) and not isthere(initrstfile):
        Bye('no '+initrstfile0)
    if initrstfile0 and initrstfile!=initrstfile0:
        os.system('/bin/cp '+initrstfile0+' '+initrstfile)
    print 'found/created initrstfile '+basename(initrstfile)
    if not isthere(reffile0) and not isthere(reffile): Bye('no '+reffile0)
    if reffile0 and reffile!=reffile0:
        os.system('/bin/cp '+reffile0+' '+reffile)
    print 'found/created reffile '+basename(reffile)
    if not isthere(refunsf): Bye('no '+refunsf)
    if not isthere(topunsf): Bye('no '+topunsf)
    if not maxframe and not maxtime:
        Bye('provide  maxframe or maxtime, please')
    if maxtime and not maxframe: maxframe=getMaxframe(infile,maxtime)
    print 'We will collect '+maxframe+' frames'
    open(maxfrfile,'w').write(maxframe)
    print 'created '+basename(maxfrfile)
    os.system('/bin/touch '+rstlistf)
    os.system('/bin/mkdir -p '+currd)
    os.system('/bin/mkdir -p '+prevd)
    os.system('/bin/cp '+initrstfile+' '+rstfile)
    for file in (crdfile,velfile,outfile):
        if isthere(file): os.system('/bin/rm '+file)
    conf_submit(submitfile)
    conf_arrange(arrangefile)
    ptraj(topfile,reffile).strip(striplist)['self'].go(outcrd=refunsf)
    os.system('qsub '+submitfile)  
elif service=='arrange':  #arrange next simulation
    #get last frame number
    rstlist=open(rstlistf,'r').readlines()
    if not rstlist:
        N=0
    else:
        N=int(re.compile('\.(\d+)\.rst').search(rstlist[-1]).group(1))
    #generate new restart file with next to last frame
    rstfilenew=rstfile+'.new'
    returnval=ptraj(topfile,crdfile,vel=velfile).create_rst(rstfilenew)
    n=returnval['N']
    if not n: Bye('ERROR: empty '+crdfile+' or '+velfile)
    M=N+n  ;print 'n=',n,'N=',N
    crd=joink(prevd,TYP+'.%05d'%(N+1)+'_'+'%05d'%M+'.crd')
    vel=joink(prevd,TYP+'.%05d'%(N+1)+'_'+'%05d'%M+'.vel')
    out=joink(prevd,TYP+'.%05d'%(N+1)+'_'+'%05d'%M+'.out')
    rst=joink(prevd,TYP+'.%05d'%M+'.rst')
    #trim trajectory if bigger than N
    if returnval['Nc']>n: #trim coordinate file
        tempcrd=junkName()
        ptraj(unstopfile,crdfile).trajin(1,N).trajout(tempcrd)
        os.system('/bin/mv '+tempcrd+' '+crdfile) #overwrite old file
    if returnval['Nv']>n: #trim velocity file
        tempvel=junkName()
        ptraj(unstopfile,velfile).trajin(1,N).trajout(tempvel)
        os.system('/bin/mv '+tempvel+' '+velfile) #overwrite old file
    #strip water and ions before storing
    for (simple,ample) in ( (crd,crdfile), (vel,velfile)):
        ptraj(topfile,ample).strip(entities=striplist)['self'].go(outcrd=simple)
    check_rms(crd,rms) #create rms file
    #compress files
    for file in (crd,vel): os.system('gzip -f '+file)
    cmd='/bin/cp '+crdfile+' '+debugd+'/'+basename(crd)+' && '+\
         '/bin/cp '+velfile+' '+debugd+'/'+basename(vel)+' && '+\
         '/bin/cp '+outfile+' '+debugd+'/'+basename(out)
    os.system(cmd)
    cmd='/bin/rm '+crdfile+' && '+\
         '/bin/rm '+velfile+' && '+\
         '/bin/mv '+outfile+' '+out+' && '+\
         '/bin/cp '+rstfilenew+' '+rst+' && '+\
         '/bin/mv '+rstfilenew+' '+rstfile
    os.system(cmd)
    for file in (out,rms,rst): os.system('gzip -f '+file)
    os.system('echo '+basename(rst)+' >> '+rstlistf) #update rstlistf
    maxfr=int(open(maxfrfile,'r').readline()) #;Bye(maxfr)
    if M>=maxfr: Bye('reached maximum frame\n')
            
sys.exit(0)
