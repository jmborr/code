#!/usr/bin/python

import os,sys,re
from utilities.can import CAN,CANl
from utilities.codedir import nodecodedir,codedir,scratchdir
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
seqL=150
joink=os.path.join
isthere=os.path.exists

pcd=joink(codedir,'projects','Hellinga')
scd=joink(scratchdir,'Hellinga/buried_Cys','out')
collf=joink(scd,'varsets.can')
canlistf=joink(scd,'list')
repod=joink(scd,'repod')

class PEPSIl(CANl):

    """inherited class from CANl, so that we can add methods to fuzz
    around the attributes of our particular CANl object
    implementation"""

    pass

#list of CAN objects for my project
cl=PEPSIl(collf,canlistf,repod)
canlist=chomp(open(cl.canlistf,'r').readlines())

#####################################################
####  OBJECT IMPLEMENTATION FOR MY PROJECT       ####
#####################################################
class PEPSI(CAN):
    
    """inherited class from CAN, so that we can add methods to fuzz
    around the attributes of our particular CAN object
    implementation"""


    def createSimdir(self):
        """create subdirectory where to put simulations output"""
        outd=joink(scd,self.id)
        os.system('/bin/mkdir -p '+outd)
        self.simdir=outd
        return True


    def addCA(self):
        """add absolute file name of target structure"""
        self.CA=joink(self.simdir,self.id+'.pdb')
        print self.CA
        return True #True will pickledump


    def createBatch(self,batch='',batchfile=''):
        """create directory structure for a batch of sequences, as well as fasta files"""
        from time import sleep
        if not 'batch' in self.__dict__.keys(): self.batch={}
        self.batch[batch]={'simdir':joink(self.simdir,batch),
                           'batchfile':batchfile}
        pt=self.batch[batch]  ;  os.system('/bin/mkdir -p '+pt['simdir'])
        bf=open(batchfile,'r')  ;  bf.readline()  ;  bf.readline()  ;  bf.readline()
        N=int(bf.readline())  ;  pt['N']=N  ;  buf=[]
        for i in range(0,N):
            name='s%04d'%(i)  ;  buf.append(name)  ;  dir=joink(pt['simdir'],name)
            print self.id,batch,name  ;  sleep(0.1)  ;  os.system('/bin/mkdir -p '+dir)
            ff=open(joink(dir,name+'.fasta'),'w')
            ff.write('>'+name+'\n'+bf.readline()+bf.readline()+bf.readline())
            ff=open(joink(dir,name+'.hongpred'),'w')
            ff.write('>'+name+'\n'+bf.readline()+bf.readline()+bf.readline())
        pt['listfile']=joink(pt['simdir'],'list')
        open(pt['listfile'],'w').write( '\n'.join(buf) )
        return True


    def storeCYSposition(self):
        """store the position of the CYS"""
        positions={'a0512':{'1':111,'2':131, '3':71},
                   'a0653':{ '1':17,'2':142, '3':71}}
        for batchid in ('1','2','3'):
            batch='seqlist_design_nnsecsel3_'+batchid
            self.batch[batch]['CYSpos']=positions[self.id][batchid]
        return True
    
        
    def sendBlast(self,batch='',task=''):
        """send groups of blast jobs"""
        pt=self.batch[batch]
        if task=='check jobs':
            for i in range(0,pt['N']):
                name='s%04d'%(i)
                logf=joink(pt['simdir'],name,'%s.%s%s.log'%(self.id,batch[-1:],name))
                if not isthere(logf):
                    sys.stderr.write('no '+logf+'\n')
                else:
                    if open(logf,'r').readlines()[-1].find('mode = 1')>=0:
                        sys.stderr.write('error '+logf+'\n')
            return False
        prefix='-n '+self.id+'.'+batch[-1:]
        groupJobs='-p 20'    ;  extraFlags=''         ; numDays='-d 0.249'  ;  memReserv='-k 900'
        listf='-a '+pt['listfile'] ;
        outdir=joink(pt['simdir'],'xxxxx')   ; fastaf=joink(outdir,'xxxxx.fasta')
        exeline='%s/python/combo_jobs/predictor/blasts.py -a %s -b %s'%\
                 (nodecodedir,fastaf,outdir)
        fake='' #'-o yes'
        cmd='%s/python/combo_jobs/generic_job_list.py %s %s %s %s %s -c %s %s -g "%s" %s'%\
             (codedir,prefix,groupJobs,extraFlags,numDays,memReserv,outdir,listf,exeline,fake)
        os.system(cmd)
        return True


    def createSeqDat(self,batch=''):
        """create seq.dat file from hongpred file"""
        from secondary.ss2seqdat import genSeqdat
        pt=self.batch[batch]
        for i in range(0,pt['N']):
            name='s%04d'%(i)  ;
            seqf=joink(pt['simdir'],name,name+'.fasta')
            ssf=joink(pt['simdir'],name,name+'.hongpred')
            seqdatf=joink(pt['simdir'],name,'seq.dat')
            open(seqdatf,'w').write( genSeqdat(seqf,ssf,lett2numb='L 1 H 2 B 4') )

        
    def ab_pair_pro_dat_tas_skp_duk(self,batch='',task='',joblist=[]):
        """do ab initio tasser plus input"""
        from jobs.job import job,genJobList
        from tasser.tasser_utilities import calc_wall_time
        if not joblist: joblist=genJobList()
        wt='%5.2lf'%( float(calc_wall_time(L=seqL))+0.50)  ;  ibatch=0  ;  pt=self.batch[batch]
        inddf=joink(pcd,'brun','in.dd')
        for i in range(0,pt['N']):
            name='s%04d'%(i)  ;  outdir=joink(pt['simdir'],name)
            fastaf=joink(outdir,name+'.fasta') ; blastout=joink(outdir,name+'.out.blast.tar.bz2')
            expdatf=joink(outdir,'exp.dat')    ; seqdatf=joink(outdir,'seq.dat')
            args='-a %s -b %s -c %s -g %s -j %s -k %s -d %s'%\
                  (fastaf,blastout,outdir,expdatf,inddf,seqdatf,name)
            Job=job(name='ab_pair_pro_dat_tas_skp_duk',exe='ab_pair_pro_dat_tas_skp_duk.py',
                    exed=pcd,args=args)
            jobname='a'+self.id+batch[-1]+name
            Job.qsub(jobname,outdir,wallt=wt,mem_limit='990',ddisk='400',joblist=joblist)
            ibatch=Job.dormant(ibatch)
        return False


    def tmscore(self,batch='',task='print'):
        """calculate TMscore of closc to native"""
        from utilities.small_utilities import unTARme        
        pt=self.batch[batch]
        if task=='calculate':
            pt['neighbors']={}
            tmscore='/library/jose/code/python/tmScore/tm_score.py'
            tmco=0.4
            for header in chomp(open(pt['listfile'],'r').readlines()):
                spkout=joink(pt['simdir'],header,header+'.out.ab_spicker.tar.bz2')
                if isthere(spkout):
                    wd=unTARme(spkout)
                    for i in range(1,11):
                        closc=joink(wd,'closc'+`i`+'.pdb')
                        tm=float(os.popen(tmscore+' -a '+closc+' -b '+self.CA).readline().strip())
                        if tm>tmco: pt['neighbors'][header+'.%02d'%(i)]=tm
                    os.system('/bin/rm -r '+wd)
            items=[(v, k) for k, v in pt['neighbors'].items()]
            items.sort()
            items.reverse()
            pt['neighbors']=items
            return True
        elif task=='calculate all':
            scores={} ; tmscore='/library/jose/code/python/tmScore/tm_score.py'
            for header in chomp(open(pt['listfile'],'r').readlines()):
                print self.id,batch,header
                spkout=joink(pt['simdir'],header,header+'.out.ab_spicker.tar.bz2')
                if isthere(spkout):
                    wd=unTARme(spkout)
                    for i in range(1,11):
                        closc=joink(wd,'closc'+`i`+'.pdb')
                        tm=float(os.popen(tmscore+' -a '+closc+' -b '+self.CA).readline().strip())
                        scores[header+'.%02d'%(i)]=tm
                    os.system('/bin/rm -r '+wd)
            items=[(v, k) for k, v in scores.items()] ; items.sort() ; items.reverse()
            buf=''
            for (tm,modl) in items:  buf+='%s %5.3lf\n'%(modl,tm)
            open(joink(scd,'analysis',self.id+'.'+batch+'.alltm.dat'),'w').write(buf)
            return False            
        elif task=='print':
            buf=''
            for (tm,modl) in self.batch[batch]['neighbors']: buf+='%s %5.3lf\n'%(modl,tm)
            open(joink(scd,'analysis',self.id+'.'+batch+'.tm.dat'),'w').write(buf)
            return False        
        

    def tmscoreHist(self):
        """calculate histograms for all tm scores stored"""
        for batch in self.batch.keys():
            hist=[0.0]*100 ; tmfile=joink(scd,'analysis',self.id+'.'+batch+'.alltm.dat')  ;  N=0
            for line in open(tmfile,'r').readlines():
                (modl,tm)=line.split()  ;  tm=int(float(tm)*100)  ;  hist[tm]+=1  ;  N+=1
                if tm>100:sys.stderr.write('error\n')
            tmfile=joink(scd,'analysis',self.id+'.'+batch+'.alltm.hist.dat')  ;  buf=''
            for tm in range(100): buf+='%4.2lf %6.4lf\n'%(float(tm)/100,float(hist[tm])/N)
            open(tmfile,'w').write(buf)
            
        
    def preparePackages(self):
        """prepare packages to send to Prof. Hellinga"""
        from utilities.small_utilities import unTARme
        os.chdir(joink(scd,'analysis'))
        dir=joink(scd,'analysis',self.id)  ;  os.system('/bin/mkdir -p '+dir)
        os.system( '/bin/cp %s %s'%(self.CA,dir) )
        for batch in self.batch.keys():
            pt=self.batch[batch]
            dir_b=joink(dir,'CYS_at_res'+`pt['CYSpos']`)  ;  os.system('/bin/mkdir -p '+dir_b)
            L=len(pt['neighbors'])  ;  N=100  ;  minN=min(L,N)  ;  n=0
            for (tm,modl) in pt['neighbors']:
                seq=modl[0:5]  ;  closcID=modl[-2:]
                dukka=joink(pt['simdir'],seq,seq+'abdukka.tar.bz2')
                if isthere(dukka):
                    wd=unTARme(dukka)
                    cmd='/bin/cp %s %s'%(joink(wd,'closc'+closcID+'.pdb.rebuilt' ),
                                         joink(dir_b, seq+'.pdb'))
                    os.system(cmd)  ;  os.system('/bin/rm -r '+wd)
                    n+=1
                if n>N: break
            os.system('/bin/tar jcf %s %s'%(self.id+'.tbz2',self.id))
        return False
                
        
#################### end of PEPSICAN  ##################
    
    
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    ojbect ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'simdir': directory to store simulation results
    |_'CA': native file
    |_'batch': dictionary containing results for different batches of sequences
      |_batch: dictionary containing results for one batch of sequences
        |_'N': number of sequences
        |_'CYSpos': position of CYS
        |_'listfile': file containing list of sequences
        |_'batchfile': file provided by Jeff containing sequences
        |_'simdir': root directory for simulations of the batch sequences
    """
    print globals()['datatree'].__doc__


def initCollection():

    """Initialize dump pickle files and collection file

    initCollection(id='')
    Initialize set of pickle files and collection file to store info
    on the splice variants. Remove previously stored pickle and
    collection file

    id: initialize only one dump file instead of all dump files and
    the collection file

    ***************************************************
    * WARNING!: this service will remove any existing *
    *           dump files and the collection file.   *
    ***************************************************"""

    print '**************************************************************'
    print '* WARNING!: this service will remove any existing dump files *'
    print '*         and the collection file. Do you want to proceed?   *'
    print '**************************************************************'
    sys.stdout.write('(y/n): ')
    reply=sys.stdin.readline().strip()
    if reply != 'y': return

    if not isthere(repod): os.system('/bin/mkdir -p '+repod)
    remaining=len(canlist)
    for header in canlist:
        c=PEPSI(id=header)
        c.genDumpf(cl.repod)
        c.pickleDump()
        print remaining,header
        remaining-=1


def updateCollection():

    """update the collection file"""
    
    cl.updColletion()
    

def createSimdir(id=''):

    """create subdirectory where to put simulations output"""

    if id:
        c=cl.loadCAN(id)
        if c.createSimdir()!=False: c.pickleDump()
    else:
        cl.cycleOverList('createSimdir',silent=True)



def view(id):

    """viewver of object properties"""

    cl.loadCAN(id).view()
    

def addCA(*kargs,**kwargs):
    """add absolute file name of target structure"""
    cl.cycleOverList('addCA',*kargs,**kwargs)


def createBatch(*kargs,**kwargs):
    """create directory structure for a batch of sequences, as well as fasta files

    Example: ... -c id=a0512,batch=the_batchname,batchfile=the_batchfile"""
    cl.cycleOverList('createBatch',*kargs,**kwargs)


def storeCYSposition(*kargs,**kwargs):
    """store the position of the CYS
    Example: ... -c id=a0512"""
    cl.cycleOverList('storeCYSposition',*kargs,**kwargs)
        

def sendBlast(*kargs,**kwargs):
    """send groups of blast jobs
    Example: ... -c id=a0512,batch=the_batchname"""
    cl.cycleOverList('sendBlast',*kargs,**kwargs)


def createSeqDat(*kargs,**kwargs):
    """create seq.dat files
    Example: ... -c id=a0512,batch=the_batchname"""
    cl.cycleOverList('createSeqDat',*kargs,**kwargs)


def ab_pair_pro_dat_tas_skp_duk(*kargs,**kwargs):
    """do ab initio tasser plus input
    Example: ... -c id=a0512,batch=the_batchname"""
    from jobs.job import genJobList    
    kwargs['joblist']=genJobList()
    cl.cycleOverList('ab_pair_pro_dat_tas_skp_duk',*kargs,**kwargs)


def number_seqs(*kargs,**kwargs):
    """calculate total number of sequences"""    
    N=0
    for id in cl.canlist:
        c=cl.loadCAN(id)
        for batch in c.batch.keys(): N+=c.batch[batch]['N']
    print N

    
def tmscore(id='',batch='',task=''):
    """calculate TMscore of closc to native
     Example: ... -c id=a0512,batch=the_batchname"""
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    for id in idlist:
        print remaining,id
        c=cl.loadCAN(id)
        if c.tmscore(batch=batch,task=task)!=False: c.pickleDump()
        remaining-=1
    

def tmscoreHist(*kargs,**kwargs):
    """calculate histograms for all tm scores stored"""
    cl.cycleOverList('tmscoreHist',*kargs,**kwargs)

    
def preparePackages(*kargs,**kwargs):
    """prepare packages to send to Prof. Hellinga
    Example: ... -c id=a0512"""
    cl.cycleOverList('preparePackages',*kargs,**kwargs)
     

def help(*kargs):

    """list available services"""
    
    services=['datatree','updateCollection','initCollection','createSimdir',
              'view']
    os.system('/usr/bin/clear')
    if not kargs:
        print 'Services='
        for service in services:
            print '  '+service+' : '+globals()[service].__doc__.split('\n')[0]
        print 'Type "workCANs.py -a help -b service" for more detail on "service"'
    else:
        sv=kargs[0]
        print globals()[sv].__doc__        


if __name__=='__main__':
    inpHand('Usage: workCANs.py [options]',
            ' -a _R_service required service. Type "help" for list of services.',              
            ' -b __servargs service required arguments in python syntax, and separate them with ",". Escape symbols \' and " (def:None)',
            ' -c __servoptargs service optional arguments (Example: "key=val" for one argument and "key1=val1, key2=val2,..." for two or more arguments (def:None)',
            ).parse(locals(),sys.argv)

    #parse required arguments [arg1,arg2,...]
    reqargs=[]
    if servargs:
        reqargs=servargs.split(',')
    #parse optional arguments {key1:val1, key2:val2,...}
    optsargs={}
    if servoptargs:
        for pair in servoptargs.split(','):
            key,val=pair.split('=')
            optsargs[key]=val
    #execute the service
    locals()[service](*reqargs,**optsargs)
