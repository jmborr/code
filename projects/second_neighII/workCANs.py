#!/usr/bin/python

import os,sys,re
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import inpHand
from jobs.job import job,pastry

"""manage poly-VAL peptides"""

#global variables
joink=os.path.join
isthere=os.path.exists
basename=os.path.basename
pcd=joink(codedir,'projects','second_neighII')
scd=joink(scratchdir,'second_neighII','out')
seqds=('/gpfs1/scratch/adrian/dat/pdb_jul07/seq','/gpfs1/scratch/adrian/dat/pdb_mar07/seq',
       '/gpfs1/scratch/adrian/dat/pdb200/seq')
pdbds=('/gpfs1/scratch/adrian/dat/pdb_jul07/CA','/gpfs1/scratch/adrian/dat/pdb_mar07/CA',
       '/gpfs1/scratch/adrian/dat/pdb200/CA') #different pdb releases
collf=joink(scd,'varsets.can')
canlistf=joink(scd,'list')
repod=joink(scd,'repod')

class PEPSIl(CANl):

    """inherited class from CANl, so that we can add methods to fuzz
    around the attributes of our particular CANl object
    implementation"""

    pass

#CAN objects for the simil library of PDB structures
cl=PEPSIl(collf,canlistf,repod)
canlist=chomp(open(cl.canlistf,'r').readlines())


#####################################################
####  OBJECT FOR SIMIL LIBRARY OF PDB STRUCTURES ####
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


    def importTASSERinput(self):

        """import TASSER imput"""

        inpdir=joink(scratchdir,'createInputsTasser',self.id[1],self.id)
        os.chdir( inpdir )
        if os.system('cp pair1.dat pair3.dat '+self.simdir): return False
        os.chdir( joink(inpdir,'abinitio') )
        if os.system('/bin/cp exp.dat CA '+self.simdir): return False
        for seqd in seqds:
            if not os.system('/bin/cp %s/%s.SEQ %s/seq.dat'%(seqd,self.id,self.simdir)): break
        os.chdir( self.simdir )
        self.L=int( os.popen('wc -l CA').readline().split()[0] )
        open('rmsinp','w').write('1 '+`self.L`+'\n'+`self.L`+'\n'+self.id)
        if os.system('tar jcf in.tasser.tbz2 pair1.dat pair3.dat exp.dat rmsinp seq.dat'):
            sys.stderr.write('ERROR %s\n'%(self.id))
            return False
        os.system('/bin/rm pair1.dat pair3.dat exp.dat rmsinp seq.dat')
        self.inpTas=joink(self.simdir,'in.tasser.tbz2')
        return True


    def submitTasSpiCa2(self,joblist='',switch='',ibatch=0):

        """submit TASSER, SPICKER, ca2py"""

        from jobs.job import job,genJobList
        from tasser.tasser_utilities import calc_wall_time
        inf=self.inpTas
        if not joblist: joblist=genJobList()
        wt='%5.2lf'%( float( calc_wall_time(L=self.L) )+0.50)
        jobname=switch[-5:]+'.'+self.id
        outd=joink(self.simdir,switch)
        Job=job(name='simulation',exe='submitTasSpkCa2.sh',exed=pcd,
                args='%s %s'%(self.id,switch))
        Job.qsub(jobname,outd,wallt=wt,mem_limit='990',ddisk='400',joblist=joblist)
        ibatch=Job.dormant(ibatch)
        return False,ibatch
        
    
    def AddSpkobj(self,switch=''):

        """add spicker object attribute"""

        from utilities.small_utilities import unTARme
        from spicker.spickerYangResultsManager import spickOut

        self.__dict__[switch]={
            'simdir':joink(self.simdir,switch),
            'tasser':{'outfile':''},
            'spicker':{'outfile':''},
            'structures':{'outfile':''},
            }
        pt=self.switch
        for task in ('tasser','spicker','structures'):
            taskoutfile=joink(pt['simdir'],'out'+task+'.tbz2')
            if isthere(taskoutfile): pt[task]['outfile']=taskoutfile
        pts=pt['spicker']
        if pts['outfile']:
            pts['spkobj']=None
            wd=joink(pt['simdir'],'spicker') #temporary directory where to unzip
            wd=unTARme(pts['outfile'],wd=wd)
            pts['spkd']=wd
            spkobj=spickOut(dir=wd,nat=self.natf)
            if spkobj.readError:
                sys.stderr.write('ERROR '+spkobj.readError+'\n')
                return False
            densestfive=spkobj.rankIDsByDens()[0:5]
            spkobj.initTMtoNat()
            idBestTM=spkobj.rankIDsByTMtoNat()[0]
            spkobj.initRmsdsToNat()
            pts['spkobj']=spkobj
            id=spkobj.rankIDsByRmsdtoNat(list=densestfive)[0]
            pts['bestByRMSD']={'id':id,'rmsd':spkobj.rmsds[id]}
            id=spkobj.rankIDsByTMtoNat(list=densestfive)[0]
            pts['bestByTM']={'id':id, 'tm':spkobj.tms[id]}
            id=spkobj.densest
            pts['bestByDens']={'id':id,'rmsd':spkobj.rmsds[id],'tm':spkobj.tms[id]}
        pts=pt['structures']
        if pts['outfile']:
            wd=joink(pt['simdir'],'spicker') #temporary directory where to unzip
            wd=unTARme(pts['outfile'],wd=wd)
            pts['strd']=wd
        if wd: os.system('/bin/rm -r '+wd)
        return True


    
########################################################################

    
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    pdbid ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'simdir': directory to store simulation results
    |_'inpTas': input file for TASSER, excluding in.dd file
    |_'L': sequence length
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

    if not isthere(cl.repod): pastry('/bin/mkdir -p '+cl.repod)
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

    """view object properties"""

    cl.loadCAN(id).view()


def importTASSERinput(id=''):
    
    """import TASSER imput"""

    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print c.id,remaining
        if c.importTASSERinput()!=False: c.pickleDump()
        remaining-=1


def submitTasSpiCa2(id='',switch=''):

    """submit TASSER, SPICKER, ca2py"""

    from jobs.job import genJobList

    joblist=genJobList() #;Bye(joblist)
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist) ; ibatch=0
    for id in idlist:
        c=cl.loadCAN(id)
        print c.id,remaining
        record,ibath=c.submitTasSpiCa2(joblist=joblist,switch=switch,ibatch=ibatch)
        if record!=False: c.pickleDump()
        remaining-=1
    

def help(*kargs):

    """list available services"""
    
    services=['createList','datatree','doitall','updateCollection','initCollection',
              'createSimdir','view','importTASSERinput']
    
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


