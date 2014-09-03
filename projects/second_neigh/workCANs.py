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
bch='/gpfs1/active/common/benchmark'
bchl=['bench100','bench200','bench300']
pcd=joink(codedir,'projects','second_neigh')
scd=joink(scratchdir,'second_neigh','out')
cad=joink(scratchdir,'second_neigh','pdb')
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


    def addOrigInfo(self):

        """add some info on original benchmark"""

        self.original={'tasser':{'infile':'','outfile':''}}
        for subd in bchl:
            inf=joink(bch,subd,'inputs',self.id+'.tasinp.bz')
            if isthere(inf):
                self.original['tasser']['infile']=inf
        for subd in bchl:
            outf=joink(bch,subd,'tasserout',self.id+'.tasout.bz')
            if isthere(outf):
                self.original['tasser']['outfile']=outf
        return True


    def addL(self):

        """find out length of target"""

        from utilities.small_utilities import unTARme
        if 'L' not in self.__dict__ and self.original['tasser']['infile']:
            wd=unTARme(self.original['tasser']['infile'],fileL=['rmsinp',],taropts='jxf')
            pr=open(joink(wd,'rmsinp'),'r')
            if pr: self.L=int(pr.readline().split()[1])
            else: sys.stderr.write('ERROR no rmsinp file for '+self.id)
            pastry('/bin/rm -r '+wd)
        return True


    def addemh(self):

        """is it easy, medium, or hard?"""

        from utilities.small_utilities import unTARme

        wd=unTARme(self.original['tasser']['infile'],fileL=['chain.dat',],taropts='jxf')
        if wd:
            self.original['tasser']['emh']=open(joink(wd,'chain.dat'),'r').readline().split()[1]
            os.system('/bin/rm -r '+wd)
        return True
    

    def submitSpk(self,joblist='',switch='',ibatch=0):
        
        """submit SPICKER simulations"""

        from jobs.job import job,genJobList

        if switch=='original':
            inTf=self.original['tasser']['infile']  ;  ouTf=self.original['tasser']['outfile']
            ouSf='out.spicker.tbz2'  ;  ibatch=0
            if inTf and ouTf:
                if not joblist: joblist=genJobList() #;Bye(joblist)
                jobname='spk.'+self.id  ;  outd=joink(self.simdir,'original')
                Job=job(name='spicker',exe='spicker',exed=joink(codedir,'f77/spicker/spicker0.0'),
                        inputs=[inTf,ouTf],
                        incmd='tar jxf %s && tar jxf %s && ls -1 rep*.tra|wc -l > tra.in &&  ls -1 rep*.tra >> tra.in'%( basename(inTf),basename(ouTf) ),
                        outcmd='tar jcf %s combo* closc* rst.dat && /bin/mv %s %s'%(ouSf,ouSf,outd)
                        )
                Job.qsub(jobname,outd,wallt='0.249',mem_limit='990',ddisk='400',joblist=joblist)
                ibatch=Job.dormant(ibatch)
            else:
                sys.stderr.write('ERROR (%s) no %s or %s\n'\
                                 %(self,id,basename(inTf),basename(ouTf))
                                 )
            return False,ibatch

        
    def submitTasSpk(self,joblist='',switch='',ibatch=0):

        """submit TASSER & SPICKER & amber minimization simulations"""

        from jobs.job import job,genJobList
        from tasser.tasser_utilities import calc_wall_time

        inf=self.original['tasser']['infile']
        ibatch=0
        if inf:
            if not joblist: joblist=genJobList() #;Bye(joblist)
            wt=calc_wall_time(L=self.L) #estimated TASSER wall time, in days
            wt='%5.2lf'%( float(wt)+0.50) #add 1/2 day for spicker plus amber minimizations
            jobname='x_'+self.id
            outd=self.simdir
            Job=job(name='tasser',exe='submitTasSpk.sh',exed=pcd,
                    args=' '+self.id+' '+inf+' '+self.simdir)
            Job.qsub(jobname,outd,wallt=wt,mem_limit='990',ddisk='400',joblist=joblist)
            ibatch=Job.dormant(ibatch)
            self.modified={'tasser':{'outfile':''},
                           'spicker':{'outfile':''},
                           'structures':{'outfile':''}
                           }
        return True


    def submitSpkStr(self,joblist='',switch='',ibatch=0):

        """submit spicker and amber minimization simulations of combo structures"""

        from jobs.job import job,genJobList
        
        if not joblist: joblist=genJobList() #;Bye(joblist)
        wt=0.249*float(self.L)/301 #300 residue protein takes six hours
        if switch=='original':
            intas=self.original['tasser']['infile']
            outtas=self.original['tasser']['outfile']
            if intas and outtas:
                outd=joink(self.simdir,'original')
                if not isthere(outd): pastry('/bin/mkdir -p '+outd)
                jobname='stro'+self.id
                self.original['spicker']={'outfile':''}
                self.original['structures']={'simdir':outd,'outfile':''}
                Job=job(name='structures',exe='submitSpkStr.sh',exed=pcd,\
                    args=' '+self.id+' '+outtas+' '+outd+' '+intas)
                Job.qsub(jobname,outd,wallt=wt,mem_limit='250',ddisk='300',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        elif switch=='modified':
            inf=joink(self.simdir,'out.spicker.tbz2')
            if isthere(inf):
                outd=self.simdir
                jobname='strm'+self.id
                Job=job(name='structuresm',exe='submitSpkStr.sh',exed=pcd,\
                    args=' '+self.id+' '+inf+' '+outd)
                Job.qsub(jobname,outd,wallt=wt,mem_limit='250',ddisk='300',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        return ibatch,True


    def submitStr(self,joblist='',switch='',ibatch=0):

        """submit spicker and amber minimization simulations of combo structures"""

        from jobs.job import job,genJobList
        
        if not joblist: joblist=genJobList() #;Bye(joblist)
        wt=0.249*float(self.L)/301 #300 residue protein takes six hours
        if switch=='original':
            pt=self.original['spicker']
            outspk=joink(self.simdir,'original','out.spicker.tbz2')
            if isthere(outspk):
                pt['outfile']=outspk
                outd=joink(self.simdir,'original')
                jobname='str'+self.id
                self.original['structures']={'simdir':outd,'outfile':''}
                Job=job(name='structures',exe='submitStr.sh',exed=pcd,\
                    args=' '+self.id+' '+outspk+' '+outd)
                Job.qsub(jobname,outd,wallt=wt,mem_limit='250',ddisk='300',joblist=joblist)
                ibatch=Job.dormant(ibatch)
                return True,ibatch
        return False,ibatch


    def createTasserSpickerStructuresDicts(self,switch=''):

        """initialize 'tasser', 'spicker', and 'structures' dictionaries"""

        from utilities.small_utilities import unTARme,chomp

        if switch=='modified':
            ptb=self.modified
            ptb['simdir']=self.simdir
            for dict in ('tasser','spicker','structures'):
                if dict not in ptb.keys(): ptb[dict]={}
                if 'outfile' not in ptb[dict].keys(): ptb[dict]['outfile']=''
                outfile=joink(ptb['simdir'],'out.'+dict+'.tbz2')
                if isthere(outfile): ptb[dict]['outfile']=outfile
            if ptb['structures']['outfile']:
                ptc=ptb['structures']
                currd=os.getcwd()
                wd=unTARme( ptc['outfile'] ) #untar structures to temporary directory
                os.chdir(wd)
                ptc['list']=chomp(os.popen('ls -1 *.pdb').readlines())
                os.chdir(currd)
                os.system('/bin/rm -r '+wd)
                for basename in ptc['list']: ptc[basename]={}
            return True
        elif switch=='original':
            ptb=self.original
            ptb['simdir']=joink(self.simdir,'original')
            for dict in ('spicker','structures'):
                if dict not in ptb.keys(): ptb[dict]={}
                if 'outfile' not in ptb[dict].keys(): ptb[dict]['outfile']=''
                outfile=joink(ptb['simdir'],'out.'+dict+'.tbz2')
                if isthere(outfile): ptb[dict]['outfile']=outfile
            if ptb['structures']['outfile']:
                ptc=ptb['structures']
                currd=os.getcwd()
                wd=unTARme( ptc['outfile'] ) #untar structures to temporary directory
                os.chdir(wd)
                ptc['list']=chomp(os.popen('ls -1 *.pdb').readlines())
                if not ptc['list']: sys.stderr.write('ERROR '+self.id+'\n')
                os.chdir(currd)
                os.system('/bin/rm -r '+wd)
                for basename in ptc['list']: ptc[basename]={}
            return True
    

    def gatherStructures(self,switch=''):

        """gather structures from the TASSER simulations"""

        from utilities.small_utilities import unTARme,chomp

        strd=joink(scd,'structures_CA')
        if not isthere(strd): os.system('/bin/mkdir -p '+strd)
        p0=self.modified['structures']
        if not isthere(p0['outfile']): return False #no minimized structures
        wd=unTARme( p0['outfile'] ) #untar structures to temporary directory
        for strf in p0['list']: #basename of structure files
            f0=joink(wd,strf)   #full atom structure file, origin
            f1=joink(strd,strf) #CA-only structure file, destination
            cmd='grep " CA " '+f0+' > '+f1 #;Bye(cmd)
            pastry(cmd)
        pastry('/bin/rm -r '+wd)
            

    def add_natf(self):

        """add path to native file"""

        natpdb=joink(cad,self.id+'.pdb')
        if isthere(natpdb):
            self.natf=natpdb
            return True
        natf=''
        for pdbd in pdbds:
            x=joink(pdbd,self.id+'.pdb')
            if isthere(x):
                natf=x
                break
        if natf:
            os.system('/bin/cp %s %s'%(self.natf,cad))
            self.natf=natpdb
        else:
            if os.system('create_single_chain.py -i %s -f %s.pdb -o %s'%(self.id,self.id,cad)):
                sys.stderr.write('ERROR no native file for %s\n'%(self.id))
            else: self.natf=natpdb
        return True
    
            
    def AddSpkobj(self,switch='original'):

        """add spicker object attribute"""

        from utilities.small_utilities import unTARme
        from spicker.spickerYangResultsManager import spickOut

        pt=self.original['spicker']
        if switch=='modified': pt=self.modified['spicker']

        if not pt['outfile']: #check if spicker output exists
            spkout=joink(self.simdir,'original','out.spicker.tbz2')
            if switch=='modified': spkout=joink(self.simdir,'out.spicker.tbz2')
            if isthere(spkout): pt['outfile']=spkout
        if not pt['outfile']:
            sys.stderr.write('ERROR no spicker output for %s %s\n'%(switch,self.id))
            return False

        if 'spkobj' not in pt.keys(): pt['spkobj']=None
#        if pt['spkobj']:
#            sys.stdout.write(self.id+' already processed\n')
#            return False #We already processed this one
        wd=joink(self.simdir,'original','spicker')
        wd=unTARme(pt['outfile'],wd=wd) #untar output file to simdir directory
        if wd and self.natf:
            pt['spkd']=wd
            spkobj=spickOut(dir=wd,nat=self.natf)
            if spkobj.readError:
                sys.stderr.write('ERROR '+spkobj.readError+'\n')
                return False
            densestfive=spkobj.rankIDsByDens()[0:5]
            spkobj.initTMtoNat()
            idBestTM=spkobj.rankIDsByTMtoNat()[0]
            spkobj.initRmsdsToNat()
            pt['spkobj']=spkobj
            id=spkobj.rankIDsByRmsdtoNat(list=densestfive)[0]
            pt['bestByRMSD']={'id':id,'rmsd':spkobj.rmsds[id]}
            id=spkobj.rankIDsByTMtoNat(list=densestfive)[0]
            pt['bestByTM']={'id':id, 'tm':spkobj.tms[id]}
            id=spkobj.densest
            pt['bestByDens']={'id':id,'rmsd':spkobj.rmsds[id],'tm':spkobj.tms[id]}
        if wd: os.system('/bin/rm -r '+wd)
        return True


    def outputCentroidMarks(self):

        """output rmsd,TM,rank for some centroids"""

        buf=self.id
        for x in ('original','modified'):
            if not self.__dict__[x]['spicker']['outfile']: return False
            a=self.__dict__[x]['spicker']['bestByDens']
            b=self.__dict__[x]['spicker']['bestByRMSD']
            c=self.__dict__[x]['spicker']['bestByTM']
            buf+= ' %2d %4.1lf %4.2lf %2d %4.1lf %2d %4.2lf %s'\
                 %(a['id'],a['rmsd'],a['tm'],b['id'],b['rmsd'],c['id'],c['tm'],
                   self.original['tasser']['emh'])
        sys.stdout.write(buf+'\n')
        return False
    
        
########################################################################

    
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    pdbid ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'natf': path to native file name
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'simdir': directory to store simulation results
    |_'original': contains info on original benchmark
    |  |_'tasser': contains info on tasser simulations
    |  |  |_'infile': input for tasser
    |  |  |_'emh': easy, medm, hard
    |  |  |_'outfile': output for tasser
    |  |_'spicker':
    |  |  |_'outfile':
    |  |  |_'spkd': directory name where to uncompress outfile
    |  |  |_'spkobj': spicker object
    |  |  |_'bestByDens': best centroid by density
    |  |     |_'id': its id
    |  |     |_'rmsd': it's rmsd to native
    |  |     |_'tm': its TM-score to native
    |  |  |_'bestByRMSD': best centroid by RMSD to native in the top five (ranked by density)
    |  |     |_'id'
    |  |     |_'rmsd'
    |  |  |_'bestByTM':  best centroid by TM to native in the top five (ranked by density)
    |  |     |_'id'
    |  |     |_'tm'
    |  |_'structures': info on amber minimized structures
    |    |_'outfile':
    |_'modified': contains info on simulations with second neighbor interaction
      |_'tasser':
      | |_'outfile'
      |_'spicker':
      | |_'outfile'
      |_'structures': info on amber minimized structures
        |_'outfile'
    """
    print globals()['datatree'].__doc__


def doitall(id):

    """a recount of everything we had done on each varset

    Useful if we screw-up a particular dump file, because we can redo
    it all"""

    c=PEPSI(id=id)
    c.genDumpf(cl.repod)

    methods=('createSimdir','addL','addOrigInfo','submitTasSpk')
    for method in methods:
        print method
        getattr(objSL,method)()
    #finally secure the beast!
    c.pickleDump()


def createList():

    """create list of pdb headers"""

    os.popen('ls -1 '+ocd+' | grep -e "\.list"| cut -d "." -f 1 > '+canlistf)

              
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


def addOrigInfo(id=''):

    """add some info on original benchmark"""

    if id:
        c=cl.loadCAN(id)
        if c.addOrigInfo()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addOrigInfo',silent=False)

def addL(id=''):

    """find out length of target"""
    
    if id:
        c=cl.loadCAN(id)
        if c.addL()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addL',silent=False)


def submitTasSpk(id='',switch=''):

    """submit TASSER & SPICKER & amber minimization simulations"""

    from jobs.job import genJobList

    joblist=genJobList() #;Bye(joblist)
    if id:
        c=cl.loadCAN(id)
        if c.submitTasSpk(joblist=joblist,switch=switch)!=False: c.pickleDump()
    else:
        ibatch=0 ; remaining=len(canlist)
        for c in cl.iterate():
            print remaining,c.id
            ibatch,record=c.submitTasSpk(joblist=joblist,switch=switch,ibatch=ibatch)
            if record!=False: c.pickleDump()
            remaining-=1


def gatherStructures(id=''):
    
    """gather structures from the TASSER simulations"""
    
    if id:
        c=cl.loadCAN(id)
        if c.gatherStructures()!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            print remaining,c.id
            if c.gatherStructures()!=False: c.pickleDump()
            remaining-=1

def view(id):

    """view object properties"""

    cl.loadCAN(id).view()


def submitSpk(id='',switch='',idlist=''):
        
    """submit SPICKER simulations"""

    from jobs.job import genJobList
    joblist=genJobList()
    if idlist:
        idlist=chomp(open(idlist,'r').readlines())
    else:
        if id: idlist=[id,]
        else: idlist=cl.canlist
    ibatch=0 ; remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        record,ibatch=c.submitSpk(joblist=joblist,switch=switch,ibatch=ibatch)
        if record!=False: c.pickleDump()
        remaining-=1


def submitStr(id='',switch='',idlist=''):

    """submit amber minimization simulations of combo structures"""

    from jobs.job import genJobList
    joblist=genJobList()

    if idlist:
        idlist=chomp(open(idlist,'r').readlines())
    else:
        if id: idlist=[id,]
        else: idlist=cl.canlist
        
    ibatch=0 ; remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        record,ibatch=c.submitStr(joblist=joblist,switch=switch,ibatch=ibatch)
        if record!=False: c.pickleDump()
        remaining-=1


def submitSpkStr(id='',switch='',idlist=''):

    """submit amber minimization simulations of combo structures"""

    from jobs.job import genJobList
    joblist=genJobList()
    
    if idlist:
        idlist=chomp(open(idlist,'r').readlines())
    else:
        if id: idlist=[id,]
        else: idlist=cl.canlist
    ibatch=0 ; remaining=len(cl.canlist)
    for c in cl.iterate():
        print remaining,c.id
        ibatch,record=c.submitSpkStr(joblist=joblist,switch=switch,ibatch=ibatch)
        if record!=False: c.pickleDump()
        remaining-=1


def filterByTassInfile():

    """find out which targets do have TASSER input in the benchmark

    PROSPECTOR did not run for some targets in the benchmark. We
    should exclude these"""

    remaining=len(cl.canlist)
    for c in cl.iterate():
        if c.original['tasser']['infile']: print c.id

        
def createTasserSpickerStructuresDicts(id='',switch=''):

    """initialize 'tasser', 'spicker', and 'structures' dictionaries"""

    if id:
        c=cl.loadCAN(id)
        if c.createTasserSpickerStructuresDicts(switch=switch)!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            print remaining,c.id
            if c.createTasserSpickerStructuresDicts(switch=switch)!=False: c.pickleDump()
            remaining-=1


def addemh(id=''):

    """is it easy, medium, or hard?"""

    if id:
        c=cl.loadCAN(id)
        if c.addemh()!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            print remaining,c.id
            if c.addemh()!=False: c.pickleDump()
            remaining-=1
            

def add_natf(id=''):

    """add path to native file"""

    if id:
        c=cl.loadCAN(id)
        if c.add_natf()!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            print remaining,c.id
            if c.add_natf()!=False: c.pickleDump()
            remaining-=1
            

def AddSpkobj(id='',switch='original',idlist=[]):

    """add spicker object attribute"""

    if idlist:
        idlist=chomp(open(idlist,'r').readlines())
    else:
        if id: idlist=[id,]
        else: idlist=cl.canlist
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.AddSpkobj(switch=switch)!=False: c.pickleDump()
        remaining-=1
            

def checkSomething():

    for c in cl.iterate():
        print c.original['tasser']['infile']


def outputCentroidMarks(id=''):

    """output rmsd,TM,rank for bestByDens, bestByRMSD, and bestByTM centroids"""

    print '#(2)(3)(4):rank,rmsd,tm for top cluster (densest cluster) in original simulation'
    print '#(5)(6):rank,rmsd for best rmsd cluster in original simulation'
    print '#(7)(8):rank,tm for best TM cluster in original simulation'
    print '#(9)(10)(11):rank,rmsd,tm for top cluster (densest cluster) in modified simulation'
    print '#(12)(13):rank,rmsd for best rmsd cluster in modified simulation'
    print '#(14)(15):rank,tm for best TM cluster in modified simulation'
    print '#(16): easy/memd/hard'
    print '# id  (2) (3) (4)  (5) (6) (7) (8) (9)(10) (11)(12)(13)(14)(15) (16)'
    if id:
        c=cl.loadCAN(id)
        if c.outputCentroidMarks()!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            #print remaining,c.id
            if c.outputCentroidMarks()!=False: c.pickleDump()
            remaining-=1


def help(*kargs):

    """list available services"""
    
    services=['createList','datatree','doitall','updateCollection','initCollection',
              'createSimdir','addL','addOrigInfo','submitTasSpk','view'
              'gatherStructures','submitSpkStr','filterByTassInfile',
              'createTasserSpickerStructuresDicts','submitSpk','addemh',
              'natf','AddSpkobj','submitStr','outputCentroidMarks']
    
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


