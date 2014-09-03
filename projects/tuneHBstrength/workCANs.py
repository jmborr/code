#!/usr/bin/python

import os,sys,re
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
joink=os.path.join
isthere=os.path.exists

pcd=joink(codedir,'projects','tuneHBstrength')
scd=joink(scratchdir,'tuneHBstrength','out')
collf=joink(scd,'varsets.can')
canlistf=joink(scd,'list') #!!!!!!!! make sure you have this file !!!!!!!
repod=joink(scd,'repod')


def structure2chaindat(structurepdbf):

    """generate a chain.dat file from one pdb file"""
    
    buf='1\n'
    lines=os.popen('cat '+pdbf+' grep ATOM | grep " CA "').readlines()
    buf+=`len(lines)`+'\n'
    for line in lines:
        n=line[6:11]
        x=line[30:38] ; y=line[38:46]; z=line[46:54]
        buf+='%5d %8.3lf %8.3lf %8.3lf\b'%(n,x,y,z)
    return buf

    
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
        return True #True will pickledump


    def addL(self):

        """add sequence length"""

        cmd='wc -l '+joink(self.simdir,'seq.dat')
        self.L=int( os.popen(cmd).readline().split()[0] ) #;print self.L
        return True


    def genChainDat(self):

        """generate chain.dat files beginning with seq.dat files"""

        cmd=joink(pcd,'genChainDat.py')+' -f no -b '+joink(pcd,'betahelix.pdb')+' -a '+joink(self.simdir,'seq.dat')+' -c '+joink(self.simdir,'chain.dat')+' -d '+joink(self.simdir,'init.pdb')
        os.system(cmd)
        return False #there's no update of objects to make
    

    def genarfavors(self,arfavors=['1.0',]):
        
        """initialize or update arfavors dictionary"""

        if not 'arfavors' in self.__dict__: self.arfavors={}

        #generate or update self.arfavors['list']
        if not 'list' in self.arfavors.keys():
            self.arfavors['list']=arfavors[:]
        else:
            list=self.arfavors['list']
            for arfavor in arfavors:
                if arfavor not in list: self.arfavors['list'].append(arfavor)

        #generate new self.arfavors[arfavor] entries and create subdirectories
        for arfavor in self.arfavors['list']:
            if arfavor not in self.arfavors.keys():
                self.arfavors[arfavor]={'simdir':joink(self.simdir,arfavor)}
                os.system('/bin/mkdir -p '+self.arfavors[arfavor]['simdir'])
        return True


    def submitGenericTasser(self,joblist='',arfavors=['1.0',]):

        """submit generic TASSER+SPICKER+AMBERMIN jobs

        submitGenericTasser(id='',arfavors=['1.0',]
    
        we submit as many jobs as elements of arfavors, which are the
        different ratios of long range to short range strengths of the
        hydrogen bond"""

        from jobs.job import job,genJobList
        from tasser.tasser_utilities import calc_wall_time
        
        if not joblist: joblist=genJobList() #;Bye(joblist)
        wt=calc_wall_time(L=self.L) #estimated TASSER wall time, in days
        wt='%5.2lf'%( float(wt)+0.50) #add 1/2 day for spicker plus amber minimizations

        for arfavor in arfavors:
            jobname='x_'+self.id+'_'+arfavor
            outd=joink(self.simdir,arfavor) #output directory
            Job=job(name='generic tasser',exe='submitGenericTasser.sh',exed=pcd,
                    args=' '+self.id+' '+arfavor+' casB_gene3 '+self.simdir)
            Job.qsub(jobname,outd,wallt=wt,mem_limit='990',ddisk='400',joblist=joblist)
        return False #There is no need to pickledump


    def submitMinimize(self,joblist='',arfavors=['1.0',]):

        """submit generic AMBERMIN jobs

        submitMinimize(id='',arfavors=['1.0',]
    
        we submit as many jobs as elements of arfavors, which are the
        different ratios of long range to short range strengths of the
        hydrogen bond"""

        from jobs.job import job,genJobList
        
        if not joblist: joblist=genJobList() #;Bye(joblist)
        wt='0.249' #almost six hours

        for arfavor in arfavors:
            jobname='min_'+self.id+'_'+arfavor
            outd=joink(self.simdir,arfavor) #output directory
            if isthere(self.arfavors[arfavor]['spicker']['outfile']):
                Job=job(name='generic tasser',exe='submitMinimize.sh',exed=pcd,
                    args=' '+self.id+' '+arfavor)
                Job.qsub(jobname,outd,wallt=wt,mem_limit='500',ddisk='300',joblist=joblist)
        return False #There is no need to pickledump


    def createTasserSpickerStructuresDicts(self):

        """initialize 'tasser', 'spicker', and 'structures' dictionaries"""

        from utilities.small_utilities import unTARme,chomp

        pta=self.arfavors
        for arfavor in pta['list']:
            ptb=pta[arfavor]
            for dict in ('tasser','spicker','structures'):
                ptb[dict]={}
                ptb[dict]['outfile']=''
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


    def gatherStructures(self):

        """gather structures from the TASSER simulations"""

        from utilities.small_utilities import unTARme,chomp
    
        for arfavor in self.arfavors['list']: #cycle over all hydrogen bond strengths
            strd=joink(scd,'structures_CA_'+arfavor) #structure repository directory
            p0=self.arfavors[arfavor]['structures']
            if not isthere(p0['outfile']): continue #no minimized structures
            wd=unTARme( p0['outfile'] ) #untar structures to temporary directory
            for strf in p0['list']: #basename of structure files
                f0=joink(wd,strf)   #full atom structure file, origin
                f1=joink(strd,strf) #CA-only structure file, destination
                cmd='grep " CA " '+f0+' > '+f1 #;Bye(cmd)
                os.system(cmd)
            os.system('/bin/rm -r '+wd)
        return False #there is nothing to update the pickle file in this method
    

    def calculateE(self):

        """calculate energies of each structure

        We use reportEab. Energies are calculated assuming arfavor=1.0
        no matter what particular arfavor we are using. The point is
        that we will be able to compare energy values for the hydrogen
        bond across structures that were generated using different
        arfavors. This way we can compare if increasing arfavor
        increases the number of hydrogen bonds."""

        from utilities.small_utilities import unTARme
        
        for arfavor in self.arfavors['list']:
            ptb=self.arfavors[arfavor]
            pta=ptb['structures']
            if pta['outfile']:
                wd=unTARme(pta['outfile']) #untar PDB structures
                currd=os.getcwd() ; os.chdir(wd)
                #bring neccessary input files to run reportEab
                os.system('/bin/cp '+joink(self.simdir,'rmsinp')+' '+
                          joink(self.simdir,'seq.dat')+' '+
                          joink(pcd,'indd',arfavor,'in.dd')+' '+
                          joink(pcd,'reportEab')+' .')
                for basename in pta['list']:
                    open('chaind.dat','w').write(structure2chaindat(basename))
                    os.system('./reportEab')
                    keystring,junk,valuestring=os.popen('energy_terms.dat').readlines()
                    keylist=keystring.split()[1:] #index==0 corresponds to '#' symbol
                    valuelist=valuestring.split()
                    pta[basename]['energy_terms']={} #contain the different energy values
                    ptc=pta[basename]['energy_terms']
                    for i in range(len(keylist)): ptc[keylist[i]]=float(valuelist[i])
                os.chdir(currd); os.system('/bin/rm -r '+wd)


    def foo(margs,mwargs):
        return False #True will pickledump


#################### end of PEPSICAN  ##################
    
    
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    object ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'simdir': directory to store simulation results
    |_'L': sequence length
    |_'arfavors': dictionary containing results on simulations with different
      |           ratios of long versus local hydrogen bond strength
      |_'list': list of arfavors
      |_arfavor: one of the ratios
        |_simdir: directory where simulation results are stored
        |_'tasser': dictionary containing TASSER derived data
        |  |_'outfile': tbz2 file containing filtered TASSER trajectories
        |_'spicker': dictionary containing SPICKER derived datat
        |  |_'outfile': tbz2 file containing some spicker output
        |_'structures': dictionary containing amber minimized structures and derive data
           |_'outfile': tbz2 file containing the structures
           |_'list': list with structure file basenames
           |_basename: dictionary holding info for one structure
             |_'energy_terms': dictionary holding the different values for the energy terms
    """
    print globals()['datatree'].__doc__


def doitall(id):

    """a recount of everything we had done on each varset

    Useful if we screw-up a particular dump file, because we can redo
    it all"""

    c=PEPSI(id=id)
    c.genDumpf(cl.repod)

    methods=('createSimdir','addL','genarfavors','submitGenericTasser',
             'createTasserSpickerStructuresDicts','gatherStructures','calculateE')
    for method in methods:
        print method
        getattr(objSL,method)()
    #finally secure the beast!
    c.pickleDump()


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
        cl.cycleOverList('createSimdir',silent=False)


def addL(id=''):

    """add sequence length"""

    if id:
        c=cl.loadCAN(id)
        if c.addL()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addL',silent=False)
    

def genChainDat(id=''):

    """generate chain.dat files beginning with seq.dat files"""
    
    if id:
        c=cl.loadCAN(id)
        if c.genChainDat()!=False: c.pickleDump()
    else:
        cl.cycleOverList('genChainDat',silent=False)

    
def genarfavors(id='',arfavors=['1.0','2.0','4.0','8.0']):
    
    """initialize or update arfavors dictionary

    genarfavors(id='',arfavors=['1.0','2.0','4.0','8.0']"""

    if id:
        c=cl.loadCAN(id)
        if c.genarfavors(arfavors=arfavors)!=False: c.pickleDump()
    else:
        remaining=len(canlist)
        for c in cl.iterate():
            print remaining,c.id                
            if c.genarfavors(arfavors=arfavors)!=False: c.pickleDump()
            remaining-=1


def submitGenericTasser(id='',arfavors=['1.0','2.0','4.0','8.0']):

    """submit generic TASSER+SPICKER+AMBERMIN jobs

    submitGenericTasser(id='',arfavors=['1.0','2.0','4.0','8.0'])
    
    we submit as many jobs as elements of arfavors, which are the
    different ratios of long range to short range strengths of the
    hydrogen bond"""

    from jobs.job import genJobList
    
    joblist=genJobList() #;Bye(joblist)
    if id:
        c=cl.loadCAN(id)
        if c.submitGenericTasser(joblist=joblist,arfavors=arfavors)!=False: c.pickleDump()
    else:
        remaining=len(canlist)
        for c in cl.iterate():
            print remaining,c.id                
            if c.submitGenericTasser(joblist=joblist,arfavors=arfavors)!=False: c.pickleDump()
            remaining-=1


def submitMinimize(id='',arfavors=['1.0','2.0','4.0','8.0']):

    """submit generic AMBERMIN jobs

    submitMinimize(id='',arfavors=['1.0','2.0','4.0','8.0'])
    
    we submit as many jobs as elements of arfavors, which are the
    different ratios of long range to short range strengths of the
    hydrogen bond"""

    from jobs.job import genJobList
    
    joblist=genJobList() #;Bye(joblist)
    if id:
        c=cl.loadCAN(id)
        if c.submitMinimize(joblist=joblist,arfavors=arfavors)!=False: c.pickleDump()
    else:
        remaining=len(canlist)
        for c in cl.iterate():
            print remaining,c.id                
            if c.submitMinimize(joblist=joblist,arfavors=arfavors)!=False: c.pickleDump()
            remaining-=1


def createTasserSpickerStructuresDicts(id=''):
    
    """create tasser dictionary"""

    if id:
        c=cl.loadCAN(id)
        if c.createTasserSpickerStructuresDicts()!=False: c.pickleDump()
    else:
        cl.cycleOverList('createTasserSpickerStructuresDicts',silent=False)


def calculateE(id=''):
    
    """calculate energies of each structure

    We use reportEab. Energies are calculated assuming arfavor=1.0 no
    matter what particular arfavor we are using. The point is that we
    will be able to compare energy values for the hydrogen bond across
    structures that were generated using different arfavors. This way
    we can compare if increasing arfavor increases the number of
    hydrogen bonds."""

    if id:
        c=cl.loadCAN(id)
        if c.calculateE()!=False: c.pickleDump()
    else:
        cl.cycleOverList('calculateE',silent=False)

    
def gatherStructures(id=''):

    """gather structures from the TASSER simulations"""

    for arfavor in ['1.0','2.0','4.0','8.0']: #init structure directories
        strd=joink(scd,'structures_CA_'+arfavor)
        os.system('/bin/mkdir -p '+strd)

    if id:
        c=cl.loadCAN(id)
        if c.gatherStructures()!=False: c.pickleDump()
    else:
        cl.cycleOverList('gatherStructures',silent=False)
    

def submitTMalign_lists():

    """Do TMalign of all against all for every hydrogen bond strength"""

    from TMalign.TMalign_lists import TMalign_lists
    
    for arfavor in ['1.0','2.0','4.0','8.0']: #init structure directories
        strd=joink(scd,'structures_CA_'+arfavor) #directory of structures
        lib='/tmp/jose/structures_CA_'+arfavor+':'+strd+'.tbz2'
        listf=joink(strd,'list')
        outd=joink(scd,'tmAll2All_'+arfavor)
        os.system('/bin/mkdir -p '+outd)
        TMalign_lists(lib,listf,outd,sizechunk='250',prefix='tm'+arfavor[0]) #qsub jobs
        

def gathertmAll2Alldat():

    """create files tmAll2All_arfavor.dat files"""

    from TMalign.TMalign_lists import gathertmAll2Alldat
    
    for arfavor in ['1.0','2.0','4.0','8.0']:
        tmAll2Alldir=joink(scd,'tmAll2All_'+arfavor)
        tmAll2Allfile=joink(scd,'tmAll2All_'+arfavor+'.dat')
        listofstruct=joink(scd,'structures_CA_'+arfavor,'list')
        tmtmregex='tm'+arfavor[0]+'??_??'
        gathertmAll2Alldat(tmAll2Alldir,tmAll2Allfile,listofstruct,tmtmregex=tmtmregex)
        

def gentmX2Allfiles():

    """create tmX2All_arfavor directories"""

    from TMalign.TMalign_lists import gentmX2Allfiles
    
    for arfavor in ['1.0','2.0','4.0','8.0']:
        tmAll2Allfile=joink(scd,'tmAll2All_'+arfavor+'.dat')
        tmX2all=joink(scd,'tmX2all_'+arfavor)
        strs=[]
        for c in cl.iterate():
            p=c.arfavors[arfavor]['structures']
            if isthere(p['outfile']):
                strs=strs+p['list'] #list of structures
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)

    
def foo(kargs,kwargs):

    """if we have to pass whatever arguments"""
    margs=kargs+largs    #figurative way of adding more arguments
    mwargs=kwargs+lwargs
    if id:
        c=cl.loadCAN(id)
        if c.foo(margs,mwargs)!=False: c.pickleDump()
    else:
        remaining=len(canlist)
        for c in cl.iterate():
            print remaining,c.id                
            if c.foo(margs,mwargs)!=False: c.pickleDump()
            remaining-=1


def help(*kargs):

    """list available services"""
    
    services=['datatree','doitall','updateCollection','initCollection','createSimdir',
              'addL','genarfavors','submitGenericTasser','submitMinimize',
              'createTasserSpickerStructuresDicts','gatherStructures','calculateE',
              'submitTMalign_lists','gathertmAll2Alldat','genChainDat']
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
