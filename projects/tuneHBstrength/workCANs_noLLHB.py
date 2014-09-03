#!/usr/bin/python

import os,sys,re
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir
from utilities.small_utilities import junkName,chomp,Bye
from inputArgs.inputArgs import inpHand
from jobs.job import pastry

"""same as workCANs, but we do TASSER simulations with casB_gene4"""

#global variables
joink=os.path.join
isthere=os.path.exists

pcd=joink(codedir,'projects','tuneHBstrength')
scd=joink(scratchdir,'tuneHBstrength','noLLHB','out')
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


    def addSecondary(self):

        """add secondary structure info"""

        from secondary.get_seq import get_seq
        self.secondary={'seqdat':joink(self.simdir,'seq.dat')}
        self.secondary['seq']=get_seq( self.secondary['seqdat'] )
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


    def submitGenericTasser(self,joblist='',arfavors=['1.0',],
                            protocol='',avHBn=0):

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
        ibatch=0
        if not protocol: protocol='generic tasser'
        if protocol=='generic tasser':
            for arfavor in arfavors:
                jobname='x_'+self.id+'_'+arfavor
                outd=joink(self.simdir,arfavor) #output directory
                Job=job(name='generic tasser',exe='submitGenericTasser.sh',exed=pcd,
                        args=' '+self.id+' '+arfavor+' casB_gene4 '+self.simdir)
                Job.qsub(jobname,outd,wallt=wt,mem_limit='990',ddisk='400',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        elif protocol=='filter by average HB number':
            arfavor='2.0'
            jobname='u_'+self.id+'_'+arfavor
            outd=joink(self.simdir,arfavor) #output directory
            Job=job(name='generic tasser',exe='submitGenericTasserHBfilter.sh',exed=pcd,
                    args=' '+self.id+' '+arfavor+' casB_gene4 '+self.simdir+' '+`avHBn`)
            Job.qsub(jobname,outd,wallt=wt,mem_limit='990',ddisk='400',joblist=joblist)
            ibatch=Job.dormant(ibatch)
        elif protocol=='second neighbors':
            arfavor='2.0' #long-range strength increase
            avHBn='48'    #filter by number of hydrogen bonds
            eh6bs=['1.0','2.0','3.0','4.0','5.0','6.0']
            for eh6b in eh6bs:
                jobname=self.id+'_'+eh6b
                outd=joink(self.simdir,arfavor,eh6b)
                args=' '.join((self.id,arfavor,eh6b,avHBn,'casB_gene5',self.simdir))
                Job=job(name='generic tasser',exe='submitGenericTasserHBfilterSecondNeig.sh',\
                     exed=pcd,args=args)
                Job.qsub(jobname,outd,wallt=wt,mem_limit='990',ddisk='400',joblist=joblist)
                ibatch=Job.dormant(ibatch)
        return False #There is no need to pickledump


    def createTasserSpickerStructuresDicts(self,switch=''):

        """initialize 'tasser', 'spicker', and 'structures' dictionaries"""

        from utilities.small_utilities import unTARme,chomp
        if not switch:
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
        elif switch=='avHBn':
            self.avHBn={}
            if self.id[0]=='b': #reference to previous simulation for beta structures
                self.avHBn=self.arfavors['2.0']
            elif self.id[0]=='c': #store info for new simulations done for alphabeta
                ptb=self.avHBn
                ptb['simdir']=joink(self.simdir,'2.0')
                for dict in ('tasser','spicker','structures'):
                    ptb[dict]={}
                    ptb[dict]['outfile']=''
                    outfile=joink(ptb['simdir'],'out.'+dict+'.avHBn.tbz2')
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
        if switch=='second neighbors':
            self.secondneig={'list':['1.0','2.0','3.0','4.0','5.0','6.0'],}
            for eh6b in self.secondneig['list']:
                self.secondneig[eh6b]={}
                pt1=self.secondneig[eh6b]
                pt1['simdir']=joink(self.simdir,'2.0',eh6b)
                for dict in ('tasser','spicker','structures'):
                    pt1[dict]={}
                    pt1[dict]['outfile']=''
                    outfile=joink(pt1['simdir'],'out.'+dict+'.tbz2')
                    if isthere(outfile): pt1[dict]['outfile']=outfile
                if pt1['structures']['outfile']:
                    pt2=pt1['structures']
                    currd=os.getcwd()
                    wd=unTARme( pt2['outfile'] ) #untar structures to temporary directory
                    os.chdir(wd)
                    pt2['list']=chomp(os.popen('ls -1 *.pdb').readlines())
                    os.chdir(currd)
                    os.system('/bin/rm -r '+wd)
                    for basename in pt2['list']: pt2[basename]={}
        return True


    def gatherStructures(self,switch=''):

        """gather structures from the TASSER simulations"""

        from utilities.small_utilities import unTARme,chomp

        if not switch:
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
        elif switch=='avHBn':
            strd=joink(scd,'structures_CA_avHBn')
            if not isthere(strd): os.system('/bin/mkdir -p '+strd)
            p0=self.avHBn['structures']
            if not isthere(p0['outfile']): return False #no minimized structures
            wd=unTARme( p0['outfile'] ) #untar structures to temporary directory
            for strf in p0['list']: #basename of structure files
                f0=joink(wd,strf)   #full atom structure file, origin
                f1=joink(strd,strf) #CA-only structure file, destination
                cmd='grep " CA " '+f0+' > '+f1 #;Bye(cmd)
                os.system(cmd)
            os.system('/bin/rm -r '+wd)
        elif switch=='second neighbors':
            arfavor='2.0' #long-range strength increase
            eh6bs=['1.0','2.0','3.0','4.0','5.0','6.0']
            for eh6b in eh6bs:
                strd=joink(scd,'structures_CA_sn'+eh6b)
                if not isthere(strd): os.system('/bin/mkdir -p '+strd)
                p1=self.secondneig[eh6b]['structures']
                if not isthere(p1['outfile']): continue #no minimized structures
                wd=unTARme( p1['outfile'] ) #untar structures to temporary directory
                for strf in p1['list']: #basename of structure files
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


    def secondaryContent(self):

        """calculate secondary content"""

        from secondary.seqdat2ss import translateSec,nn2ss
        self.secondary={'str':'','content':{'H':0.0,'E':0.0,'C':0.0}}
        self.secondary['str']=translateSec(joink(self.simdir,'seq.dat'))
        str=self.secondary['str']
        patt={'H':re.compile('H'),'E':re.compile('E')}
        content=self.secondary['content']
        buf=self.id
        for X in ('H','E'):
            content[X]=float(len(patt[X].findall(str)))/self.L
            buf+='  %4.2lf'%content[X]
        content['C']=1-content['H']-content['E']
        buf+='  %4.2lf'%content['C']
        buf+='  %4.2lf'%(1-content['C'],)
        print buf


    def hbnet(self):

        """calculate hydrogen bonds for structures created"""

        hbexe='/gpfs1/active/jose/code/f77/tertiary/hbnet.x'
        for arfavor in self.arfavors['list']:
            pta=self.arfavors[arfavor]['structures']
            for basename in pta['list']:
                pta[basename]['hbnetf']=''
                pdbf=joink(scd,'structures_CA_'+arfavor,basename)
                if isthere(pdbf):
                    secf=joink(self.simdir,'seq.dat')
                    outf=joink(self.simdir,arfavor,basename[0:-3]+'hb')
                    cmd=hbexe+' -pdbf '+pdbf+' -secf '+secf+' -outf '+outf
                    os.system(cmd)
                    pta[basename]['hbnetf']=outf
        return True


    def nhb(self):

        """calculate number of hydrogen bonds and total energy per structure"""

        for arfavor in self.arfavors['list']:
            pta=self.arfavors[arfavor]['structures']
            for basename in pta['list']:
                if pta[basename]['hbnetf']:
                    pta[basename]['hb']={}
                    ptb=pta[basename]['hb']
                    hblines=open(pta[basename]['hbnetf'],'r').readlines()
                    ptb['N']=len(hblines)
                    ptb['E']=0.0
                    for hbline in hblines: ptb['E']+=float(hbline.split()[2])
                    ptb['e']=ptb['E']/ptb['N']
                    print '%s %s %12s %4d %7.2lf %6.3lf'%\
                          (self.id,arfavor,basename,ptb['N'],ptb['E'],ptb['e'])
        return True

    
    def foo(margs,mwargs):
        return False #True will pickledump


#################### end of PEPSICAN  ##################
    
    
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    ojbect ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'simdir': directory to store simulation results
    |_'L': sequence length
    |_'secondary': info on secondary structure
    | |_'seqdat': seq.dat file
    | |_'seq': secondary structure asignment in a string (H,E,C)
    |_'arfavors': dictionary containing results on simulations with different
    | |           ratios of long versus local hydrogen bond strength
    | |_'list': list of arfavors
    | |_arfavor: one of the ratios
    |   |_simdir: directory where simulation results are stored
    |   |_'tasser': dictionary containing TASSER derived data
    |   |  |_'outfile': tbz2 file containing filtered TASSER trajectories
    |   |_'spicker': dictionary containing SPICKER derived datat
    |   |  |_'outfile': tbz2 file containing some spicker output
    |   |_'structures': dictionary containing amber minimized structures and derive data
    |     |_'outfile': tbz2 file containing the structures
    |     |_'list': list with structure file basenames
    |     |_basename: dictionary holding info for one structure
    |       |_'hbnet': hydrogen bonds file
    |       |_'hb': dictionary related to hydrogen bonding
    |       | |_'N': number of hydrogen bonds (but they have different energies)
    |       | |_'E': total hydrogen bond energy
    |       | |_'e': hydrogen bond energy per bond
    |       |_'energy_terms': dictionary holding the different values for the energy terms
    |_'avHBn': results from simulations with filtering of number of hydrogen bonds. Done
    |  |        with twice strength for long-range hydrogen bond
    |  |_'tasser': dictionary containing TASSER derived data
    |  |  |_'outfile': tbz2 file containing filtered TASSER trajectories
    |  |_'spicker': dictionary containing SPICKER derived datat
    |  | |_'outfile': tbz2 file containing some spicker output
    |  |_'structures': dictionary containing amber minimized structures and derive data
    |    |_'outfile': tbz2 file containing the structures
    |    |_'list': list with structure file basenames
    |    |_basename: dictionary holding info for one structure
    |_'secondneig': results from simulations with list of hydrogen bonds and second
      |        neighbor interactions. It also includes trajectory filtering for
      |        number of hydrogen bonds
      |_'list': list of second neighbor interaction weight values
      |_eh6b: one of the second neighbor interaction weight values
        |_'simdir': directory where simulation results are stored
        |_'tasser': dictionary containing TASSER derived data
        |  |_'outfile': tbz2 file containing filtered TASSER trajectories
        |_'spicker': dictionary containing SPICKER derived datat
        | |_'outfile': tbz2 file containing some spicker output
        |_'structures': dictionary containing amber minimized structures and derive data
          |_'outfile': tbz2 file containing the structures
          |_'list': list with structure file basenames
          |_basename: dictionary holding info for one structure
          
    """
    print globals()['datatree'].__doc__


def doitall(id):

    """a recount of everything we had done on each varset

    Useful if we screw-up a particular dump file, because we can redo
    it all"""

    c=PEPSI(id=id)
    c.genDumpf(cl.repod)

    methods=('createSimdir','addL','genarfavors','submitGenericTasser',
             'createTasserSpickerStructuresDicts','gatherStructures','calculateE',
             'secondaryContent','hbnet','nhb')
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


def submitGenericTasser(id='',arfavors=['1.0','2.0','4.0','8.0'],
                        protocol='',avHBn=''):

    """submit generic TASSER+SPICKER+AMBERMIN jobs

    submitGenericTasser(id='',arfavors=['1.0','2.0','4.0','8.0'])
    
    we submit as many jobs as elements of arfavors, which are the
    different ratios of long range to short range strengths of the
    hydrogen bond"""

    from jobs.job import genJobList

    if not protocol: protocol='generic tasser'
    joblist=genJobList() #;Bye(joblist)
    if id:
        c=cl.loadCAN(id)
        if c.submitGenericTasser(joblist=joblist,arfavors=arfavors)!=False: c.pickleDump()
    else:
        remaining=len(canlist)
        if protocol=='generic tasser':
            for c in cl.iterate():
                print remaining,c.id                
                if c.submitGenericTasser(joblist=joblist,arfavors=arfavors)!=False: c.pickleDump()
                remaining-=1
        elif protocol=='filter by average HB number':
            if not avHBn:
                sys.stderr.write('ERROR: must pass argument avHBn\n')
            else: avHBn=int(avHBn)
            for id in canlist:
                if id[0]=='c': #only for alphabeta targets
                    c=cl.loadCAN(id)
                    if c.submitGenericTasser(joblist=joblist,arfavors=['2.0',],
                                             protocol=protocol,avHBn=avHBn)!=False: c.pickleDump()
                remaining-=1
        elif protocol=='second neighbors':
            for c in cl.iterate():
                if c.submitGenericTasser(joblist=joblist,protocol=protocol)!=False: c.pickleDump()
                remaining-=1


def createTasserSpickerStructuresDicts(id='',switch=''):
    
    """create tasser dictionary"""

    if id:
        c=cl.loadCAN(id)
        if c.createTasserSpickerStructuresDicts(switch=switch)!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            print remaining,c.id
            if c.createTasserSpickerStructuresDicts(switch=switch)!=False: c.pickleDump()
            remaining-=1


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

    
def gatherStructures(id='',switch=''):

    """gather structures from the TASSER simulations"""

    if not switch:
        for arfavor in ['1.0','2.0','4.0','8.0']: #init structure directories
            strd=joink(scd,'structures_CA_'+arfavor)
            os.system('/bin/mkdir -p '+strd)

    if id:
        c=cl.loadCAN(id)
        if c.gatherStructures(switch=switch)!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            if c.gatherStructures(switch=switch)!=False: c.pickleDump()
            print remaining,c.id
            remaining-=1
    

def submitTMalign_lists(switch='',wallt=None):

    """Do TMalign of all against all for every hydrogen bond strength"""

    from TMalign.TMalign_lists import TMalign_lists,TMalign_listsII
    
    if not switch:
        for arfavor in ['1.0','2.0','4.0','8.0']: #init structure directories
            strd=joink(scd,'structures_CA_'+arfavor) #directory of structures
            lib='/tmp/jose/structures_CA_'+arfavor+':'+strd+'.tbz2'
            listf=joink(strd,'list')
            outd=joink(scd,'tmAll2All_'+arfavor)
            os.system('/bin/mkdir -p '+outd)
            TMalign_lists(lib,listf,outd,prefix='tm'+arfavor[0],wallt=wallt) #qsub jobs
    elif switch=='avHBn':
        strd=joink(scd,'structures_CA_avHBn') #directory of structures
        lib='/tmp/jose/structures_CA_avHBn:'+strd+'.tbz2'
        listf=joink(strd,'list')
        outd=joink(scd,'tmAll2All_avHBn')
        os.system('/bin/mkdir -p '+outd)
        TMalign_lists(lib,listf,outd,prefix='tmav',exe_type='shashi',wallt=wallt) #qsub jobs
    elif switch=='alpha_too':
        strd1=joink(scratchdir,'tuneHBstrength','structures_CA_alpha')
        lib1='/tmp/jose/structures_CA_alpha:'+strd1+'.tbz2'
        listf1=joink(strd1,'list')
        strd2=joink(scd,'structures_CA_avHBn')
        lib2='/tmp/jose/structures_CA_avHBn:'+strd2+'.tbz2'
        listf2=joink(strd2,'list')
        outd=joink(scd,'tmAll2All_alpha_too') ;os.system('/bin/mkdir -p '+outd)
        TMalign_listsII(lib1,listf1,lib2,listf2,outd,prefix='tmx',exe_type='shashi',wallt=wallt)
    elif switch=='second neighbors alphabeta only':
        for eh6b in ['1.0','2.0','3.0','4.0','5.0','6.0']:
            strd=joink(scd,'structures_CA_sn'+eh6b)
            lib='/tmp/jose/structures_CA_sn'+eh6b+':'+strd+'.tbz2'
            junkf=junkName()
            listf=joink(strd,'c.list')
            os.system('ls -1 '+strd+'|grep ^c >'+junkf+' && /bin/mv '+junkf+' '+listf)
            outd=joink(scd,'tmAll2All_sn'+eh6b) ;os.system('/bin/mkdir -p '+outd)
            TMalign_lists(lib,listf,outd,prefix='tm'+eh6b,exe_type='shashi',wallt=wallt)
    elif switch=='abc':
        strd=joink(scd,'structures_CA_abc')
        lib='/tmp/jose/structures_CA_abc:'+strd+'.tbz2'
        listf=joink(strd,'list')
        outd=joink(scd,'tmAll2All_sn') ;os.system('/bin/mkdir -p '+outd)
        TMalign_lists(lib,listf,outd,prefix='tm',exe_type='fragment',wallt=wallt)
    elif switch=='origPDB': #compare original set of structures with PDB subset
        if not wallt: wallt=0.249
        strd=joink(scratchdir,'Hellinga/distrByLengthDecreasedLoopLoopInteraction',\
                   'linker42sum/structures_CA_only')
        lib='/tmp/jose/structures_CA_only:'+strd+'.tbz2'
        listf=joink(strd,'abc.list')
        strd2=joink(scratchdir,'structural_neighbors/pdb')
        lib2='/tmp/jose/pdb:'+strd2+'.tbz2'
        listf2=joink(strd2,'100_200.list')
        outd=joink(scd,'tmAll2All_origPDB') ;os.system('/bin/mkdir -p '+outd)
        TMalign_listsII(lib,listf,lib2,listf2,outd,\
                        prefix='tmo',exe_type='fragment',sizechunk=100,wallt=wallt)
    elif switch=='abcPDB': #compare original set of structures with PDB subset
        if not wallt: wallt=0.249
        strd=joink(scd,'structures_CA_abc')
        lib='/tmp/jose/structures_CA_abc:'+strd+'.tbz2'
        listf=joink(strd,'list')
        strd2=joink(scratchdir,'structural_neighbors/pdb')
        lib2='/tmp/jose/pdb:'+strd2+'.tbz2'
        listf2=joink(strd2,'100_200.list')
        outd=joink(scd,'tmAll2All_abcPDB') ;os.system('/bin/mkdir -p '+outd)
        TMalign_listsII(lib,listf,lib2,listf2,outd,\
                        prefix='tma',exe_type='fragment',sizechunk=100,wallt=wallt)
        
        
def gathertmAll2Alldat():

    """create files tmAll2All_arfavor.dat files"""

    from TMalign.TMalign_lists import gathertmAll2Alldat
    
    for arfavor in ['1.0','2.0','4.0','8.0']:
        tmAll2Alldir=joink(scd,'tmAll2All_'+arfavor)
        tmAll2Allfile=joink(scd,'tmAll2All_'+arfavor+'.dat')
        listofstruct=joink(scd,'structures_CA_'+arfavor,'list')
        tmtmregex='tm'+arfavor[0]+'??_??'
        gathertmAll2Alldat(tmAll2Alldir,tmAll2Allfile,listofstruct,tmtmregex=tmtmregex)
        

def gentmX2Allfiles(switch=''):

    """create tmX2All_arfavor directories"""

    from TMalign.TMalign_lists import gentmX2Allfiles
    from utilities.small_utilities import chomp
    
    tmco=0.4
    if not switch:
        for arfavor in ['1.0','2.0','4.0','8.0']:
            tmAll2Allfile=joink(scd,'tmAll2All_'+arfavor+'.dat')
            tmX2all=joink(scd,'tmX2all_'+arfavor)
            strs=[]
            for c in cl.iterate():
                p=c.arfavors[arfavor]['structures']
                if isthere(p['outfile']): strs=strs+p['list'] #list of structures
            gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
    elif switch=='avHBn':
        tmAll2Allfile=joink(scd,'tmAll2All_avHBn.dat')
        tmX2all=joink(scd,'tmX2all_avHBn')
        strs=[]
        for c in cl.iterate():
            p=c.avHBn['structures']
            if isthere(p['outfile']): strs=strs+p['list'] #list of structures
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
    elif switch=='avHBn_and_alpha':
        tmAll2Allfile=joink(scd,'tmAll2All_avHBn_alpha.dat') #all TM-TM comparisons
        listfile=joink(scd,'avHBn_and_alpha.list')
        strs=chomp(open(listfile,'r').readlines())
        tmX2all=joink(scd,'tmX2all_avHBn_alpha')  #output directory
        strs=chomp(open(joink(scd,'avHBn_and_alpha.list'),'r').readlines()) #list of structures
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
    elif switch=='second neighbors alphabeta only':
        for eh6b in ['1.0','2.0','3.0','4.0','5.0','6.0']:
            tmAll2Allfile=joink(scd,'tmAll2All_sn'+eh6b+'.c.dat') #all TM-TM comparisons
            tmX2all=joink(scd,'tmX2all_sn'+eh6b+'.c')  #output directory
            if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
            strs=[]
            for c in cl.iterate():
                if c.id[0]=='c': #do only for alphabeta structures
                    p=c.secondneig[eh6b]['structures']
                    if isthere(p['outfile']): strs=strs+p['list'] #list of structures
            open(joink(tmX2all,'c.list'),'w').write('\n'.join(strs))
            gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
    elif switch=='abc':
        tmAll2Allfile=joink(scd,'tmAll2All_sn.dat') #all TM-TM comparisons
        tmX2all=joink(scd,'tmX2all_sn')  #output directory
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strs=chomp(open(joink(scd,'structures_CA_abc/list'),'r').readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
    elif switch=='abcPDB':
        tmAll2Allfile=joink(scd,'tmAll2All_abcPDB.dat') #all TM-TM comparisons
        tmX2all=joink(scd,'tmX2all_abcPDB')  #output directory
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strd=joink(scd,'structures_CA_abc')
        strd2=joink(scratchdir,'structural_neighbors/pdb')        
        strs=chomp(open(joink(strd,'list'),'r').readlines())+\
              chomp(open(joink(strd2,'100_200.list'),'r').readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
    elif switch=='origPDB':
        tmAll2Allfile=joink(scd,'tmAll2All_origPDB.dat') #all TM-TM comparisons
        tmX2all=joink(scd,'tmX2all_origPDB')  #output directory
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strd=joink(scratchdir,'Hellinga/distrByLengthDecreasedLoopLoopInteraction',\
                   'linker42sum/structures_CA_only')
        strd2=joink(scratchdir,'structural_neighbors/pdb')        
        strs=chomp(open(joink(strd,'abc.list'),'r').readlines())+\
              chomp(open(joink(strd2,'100_200.list'),'r').readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)

        
def secondaryContent(id=''):

    """calculate secondary content"""

    print '# id     H     E     C    H+E'
    if id:
        c=cl.loadCAN(id)
        if c.secondaryContent()!=False: c.pickleDump()
    else:
        cl.cycleOverList('secondaryContent',silent=False)


def mixsecondaryneighbors():

    """add secondary info to file nneighbors.dat"""
    
    nneig=joink(scd,'tmX2all_2.0','nneighbors.dat')
    nneigmore=joink(scd,'tmX2all_2.0','nneighbors_sec.dat')
    outp=open(nneigmore,'w')
    outp.write('# structure   nneig  H     E     C     H+E\n')
    for line in chomp(open(nneig,'r').readlines()):
        items=line.split()
        id=items[0][0:5]
        content=cl.loadCAN(id).secondary['content'] #secondary content
        buf=''
        for X in ('H','E','C'): buf+='  %4.2lf'%content[X]
        buf+='  %4.2lf'%(1-content['C'])
        outp.write(line+buf+'\n')


def is_structures_CA_2_in_PDB():

    """Check if out/structures_CA_2.0 is in PDB /gpfs1/scratch/jose/db/pdb_mar07.tbz2

    Jeff wonders if singletons of out/structures_CA_2.0 have hits in the PDB"""
    
    from TMalign.TMalign_lists import TMalign_listsII

    arfavor='2.0'
    strd=joink(scd,'structures_CA_'+arfavor) #directory of structures
    lib='/tmp/jose/structures_CA_'+arfavor+':'+strd+'.tbz2'
    listf=joink(strd,'list')

    pdblib='/tmp/jose/pdb_mar07:/gpfs1/scratch/jose/db/pdb_mar07.tbz2'
    pdblist='/gpfs1/scratch/jose/db/pdb_mar07.list'

    outd=joink(scd,'structures_CA_'+arfavor+'_in_PDB')
    os.system('/bin/mkdir -p '+outd)
    
    TMalign_listsII(lib,listf,pdblib,pdblist,outd,prefix='',exe_type='yang02')
    

def filter_singletons_2():

    """do Shashi's TMalign on singletons in structures_CA_2.0

    We obtained a list of singletons using Yang's TMalign in
    structures_CA_2.0, but maybe we can elminate some of them using
    Shashi's TMalign"""

    from TMalign.TMalign_lists import TMalign_listsII

    arfavor='2.0'
    strd=joink(scd,'structures_CA_'+arfavor) #directory of structures
    lib='/tmp/jose/structures_CA_'+arfavor+':'+strd+'.tbz2'
    listf=joink(strd,'list')

    lib2=lib
    list2=joink(strd,'yang_singletons.list')

    outd=joink(scd,'tmAll2All_'+arfavor+'_ys')
    os.system('/bin/mkdir -p '+outd)
    
    TMalign_listsII(lib,listf,lib2,list2,outd,prefix='stm',exe_type='shashi')
    

def hbnet(id=''):

    """calculate hydrogen bonds for structures created"""

    if id:
        c=cl.loadCAN(id)
        if c.hbnet()!=False: c.pickleDump()
    else:
        cl.cycleOverList('hbnet',silent=False)


def nhb(id=''):

    """calculate number of hydrogen bonds and total energy per structure"""
    print '#N: number of hydrogen bonds'
    print '#E: total enery of hydrogen bonds'
    print '#n: hydrogen bond energy per bond'
    print '#id arfavor  structure    N    E       n'
    print 'b0034 8.0 b0034.10.pdb   44  -25.31 -0.575'
    if id:
        c=cl.loadCAN(id)
        if c.nhb()!=False: c.pickleDump()
    else:
        cl.cycleOverList('nhb',silent=True)


def browser(id):

    """ browse attributes of a particular object """

    from dataStructures.browser import Browser
    import gtk
    
    window = gtk.Window()
    window.set_title('Interactive GTK Object Browser')
    window.set_default_size (350,400)
    window.set_border_width (8)
    window.set_name (id)
    vbox = gtk.VBox(False, 8)
    window.add(vbox)
    label = gtk.Label ("Interactive GTK Object Browser. Copyright (c) 2006 Nicolas Rougier")
    vbox.pack_start (label, False, False)

    # Dummy selection callback
    def print_selection (child, name, fullname):
        if child is not None:
            print 'New selection :', name, fullname
        else:
            print 'No selection'

    def predicate (obj):
        return isinstance (obj,PEPSI)

    browser = Browser ((cl.loadCAN(id), ), predicate, print_selection)
	
    vbox.pack_start (browser)
    window.connect ('delete_event', lambda w,x: gtk.main_quit())
    window.show_all()
    gtk.main()


def view(id):

    """primitive viewer of object properties

    workCANs_noLLHB.py -a view -b id
      id: one PDB id"""
    
    cl.loadCAN(id).view()


def singleton_content_histogram(switch=''):

    """cummulative histogram of singleton content in a given secondary structure assignment

    There are ten clusters per secondary structure assignment. The singleton content is
    the number of singletons out of these ten clusters."""

    nneighbors_0f=''
    if switch=='original':
        nneighbors_0f=joink(scratchdir,'Hellinga/distrByLengthDecreasedLoopLoopInteraction',\
                   'linker42sum/tmAll2All/tmX2all/nneighbors.abc.0.dat')
    elif switch=='abc':
        nneighbors_0f=joink(scd,'tmX2all_sn/nneighbors_0.dat')
    #fill dictionary of secondary structurea assignments
    secasgs={}  ; nsingl=0
    for line in open(nneighbors_0f,'r').readlines():
        secasg=line[0:5]
        if secasg not in secasgs: secasgs[secasg]=1
        else: secasgs[secasg]+=1
        nsingl+=1
    #create cummulative singleton content histogram
    maxcontent=10 ; hist=[0]*(maxcontent+1)
    for content in secasgs.values(): hist[content]+=content
    for content in range(1,maxcontent): hist[maxcontent-content]+=hist[maxcontent-content+1]
    #output
    for content in range(1,maxcontent):
        print '%02d %4.2lf'%(content,float(hist[content])/nsingl)


def edge_ratio(switch=''):

    """ratio of internal edges to combined internal plus external edges
    
    From the set of singletons in the original set, we can obtain the set of secondary
    structural classes (SSC) they belong to. Now, for one of these SSC's, we can span a
    network where nodes are the ten clusters, we connect two nodes if their
    TM-score>0.4. These are edges internal to the SSC. We can think of external edges linking
    the 10 clusters in the SSC to the remaining structures of the original set. Now, after
    inclusion of the second order interaction in the hydrogen bond, we can again compute the
    number of such internal and external edges for the same SSC. We report the ratio of the
    change in the number of internal nodes for all SSC's in the set to the combined change in
    the number of internal and external edges for all SSC's in the set.
    $CDP/workCANs_noLLHB.py -a edge_ratio"""


    dirs={'original':joink(scratchdir,'Hellinga/distrByLengthDecreasedLoopLoopInteraction',\
                              'linker42sum/tmAll2All'),
          'abc':scd }
    tmAll2All={'original':joink(dirs['original'],'tmAll2All.dat'),
               'abc':joink(dirs['abc'],'tmAll2All_sn.dat')}
    nneighbors_0f={'original': joink(dirs['original'],'tmX2all/nneighbors.abc.0.dat'),
                   'abc':joink(dirs['abc'],'tmX2all_sn/nneighbors_0.dat')
                   }
    ninternal={'original':0,'abc':0}
    nall={'original':0,'abc':0}

    tmco=0.4
    #obtain set of SSC's
    sscs=[]
    for line in open(nneighbors_0f['original'],'r').readlines():
        ssc=line[0:5]
        if ssc not in sscs: sscs.append(ssc)
    #obtain number of edges for structures belonging to the SSC's in sscs list
    for run in ('original','abc'):
        if not isthere(tmAll2All[run]): Bye('ERROR file '+tmAll2All[run]+'not found')
        pin=open(tmAll2All[run],'r') ; line=pin.readline()
        while line:
            items=line.split() ; ssc0=items[0][0:5] ; ssc1=items[1][0:5] ; tm=float(items[2])
            if tm>tmco and ssc0[0]!='d' and ssc1[0]!='d': #structurally similar, not "d" category
                if ssc0==ssc1: #internal node
                    if ssc0 in sscs:
                        ninternal[run]+=1 ; nall[run]+=1
                else:
                    if (ssc0 in sscs) or (ssc1 in sscs): nall[run]+=1
            line=pin.readline()
        print ninternal[run],nall[run]
    

def addSecondary(id=''):

    """add secondary structure info"""

    if id:
        c=cl.loadCAN(id)
        if c.addSecondary()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addSecondary',silent=False)


def checkSomething(id=''):

    idlist=cl.canlist
    if id: idlist=[id,]
    for id  in idlist:
        c=cl.loadCAN(id)
        try:
            print c.id,len(c.arfavors['list']['structures']['list'])
        except:
            pass  #            sys.stderr.write(c.id+'\n')

        
def help(*kargs):

    """list available services"""
    
    services=['datatree','doitall','updateCollection','initCollection','createSimdir',
              'addL','genarfavors','submitGenericTasser','createTasserSpickerStructuresDicts',
              'gatherStructures','calculateE','submitTMalign_lists','gathertmAll2Alldat',
              'gentmX2Allfiles','secondaryContent','is_structures_CA_2_in_PDB',
              'nhb','genChainDat','browser','view','singleton_content_histogram',
              'addSecondary','checkSomething']
    os.system('/usr/bin/clear')
    if not kargs:
        print 'Services='
        for service in services:
            print '  '+service+' : '+globals()[service].__doc__.split('\n')[0]
        print 'NOTE: type "workCANs.py -a help -b service" for more detail on "service"'
    else:
        sv=kargs[0]
        print globals()[sv].__doc__        


if __name__=='__main__':
    inpHand('Usage: workCANs.py [options]',
            ' -a _R_service required service. Type "-a help" for list of services.',              
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
    locals()[service](*reqargs,**optsargs
)
