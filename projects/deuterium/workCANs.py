#!/usr/bin/python

import os,sys,re
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir,projectsdir
from utilities.small_utilities import chomp,Bye,junkName
from inputArgs.inputArgs import inpHand
from time import sleep
from random import random
"""project description goes here"""

atomic_weight={'H':'1', 'D':'2', 'C':'12', 'N':'14', 'O':'16', 'S':'32'}
def genList():
    """create inital list file"""
    from amber.amber9 import tleap
    open(joink(pjd,'list'),'w').write('\n'.join(tleap().AAA))


#global variables joink=os.path.join
isthere=os.path.exists
joink=os.path.join

my_project='deuterium'
pcd=joink(codedir,'projects',my_project)
pjd=joink(projectsdir,my_project)
pjdout=joink(projectsdir,my_project,'out')
scd=joink(scratchdir,my_project,'out')
scdout=joink(scratchdir,my_project,'out')
collf=joink(pjd,'varsets.can')
canlistf=joink(pjd,'list')
if not isthere(canlistf): genList()
repod=joink(pjd,'repod')

for dir in (pcd,pjdout,scdout,repod):
    if not isthere(dir): os.system('/bin/mkdir -p '+dir)
    

class PEPSIl(CANl):

    """inherited class from CANl, so that we can add methods to fuzz
    around the attributes of our particular CANl object
    implementation"""

    pass

#list of CAN objects for my project
if not isthere(repod): os.system('/bin/mkdir -p '+repod)
cl=PEPSIl(collf,canlistf,repod)
canlist=chomp(open(cl.canlistf,'r').readlines())

def zPdbConvoluteAtomName(name,unit):
    """amber9 to PDB atom name mapping for hydrogens. Transcribed from
    within $AMBERHOME/src/leap/src/leap/pdbFile.c"""
    if len(name)==4: return name[-1]+name[:-1]
    if name[-1] in ('1','2','3') and \
           (name[1]=='B' or \
            (not unit in ('PHE','TYR','TRP') and unit[0:2]!='HI') \
            ):
        return name[-1]+name[:-1]
    return name
    
def name2type(unit):
    pairs={}
    pat=re.compile('H\w\d')
    amino94f=joink(os.environ['AMBERHOME'],'dat/leap/lib/all_amino94.lib')
    amino94=open(amino94f,'r') ; line=amino94.readline()
    while '!entry.'+unit+'.unit.atoms ' not in line:line=amino94.readline()
    line=amino94.readline() #;Bye(line)
    while '!entry.'+unit+'.unit.atomsp' not in line:
        pair=line.split()[:3] ; name=pair[0][1:-1] ; type=pair[1][1:-1]
        match=pat.match(name)
        if match: name=zPdbConvoluteAtomName(name,unit)
        if len(type)==1: type=type+' '
        pairs[name]=type
        line=amino94.readline()
    return pairs

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

    def gensimdir(self):
        """mkdir every simdir"""

        self.psimdir=joink(pjdout,self.id)   #projects directory
        self.ssimdir=joink(scdout,self.id)   #scratch directory
        for dir in (self.psimdir,self.ssimdir):
            os.system('/bin/mkdir -p '+dir)
        return True
        
    def isstandarAAA(self):
        """define if molecule is an amber9 amino acid"""
        from amber.amber9 import leap
        lp=leap()
        self.standarAAA=False
        if self.id in lp.AAA: self.standarAAA=True
        return True

    def genhydrogenated(self):
        """create hydrogenated directory"""
        self.hydrd={'simdir':joink(self.psimdir,'hydrogenated'),
                    'optPDB':None}
        os.system('/bin/mkdir -p '+self.hydrd['simdir'])
        return True

    def gendeuterated(self):
        """create deuterated directory"""
        self.deutd={'simdir':joink(self.psimdir,'deuterated')}
        os.system('/bin/mkdir -p '+self.deutd['simdir'])
        return True

    def geninitPDB(self,switch=''):
        """print PDB, TOP, and CRD of the molecule as initial geometry
        switch='spit script' will spit the tleap script to STDOUT"""
        if self.standarAAA:
            pdbf=joink(self.psimdir,self.id+'.pdb');self.initPDB=None
            topf=joink(self.psimdir,self.id+'.top');self.initTOP=None
            crdf=joink(self.psimdir,self.id+'.crd');self.initCRD=None
            buf='source leaprc.ff98\n'+\
                 'mdl=sequence{ACE '+self.id+' NME}\n'+\
                 'savepdb mdl '+pdbf+'\n'\
                 'saveAmberParm mdl '+topf+' '+crdf+'\n'+\
                 'quit\n'
            if switch=='spit script':
                print '\n'.join(cmds)
                return False
            junkf=junkName() ; open(junkf,'w').write(buf)
            os.system('tleap -f '+junkf+' &>/dev/null && /bin/rm '+junkf)
            self.initPDB=pdbf ; self.initTOP=topf ; self.initCRD=crdf
        elif self.id in ('SPC',):
            pdbf=joink(self.psimdir,self.id+'.pdb');self.initPDB=None
            topf=joink(self.psimdir,self.id+'.top');self.initTOP=None
            crdf=joink(self.psimdir,self.id+'.crd');self.initCRD=None
            buf='source leaprc.ff98\n'+\
                 'mdl=WAT\n'+\
                 'savepdb mdl '+pdbf+'\n'\
                 'saveAmberParm mdl '+topf+' '+crdf+'\n'+\
                 'quit\n'
            junkf=junkName() ; open(junkf,'w').write(buf)
            os.system('tleap -f '+junkf+' &>/dev/null && /bin/rm '+junkf)
            self.initPDB=pdbf ; self.initTOP=topf ; self.initCRD=crdf
        elif self.id=='H2':
            pdbf=joink(self.psimdir,self.id+'.pdb');self.initPDB=None
            buf="""ATOM      1  H1  HY2     1       0.000   0.000   0.000  1.00  0.00
ATOM      2  H2  DE2     1       1.000   0.000   0.000  1.00  0.00
TER
END"""
            open(pdbf,'w').write(buf) ; self.initPDB=pdbf 
        return True

    def rotatePDB(self,deg=1):
        """rotate around the CA-CB bond"""
        pk='/'+self.id+'///'+self.id+'`2/'
        buf=['from pymol import cmd',
             'cmd.load(\''+self.initTOP+'\',format=\'top\')',
             'cmd.load(\''+self.initPDB+'\',format=\'pdb\')']
        if self.id=='THR':
            buf+=['cmd.unbond(\'/THR///ACE`1/O\',\'/THR///THR`2/1HG2\')',]
        buf+=['cmd.bond(\''+pk+'CA\',\''+pk+'CB\')',
             'cmd.edit(\''+pk+'CA\',\''+pk+'CB\')',
             'cmd.torsion('+`deg`+')',
             'cmd.save(\''+self.initPDB+'\',format=\'pdb\')']
        junkf=junkName()+'.py' ; open(junkf,'w').write('\n'.join(buf))
        os.system('pymol -qc '+junkf+' &>/dev/null && /bin/rm '+junkf)
        return False
        
    def checkForClashes(self,switch='init'):
        """check a PDB"""
        pdbdict={'init':self.initPDB,'opt':self.hydrd['optPDB']}
        pdb=pdbdict[switch]  ;  isclash=False
        if pdb:
            cmds=['source leaprc.ff98',
                  'mdl=loadPDB '+pdb,
                  'check mdl',
                  'quit']
            junkf=junkName()  ;  open(junkf,'w').write('\n'.join(cmds))
            lines=os.popen('tleap -f '+junkf+';/bin/rm '+junkf).readlines()
            if 'Close contact' in ''.join(lines): isclash=True
        return {'dump':False,'isclash':isclash}

    def resolveClashes(self,deg=1):
        """twist sidechain around CA-CB bond to remove clashes"""
        angle=0
        while self.checkForClashes()['isclash']:
            self.rotatePDB(deg=deg) ; angle+=deg #; print 'angle=',angle
            if angle>360:
                print 'we could not resolve the clashes'
                return False
        print 'rotated angle:',angle
        return False

    def amberMinimizeInit(self):
        """minimize a PDB structure"""
        input='minimize\n&cntrl\nimin = 1, maxcyc = 500, ntpr = 25,\n'+\
               'igb = 0, ntb = 0, cut = 10.0,\n&end\n'
        inputf=junkName() ; open(inputf,'w').write(input)
        outf=junkName() ; minf=junkName()
        cmd='sander -O -i '+inputf+' -o '+outf+'  -p '+self.initTOP+\
             ' -c '+self.initCRD+' -r '+minf
        os.system(cmd) #execute minimization in vacuum
        cmd='ambpdb -p '+self.initTOP+' < '+minf+' > '+self.initPDB
        os.system(cmd) #overwrite initPDB
        os.system('/bin/rm '+inputf+' '+outf+' '+minf) #rm temp files
        return False
        
    def genMinimizedInit(self):
        """generate a minimized structure in vacuum
        Composed of three steps:
        (1) generate initial structure with tleap
        (2) remove clashes as defined by tleap
        (3) minimize in vacuum"""
        if not self.standarAAA:
            if self.id in ('SPC','H2'):
                self.geninitPDB()     #create initial geometry
                return True
            return False
        
        self.geninitPDB()     #create initial geometry
        self.resolveClashes() #first resolve clashes with tleap
        self.amberMinimizeInit() #minimize resulting structure
        return True

    def setcharge(self,chargelist):
        """enter residue charge"""
        if not self.standarAAA and\
           self.id not in ('SPC','H2'): return False
        self.charge=0
        if self.id in chargelist: self.charge=1
        return True

    def setmultiplicity(self,doubletlist):
        """enter spin multiplicity info"""
        if not self.standarAAA and\
           self.id not in ('SPC','H2'): return False
        self.multiplicity=1
        if self.id in doubletlist: self.multiplicity=2
        return True


    def genoptPDB(self):
        """generate PDB from gaussian optimization run"""
        if not self.standarAAA and\
           self.id not in ('SPC','H2'): return False
        self.hydrd['optPDB']=None
        name=self.id+'_hydrd_opt'
        optfile=joink(self.hydrd['simdir'],name+'.out')
        if os.system('tail -1 '+optfile+' |grep "Normal termination"'):
            sys.stderr.write('incorrect '+optfile+'\n')
            return True
        try:
            if not isthere(optfile):
                sys.stderr.write('Error: no '+optfile)
                return True
            cmd='grep -b "Standard orientation" '+optfile #;Bye(cmd)
            lines=os.popen(cmd).readlines() #;Bye(lines)
            if not lines:
                sys.stderr.write('Error: no Standard orientation')
                return True
            byte=int(lines[-1].split(':')[0]) #byte position
            fin=open(optfile,'r') ; fin.seek(byte)
            line=fin.readline() #;Bye(line)
            if line:
                fin2=open(self.initPDB,'r')  ;  line2=fin2.readline()
                while(line2[0:4]!='ATOM'): line2=fin2.readline()
                buf=''
                for i in range(5): line=fin.readline() #advance four lines
                while '------' not in line:
                    c,a,a,x,y,z=line.split()
                    x=float(x) ; y=float(y) ; z=float(z)
                    buf+=line2[0:30]+'%8.3f%8.3f%8.3f'%(x,y,z)+line2[54:]
                    line=fin.readline()  ;  line2=fin2.readline()
                buf+='TER\nEND\n'
                self.hydrd['optPDB']=joink(self.hydrd['simdir'],name+'.pdb')
                open(self.hydrd['optPDB'],'w').write(buf)
        except:
            sys.stderr.write('Error except\n')
        return True
            

    def g98in(self,job='opt',HD='hydrd',conv=''):
        """generate a g98 input file for optimization or
        frequency calculation for hydrogenated or
        deuterated proteins"""
        if not self.standarAAA and\
           self.id not in ('SPC','H2'): return False
        types=[]  ;  name=self.id+'_'+HD+'_'+job  ;
        title=self.id+'  '+HD+'  '+job  ;  hf='rhf'
        if self.multiplicity>1: hf='uhf'
        buf='%Mem=64MW\n%chk='+name+'.chk\n'
        if job=='opt':
            buf+='#'+hf+'/6-31G* opt SCF=QC\n\n'
            buf+=title+'\n\n'+\
                  `self.charge`+' '+`self.multiplicity`+'\n'
            if self.standarAAA or self.id in ('SPC','H2'):
                for line in open(self.initPDB,'r').readlines():
                    if 'ATOM'==line[0:4]:
                        type=line[13]
                        xyz =line[30:54]
                        buf+=type+'   '+xyz+'\n'
                        types.append(type)            
        elif job=='freq':
            if not self.hydrd['optPDB']:
                sys.stderr.write('no optimed geometry\n')
                self.__dict__[HD]['in'+job]=None
                return True
            buf+='#'+hf+'/6-31G* freq=readisotopes SCF=QC\n\n'+\
                  title+'\n\n'+\
                  `self.charge`+' '+`self.multiplicity`+'\n'
            if self.standarAAA or self.id in ('SPC','H2'):
                for line in open(self.hydrd['optPDB'],'r').readlines():
                    if 'ATOM'==line[0:4]:
                        type=line[13]
                        xyz =line[30:54]
                        buf+=type+'   '+xyz+'\n'
                        types.append(type)            
            buf=buf+'\n298.15 1.0\n'
            for type in types:
                if HD=='deutd' and type=='H': type='D'
                buf+=atomic_weight[type]+'\n'
        buf+='\n'
        ing98file=joink(self.__dict__[HD]['simdir'],name+'.com')
        open(ing98file,'w').write(buf)
        self.__dict__[HD]['in'+job]=ing98file
        return True

    def rmsdInitOpt(self):
        """find RMSD between initPDB and optPDB"""
        if not self.standarAAA: return False
        self.hydrd['rmsdOptInit']=None
        if not self.hydrd['optPDB'] or not isthere(self.hydrd['optPDB']):
            print self.hydrd['rmsdOptInit']
            return False#True
        cmd='rmsd.x '+self.initPDB+' '+self.hydrd['optPDB']+' -b all'
        rmsd=float(os.popen(cmd).readline()) ; print 'rmsd=',rmsd
        self.hydrd['rmsdOptInit']=rmsd
        return True


    def overlapFrequencies(self):
        """for every hydrd frequency, print overlap to deutd ones"""
        if not self.standarAAA and\
           self.id not in ('SPC','H2'): return False
        from gaussian.g98 import gaussian
        import numpy
        for type in ('hydrd','deutd'):
            self.__dict__[type]['outfreq']=None,''
            name=self.id+'_'+type+'_freq.out'
            outf=joink(self.__dict__[type]['simdir'],name)
            if isthere(outf): self.__dict__[type]['outfreq']=outf
        h_g=gaussian(outf=self.hydrd['outfreq']) #gaussian object
        d_g=gaussian(outf=self.deutd['outfreq'])
        ovl=numpy.array([-999.9]*(1+h_g.Job.n)) #store overlaps
        buf='#fh.fn  fd.fn  ovlerlap dot_product fh.f/fd.f  fh.fc/fd.fc\n'
        N=10 #print the five highest overlaps
        for fh in h_g.Job.freqs: #iterate over hydrated frequencies
            for fd in d_g.Job.freqs: ovl[fd.fn]=abs(fh.dot(fd))
            #report indixes of sorted overlaps in decreasing overlap
            sorted=list(numpy.argsort(ovl, axis=0, kind='quicksort'))
            sorted.reverse() ; buf+='\n'
            #print the first N highest overlaps
            for i in sorted[0:N]:
                fd=d_g.Job[i]
                x=0.0 ; y=0.0
                if(fd.f>0.0001): x=fh.f/fd.f
                if(fd.fc>0.0001): y=fh.fc/fd.fc
                props=(fh.fn,fd.fn,ovl[i],fh.overlap(fd),x,y)
                buf+='%03d %03d %4.2f %5.2f %4.2f %4.2f\n'%props
        self.freqoverlap=joink(self.psimdir,'freq_overlap.dat')
        open(self.freqoverlap,'w').write(buf)
        return True

    def printDAAA(self):
        """Print mapping of deuterated AA residue"""
        print '"'+self.id+'" "'+self.deutd['libname']+'"'
        
    def genAAAlib(self):
        """generate topology for deuterated amino acids and their termini
        versions"""
        from amber.amber9 import leap
        from seq.letters import three2one
        pat=re.compile('"H\w*\d*"\s"H\w*"')
        lp=leap()
        libname=''
        if self.standarAAA:
            libname='D'+three2one[ lp.basenameAAA[self.id] ]
            if lp.isstandarAAAname[self.id]: libname+='O'
            else: libname+=self.id[-1]
        else:
            libname='D'+self.id[0:2]
        for p in ('N','','C'):
            if self.id in ('ASH','CYM','GLH','ACE','NME','SPC')\
               and p: continue #exceptions, no termini versions
            libfile=joink(self.deutd['simdir'],p+libname+'.lib')
            if isthere(libfile): os.system('/bin/rm '+libfile)
            buf='source leaprc.ff98\nsaveOff '+p+self.id+\
                 ' '+libfile+'\nquit'
            scriptf=junkName()
            open(scriptf,'w').write(buf)
            os.system('tleap -f '+scriptf+' && /bin/rm '+scriptf)
            buf=''
            for l in  open(libfile,'r').readlines():
                if p+self.id in l:  l=l.replace(p+self.id,p+libname)
                elif pat.search(l): l=l.replace('" "H','" "D')
                buf+=l
            open(libfile,'w').write(buf)
            self.deutd[p+'libname']=p+libname
            self.deutd[p+'libfile']=libfile
            #os.system('emacs -nw -r '+libfile)
        return True

    def genTypeBonds(self,verbose=False):
        """generate type-type and type-type-type bonds"""
        if not self.standarAAA: return False
        csvf=joink(self.psimdir,self.id+'_fc.csv') #;Bye(csvf)
        if not isthere(csvf):
            sys.stderr.write('ERROR no csv file\n')
            return False
        pairs=name2type(self.id) #;Bye(pairs)#name to type mapping
        stretch=[] ; bend=[]
        self.g98comp=csvf #gaussian 98 comparison of force constants
        for line in open(csvf,'r').readlines()[1:]: #first line is comment
            items=line.split()
            if len(items)==8 or len(items)==9: #stretching or bending
                r=items[5-len(items)] #;Bye(items[6-len(items):])#;Bye(r)
                types=[ pairs[name] for name in items[6-len(items):] ]
                if len(types)==2: types.sort() #alphabetic sort
                else:
                    tmp=[types[0],types[-1]] ; tmp.sort()
                    types=[tmp[0],types[1],tmp[-1]]
                x='-'.join(types)+' '+r #;Bye(x)
                if verbose:
                    x+=' '+self.id+' '+line.strip()
                if len(items)==8: stretch.append(x)
                else: bend.append(x)
        return {'dump':False,'stretch':stretch,'bend':bend}
            

    def checkCharge(self):
        """compare charge of hydrogenated and deuterated amino acids"""
        if not self.standarAAA: return False
        buf='source leaprc.ff98D\n'
        buf+='charge '+self.id+'\n'
        buf+='charge '+self.deutd['libname']+'\n'
        buf+='quit\n'
        script=junkName() ; open(script,'w').write(buf)
        tlout=junkName()
        os.system('tleap -f '+script+' > '+tlout)
        pattern=re.compile('Total unperturbed charge:\s*(-*\d+\.\d+)')
        buf=''.join(open(tlout,'r').readlines())
        outl=self.id+' '+self.deutd['libname']
        for match in re.finditer(pattern,buf): outl+=' '+match.group(1)
        print outl
        os.system('/bin/rm '+script+' '+tlout)
        return False

    def overrideAAAtoDDD(self):
        """print tleap command that override the amino acid with its
        deuteraed form"""
        if not self.isstandarAAA: return False
        for p in ('N','','C'):
            if self.id in ('ASH','CYM','GLH','ACE','NME','SPC')\
               and p: continue #exceptions, no termini versions
            print p+self.id+'='+self.deutd[p+'libname']
            print 'set '+p+self.id+'.1 name "'+p+self.id+'"'
        return False
        
#################### end of PEPSICAN  ##################
    
    
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    ojbect ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'psimdir': projects directory to store simulation results
    |_'ssimdir': scratch directory to store simulation results
    |_'standarAAA': boolean flag defines a standar AMBER amino acid
    |_'multiplicity': spin multiplicity
    |_'hydrd': dictonary with info on hydrogenated molecule
    | |_'simdir': directory to store info on hydrogenated molecule
    | |_'inopt':  infile for gaussian optimization
    |_'deutd': dictionary with info on deuterated molecule
      |_'simdir': directory to store info on deuterated molecule
    """
    print globals()['datatree'].__doc__


def doitall(id):

    """a recount of everything we had done on each varset

    Useful if we screw-up a particular dump file, because we can redo
    it all"""

    c=PEPSI(id=id)
    c.genDumpf(cl.repod)

    methods=('createSimdir',)
    for method in methods:
        print method
        getattr(objSL,method)()
    #finally secure the beast!
    c.pickleDump()

def initCollection(id='',idlist=[],idlistf=None):

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
    if id: idlist=[id,]
    if idlistf: idlist=chomp(open(idlistf,'r').readlines())
    if not idlist: idlist=canlist
    remaining=len(idlist)
    for header in idlist:
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

def foo2(*kargs,**kwargs):
    """create seq.dat files
    Example: ... -c id=a0512,switch=hello"""
    cl.cycleOverList('foo',*kargs,**kwargs)

def isstandarAAA(*kargs,**kwargs):
    """define if molecule is an standar amino acid
    Example: ... -c id=ALA"""
    cl.cycleOverList('isstandarAAA',*kargs,**kwargs)
    
def geninitPDB(*kargs,**kwargs):
    """print pdb of the molecule as initial geometry"""
    cl.cycleOverList('geninitPDB',*kargs,**kwargs)

def checkForClashes(*kargs,**kwargs):
    """check a PDB"""
    cl.cycleOverList('checkForClashes',*kargs,**kwargs)
    
def gensimdir(*kargs,**kwargs):
    """mkdir every simdir"""
    cl.cycleOverList('gensimdir',*kargs,**kwargs)

def gendeuterated(*kargs,**kwargs):
    """create deuterated directory"""
    cl.cycleOverList('gendeuterated',*kargs,**kwargs)

def genhydrogenated(*kargs,**kwargs):
    """create hydrogenated directory"""
    cl.cycleOverList('genhydrogenated',*kargs,**kwargs)

def setcharge(*kargs,**kwargs):
    """enter residue charge
    Example: -a setcharge -b doublet.list"""
    kargs=[chomp(open(kargs[0],'r').readlines()),]
    cl.cycleOverList('setcharge',*kargs,**kwargs)

def setmultiplicity(*kargs,**kwargs):
    """enter double spin multiplicity
    Example: -a setmultiplicity -b doublet.list"""
    kargs=[chomp(open(kargs[0],'r').readlines()),]
    cl.cycleOverList('setmultiplicity',*kargs,**kwargs)

def g98in(*kargs,**kwargs):
    """generate a g98 input file for optimization or frequency
    calculation for hydrogenated or deuterated proteins"""
    cl.cycleOverList('g98in',*kargs,**kwargs)

def g98run(id="",job='opt',HD='hydrd'):
    import popen2
    idlist=cl.canlist
    if id: idlist=[id,]
    running=[]  ; n=0  ;  nproc=3  ;  lapse=600  ;  left=len(idlist)
    for id in idlist:
        c=cl.loadCAN(id)  ;  print left,c.id
        if not (c.standarAAA or c.id in ('SPC','H2')) : continue
        name=c.id+'_'+HD+'_'+job
        #check if input file exists
        inf=c.__dict__[HD]['in'+job]
        if not inf:
            sys.stderr.write('no input file\n')
            continue #no input, go to next id
        #check if out file exists with error
        outf=joink(c.__dict__[HD]['simdir'],name+'.out')
        if isthere(outf):
            if os.system('tail -1 '+outf+' |grep "Normal termination"'):
                os.system('/bin/rm '+outf)
                chkf=joink('./',name+'.chk')#  ;  print outf,chkf
                if isthere(chkf): os.system('/bin/rm '+chkf)
            else:
                sys.stdout.write('completed\n')
                continue #we successfully did this job
        #cmd='sleep '+`int(20*random())` ; print cmd #for debugging
        cmd= joink(os.environ['HOME'],'bin','g98.csh')+'  '+inf
        child=popen2.Popen3(cmd)  ;  sleep(60)  #;  print cmd
        left-=1  ;  running.append(child)
        while(float(os.popen('uptime').readline().split()[-3][0:-1]))>nproc:
            sleep(lapse)

def genoptPDB(*kargs,**kwargs):
    """generate PDB from gaussian optimization run"""
    cl.cycleOverList('genoptPDB',*kargs,**kwargs)

def rotatePDB(*kargs,**kwargs):
    """rotate around the CA-CB bond"""
    cl.cycleOverList('rotatePDB',*kargs,**kwargs)

def resolveClashes(*kargs,**kwargs):
    """twist sidechain around CA-CB bond to remove clashes"""
    cl.cycleOverList('resolveClashes',*kargs,**kwargs)

def amberMinimizeInit(*kargs,**kwargs):
    """minimize structur in vacuum"""
    cl.cycleOverList('amberMinimizeInit',*kargs,**kwargs)

def genMinimizedInit(*kargs,**kwargs):
    """generate a minimized structure in vacuum
    Composed of three steps:
    (1) generate initial structure with tleap
    (2) remove clashes as defined by tleap
    (3) minimize in vacuum"""
    cl.cycleOverList('genMinimizedInit',*kargs,**kwargs)

def rmsdInitOpt(*kargs,**kwargs):
    """find RMSD between initPDB and optPDB"""
    cl.cycleOverList('rmsdInitOpt',*kargs,**kwargs)

def overlapFrequencies(*kargs,**kwargs):
    """for every hydrd frequency, print overlap to deutd ones"""
    cl.cycleOverList('overlapFrequencies',*kargs,**kwargs)

def genAAAlib(*kargs,**kwargs):
    """generate topology for deuterated amino acid"""
    cl.cycleOverList('genAAAlib',*kargs,**kwargs)

def printDAAA(*kargs,**kwargs):
    """Print mapping of deuterated AA residue"""
    cl.cycleOverList('printDAAA',*kargs,**kwargs)

def genTypeBonds(*kargs,**kwargs):
    """generate type-type and type-type-type bonds"""
    stretch=[] ; bend=[]
    for (id,x) in cl.cycleOverList('genTypeBonds',*kargs,**kwargs).items():
        stretch+=x['stretch'] ; bend+=x['bend']
    stretch.sort() ; bend.sort()
    print '\n'.join(stretch)
    print '\n'.join(bend)

def genParm98Bonds():
    """for every stretch and bend bond of
    analysis/genTypeBonds_trim.dat, generate appropriate line for the
    deuterated version by looking for the line in parm98.dat"""
    blen={'stretch':5,'bend':8}
    ffmt={'stretch':'%5.1f','bend':'%4.1f'} #force constant format
    fpos={'stretch':(7,12),'bend':(12,16)} #force constant position
    parm98f=joink(os.environ['AMBERHOME'],'dat/leap/parm/parm98.dat')
    parm98=open(parm98f,'r').readlines()
    bondsf=joink(pjd,'analysis/genTypeBonds_trim.dat')
    buf=''
    for line_bondsf in open(bondsf,'r').readlines():
        bondMissing=True
        bond_type='stretch'
        if len(line_bondsf.split('-'))==3: bond_type='bend'
        bond=line_bondsf[0:blen[bond_type]]
        r=float(line_bondsf[blen[bond_type]:]) #;Bye(r)
        for line_parm98 in parm98: #search bond in parm98
            if len(line_parm98)<2: continue
            putatitebond=line_parm98[0:1+blen[bond_type]]
            if bond+' '==putatitebond: #bond found->change bond and ff
                bondMissing=False
                ffH=line_parm98[fpos[bond_type][0]:fpos[bond_type][1]]
                ffD=ffmt[bond_type]%( float(ffH)*float(r) )
                bond=bond.replace('H','D')
                bond=bond.replace('SD','SH') #take it back
                bond=bond.replace('OD','OH') #take it back
                line=bond+\
                      line_parm98[ blen[bond_type]: fpos[bond_type][0]]+\
                      ffD+\
                      line_parm98[fpos[bond_type][1]:]
                buf+=line
                break
        if bondMissing: sys.stderr.write(bond+'\n')
    print buf
                
def deuterateParm98(task='torsion'):
    """change H to D in selected entries of parm98.dat"""
    parm98f=joink(os.environ['AMBERHOME'],'dat/leap/parm/parm98.dat')
    parm98=open(parm98f,'r').readlines()
    buf='' #store deuterated lines of parm98.dat
    if task=='torsion':
        """change H to D in torsion entries """
        unit='\w[\w\s\*]-' #atom type
        regexp=unit*3+'\w[\w\s\*]' #torsion is four of them
        pattern=re.compile('('+regexp+')') #pattern to match
        for line in parm98:
            match=pattern.match(line)
            if match:
                mH=match.group(1)
                if 'H' in (mH[0],mH[3],mH[6],mH[9]):
                    mD=mH.replace('H','D')
                    buf+=line.replace(mH,mD)
    elif task=='VdW':
        """output VdW parameters for the deuterated versions"""
        parm98p=open(parm98f,'r') ; line=parm98p.readline()
        while line.find('MOD4')<0:
            line=parm98p.readline()
        buf='' ; line=parm98p.readline() #;print line
        while len(line)>2:
            if line[2]=='H':
                d='%6.4f'%(0.995*float(line[14:20]))
                buf+=line[0:2]+'D'+line[3:14]+d+line[20:]
            line=parm98p.readline()
    elif task=='HP':
        print 'DP-CT-C     50.0      109.50    AA zwitterion  changed based on NMA nmodes'

def deuteratePDB(pdbf):
    """substitute of H for D in a PDB file"""
    from amber.amber9 import leap
    from seq.letters import three2one
    lp=leap()
    buf=''
    for l in open(pdbf,'r').readlines():
        if buf[0:5]=='ATOM ':
            if buf[13]=='H':l=l[0:13]+'H'+l[14:]
            elif buf[13:15]=='SH': l=l[0:13]+'SD'+l[16:]
            elif buf[13:15]=='OH': l=l[0:13]+'OD'+l[16:]
            AAA='D'+three2one[ lp.basenameAAA[ l[17:20] ] ]
            l=l[0:17]+AAA+l[20:]
        buf+=l
    print buf        

def checkCharge(*kargs,**kwargs):
    """compare charge of hydrogenated and deuterated amino acids"""
    cl.cycleOverList('checkCharge',*kargs,**kwargs)

def overrideAAAtoDDD(*kargs,**kwargs):
    """print tleap command that override the amino acid with its
    deuteraed form"""
    cl.cycleOverList('overrideAAAtoDDD',*kargs,**kwargs)

def help(*kargs):

    """list available services"""
    
    services=['datatree','doitall','updateCollection','initCollection','createSimdir',
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
