import os,sys,re,numpy
from utilities.small_utilities import Bye
import math

class atom:
    def __init__(self,*kargs,**kwargs):
        self.type=None              #atom type (H,C,N,..)
        self.m=None                 #mass (AMU)
        self.an=None                #atomic number
        self.q=None                 #charge
        self.s=None                 #spin
        self.n=None                 #atom number
        self.xyz=numpy.zeros(3)     #atomic coordinates
        self.__dict__.update(kwargs)

#atom templates
attpl={'H':atom(type='H',an=1),
       'N':atom(type='N',an=7),
       'O':atom(type='O',an=8),
       'C':atom(type='C',an=12),
       'S':atom(type='S',an=16)
       }


class Conf:
    """system configuration"""
    def __init__(self):
        self.charge=None  #charge
        self.mult=None    #sping multiplicity
        self.ats=[]       #list of atoms
        self.n=0

    def __getitem__(self,key):
        return self.ats[key-1]  #gaussian indexing of atoms begins with 1
        
    def readout(self,outf):
        """read system configuration"""
        outp=open(outf,'r') ; l=outp.readline() #;Bye(outf)
        while l[0:9] != ' Charge =': l=outp.readline()
        pat=re.compile(' Charge =\s+(-*\d+)\s+Multiplicity =\s+(\d+)')
        match=pat.search(l)
        self.charge=int(match.group(1))
        self.mult=int(match.group(2))
        n=0
        ats=[] #list of atoms
        items=outp.readline().split()
        while items:
            n+=1 ; type=items[0]
            at=attpl[type] #template atom
            x=float(items[1]) ; y=float(items[2]) ; z=float(items[3])
            at.xyz= numpy.array( [x,y,z] ) #insert coordinates
            at.n=n
            ats.append(at)
            items=outp.readline().split()
        self.ats=ats
        self.n=n

class frequency:
    def __init__(self):
        self.f=None       #frequency (cm****-1)
        self.orb=None     #orbital
        self.rm=None      #reduced mass (AMU)
        self.fc=None      #force constant (mDyne/A)
        self.IR=None      #IR Inten
        self.raman=None   #Raman Activity
        self.depolar=None #Depolar
        self.fn=None      #frequency number
        self.atd=numpy.array([])  #array of atomic displacements
        self.n=0          #number of atomic displacements (number atoms)

    def __repr__(self):
        buf=''
        buf+='                           %3d\n'%self.fn
        buf+='                           %s\n'%self.orb
        buf+='Frequencies --%11.4f\n'%self.f
        buf+='Red. masses --%11.4f\n'%self.rm
        buf+='Frc consts  --%11.4f\n'%self.fc
        buf+='IR Inten    --%11.4f\n'%self.IR
        buf+='Raman Activ --%11.4f\n'%self.raman
        buf+='Depolar     --%11.4f\n'%self.depolar
        buf+='Atom AN      X      Y      Z\n'
        i=0
        for atd in self.atd:
            i+=1
            buf+='%4d%4d  %7.2f%7.2f%7.2f\n'%(i,0,atd[0],atd[1],atd[2])
        return buf

    def __getitem__(self,key):
        return self.atd[key-1]  #gaussian indexing of atoms begins with 1
    
    def dot(self,right):
        """dot product of two frequencies is dot product of their
        atomic displacements"""
        if self.n != right.n:
            sys.stderr.write('unequal number of atomic displacements\n')
            return None
        return numpy.dot(self.atd.reshape(3*self.n),\
                         right.atd.reshape(3*right.n))

    def overlap(self,right):
        """overlap=sqrt( sum( dot(self.atd[i],right.atd[i])^2  ) )"""
        return math.sqrt(numpy.sum(\
            numpy.dot(self.atd,right.atd.T).diagonal()**2))
        
subclasses_types={'Job':('freq',),
                  'Method':(),
                  'Basis':()
                  }

class Job:
    def __init__(self):
        self.type=None
        self.parent=None      #reference to parent gaussian object

class freq(Job):
    def __init__(self):
        self.type='freq'
        self.freqs=[]         #list   of frequencies
        self.n=0              #number of frequencies
        
    def __getitem__(self,key):
        return self.freqs[key-1]
    
    def read_triad(self,outp):
        """read three frequencies"""
        triad=[frequency(),frequency(),frequency()]
        for prop in ('fn','orb','f','rm','fc','IR','raman','depolar'):
            items=outp.readline().split()
            if not items: return None
            for i in range(3):
                x=items[-3+i]
                if prop in ('fn',):x=int(x)
                elif prop in ('f','rm','fc','IR','raman','depolar'):
                    x=float(x)
                triad[i].__dict__[prop]=x
        outp.readline() ; n=0 ; items=outp.readline().split()
        lists=[ [], [], [] ]
        while len(items)==11:
            for i in range(3):
                x=float(items[i*3+2])
                y=float(items[i*3+3])
                z=float(items[i*3+4]) 
                lists[i].append([x,y,z])
            position=outp.tell() ; items=outp.readline().split() ; n+=1
        outp.seek(position)
        for i in range(3):
            triad[i].atd=numpy.array(lists[i])
            triad[i].n=n
        return triad
            
    def read_frequencies(self,outf):
        outp=open(outf,'r') ; l=outp.readline()
        while 'reduced masses (AMU), force' not in l: l=outp.readline()
        triad=self.read_triad(outp)
        while triad:
            self.freqs+=triad
            triad=self.read_triad(outp)
        self.n=len(self.freqs)
            
    def readout(self,outf):
        """read related info to frequencies in the output file"""
        self.read_frequencies(outf)  #read frequencies
        

class Method:
    def __init__(self):
        self.type=None
        self.parent=None   #reference to parent gaussian object
        
class Basis:
    def __init__(self):
        self.type=None
        self.parent=None   #reference to parent gaussian object
        
class out_g:
    """gaussian output"""
    def __init__(self,parent=None):
        self.outf=None     #gaussian output file
        self.parent=parent #link to parent gaussian object
        
    def present(self,parentprops ):
        """read a list of properties to present to parent gaussian obj

        properties are  names of base classes like Job, Method, and
        Basis"""
        buf=''.join(open(self.outf,'r').readlines())
        patterns={'Method':'#(\w+)/', 'Basis':'#\w+/(\S+)\s+',
                  'Job':'#\w+/\S+\s+(\w+)'}
        if type(parentprops) not in (tuple,list):
            parentprops=(parentprops,)
        for parentprop in parentprops:
            if parentprop not in patterns.keys(): continue
            match=re.compile(patterns[parentprop]).search(buf)
            if not match: continue
            prop=match.group(1).lower() #e.g: parentprop=Job, prop=Opt
            if prop not in subclasses_types[parentprop]:
                prop=parentprop
            self.parent.__dict__[parentprop]=globals()[prop]()
            self.parent.__dict__[parentprop].parent=self.parent

    def readout(self,outf):
        """read an output file"""
        try:
            outp=open(outf,'r')
            self.outf=outf
            self.parent.Conf.readout(self.outf)
            self.present( ('Job','Method','Basis'))
            self.parent.Job.readout(self.outf) #read particulars of Job
        except IOError:
            sys.stderr.write('No output file!\n')
            
        
class gaussian:
    def __init__(self,outf=None):
        self.Conf=Conf() #system configuration
        self.Job=None    #Job type
        self.Method=None #method typ
        self.Basis=None  #basis set
        self.out=out_g(parent=self) #output object
        if outf: self.out.readout(outf)

    def read_outf(self,outf):
        """read output file"""
        self.out.readout(outf) #pass the reading task to out object
        
