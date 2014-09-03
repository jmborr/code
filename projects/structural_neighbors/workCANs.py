#!/usr/bin/python

import os,sys,re
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import inpHand
from jobs.job import pastry

"""manage poly-VAL peptides"""

#global variables
joink=os.path.join
isthere=os.path.exists
ocd='/gpfs1/scratch/adrian/protsup/netdist0.4' #results from the real sequences
scop='/gpfs1/scratch/adrian/protsup/scop'
pcd=joink(codedir,'projects','structural_neighbors')
scd=joink(scratchdir,'structural_neighbors','out')
cad=joink(scratchdir,'structural_neighbors','pdb') #directory with CA PDB's for natives
pdb_jul07=joink(scratchdir,'structural_neighbors','pdb_jul07')
collf=joink(scd,'varsets.can')
canlistf=joink(scd,'list')
repod=joink(scd,'repod')
alignx=codedir+'/f77/seq/align'
pdbrootds=('/gpfs1/scratch/adrian/dat/pdb_jul07','/gpfs1/scratch/adrian/dat/pdb_mar07',
       '/gpfs1/scratch/adrian/dat/pdb200')
inputds=('/gpfs1/scratch/adrian/dat/pdb_jul07/input','/gpfs1/scratch/adrian/dat/pdb_mar07/input',
       '/gpfs1/scratch/adrian/dat/pdb200/input')
pdbds=('/gpfs1/scratch/adrian/dat/pdb_jul07/CA','/gpfs1/scratch/adrian/dat/pdb_mar07/CA',
       '/gpfs1/scratch/adrian/dat/pdb200/CA') #different pdb releases
secmaster=joink(cad,'secmaster')
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


    def addSeq(self):

        """add sequence of native PDB protein"""

        from seq.fastaManager import importFastaEntry
        
        if 'pdb300' not in self.__dict__: self.pdb300={}
        for inputd in inputds:
            inputf=joink(inputd,self.id)
            if isthere(inputf):
                self.pdb300['seq']=importFastaEntry(open(inputf,'r'))[1]
                break
        return True

        
    def extractSS(self,ss=''):

        """access protein associated with PDB id, then extract its secondary structure

        Extract sec. str. according to two methods: linker42sum and
        pdb2ss.py.
        linker42sum is not very good with helixes (see 153l_ for
        instance), and pdb2ss is not very good with strands (see 117eA
        for instance). Thus, we mix them"""

        from secondary.pdb2ss import outputSec
        from secondary.dssp import dssp2ss
        if 'secondary' not in self.__dict__: self.secondary={}
        if ss:
            self.secondary['HE-']=ss
        else:
            pin=open(secmaster,'r')
            header=pin.readline()[1:6]
            while header:
                if header==self.id:
                    self.secondary['HE-']=pin.readline().strip()
                    break
                header=pin.readline()
        pdbf=joink(cad,self.id)
        self.secondary['HEC']=outputSec(pdbf).replace('T','C')

        mix=''
        for i in range( len(self.secondary['HEC']) ):
            x=self.secondary['HEC'][i]
            if self.secondary['HE-'][i]=='E': x='E' #overwrite strand assignment
            mix+=x
        self.secondary['mix']=mix


    def extracSSfromDSSP(self,sspt,indexl):
        """extract secondary structure from DSSP asignments"""
        from seq.fastaManager import importFastaEntry
        from secondary.dssp import dssp2ss,paa
        id=self.id[0:4]+':'  ;  id=id.upper()
        if self.id[4]=='_': id+='A'
        else: id+=self.id[4]
        if id not in indexl.keys(): self.secondary['dssp']=self.secondary['mix']
        else:
            index=indexl[id]  ;
            sspt.seek(index,0)  ;
            buf=''
            for x in importFastaEntry(sspt,paa=paa)[1]:
                if x in dssp2ss.keys(): buf+=dssp2ss[x]
                else: buf+='C'
            buf=buf.replace('HHHHHCHHHH','HHHHCCHHHH') #fix sudden orient change in two helixes
            self.secondary['dssp']=buf
            if len(buf)!=self.L: self.secondary['dssp']=self.secondary['mix']
        return True
    

    def megamixSS(self):
        """We prepare a megamix assignment for secondary structure, whereby we use helical
        assignments from dssp and strand assignments from mix"""
        buf=''
        for i in range(self.L):
            if self.secondary['dssp'][i]=='H': buf+='H'
            elif self.secondary['dssp'][i]=='E': buf+='E'
            elif self.secondary['mix'][i]=='E': buf+='E'
            else: buf+='C'
        self.secondary['megamix']=buf
        return True

        
    def genSeqDatFile(self,switch=''):

        """generate seq.dat file"""

        ss2n={'H':2,'E':4,'C':1}
        ss=self.secondary['mix']  ;  outdir=self.simdir
        if switch=='megamix':
            ss=self.secondary['megamix']
            if not 'megamix' in self.__dict__.keys():
                self.megamix={'simdir':joink(self.simdir,'megamix')}
            outdir=self.megamix['simdir']
            if not isthere(outdir): os.system('/bin/mkdir -p '+outdir)
        n=1
        lines=[]
        for x in ss:
            lines.append( '%3d VAL %s 10'%(n,ss2n[x]) )
            n+=1
        open( joink(outdir,'seq.dat'),'w' ).write( '\n'.join(lines) )
        self.megamix['seqdat']=joink(self.megamix['simdir'],'seq.dat')
        return True


    def genrmsinp(self):

        """generate rmsinp file and sequence length"""

        self.L=len( self.secondary['HE-'] )
        open(joink(self.simdir,'rmsinp'),'w').write('1 '+`self.L`+'\n'+`self.L`+'\n'+self.id)
        return False


    def genChainDatFile(self,switch=''):
        """generate chain.dat file

        We only output coordinates for helixes. Strands are treated
        on-lattice"""
        seqdat=joink(self.simdir,'seq.dat')
        chaindat=joink(self.simdir,'chain.dat')
        initpdb=joink(self.simdir,'init.pdb')
        if switch=='megamix':
            seqdat=self.megamix['seqdat']
            chaindat=joink(self.megamix['simdir'],'chain.dat')
            self.megamix['chaindat']=chaindat
            initpdb=joink(self.megamix['simdir'],'init.pdb')
        cmd=joink(pcd,'genChainDat.py '+
                  ' -a '+seqdat                     +
                  ' -b '+joink(pcd,'betahelix.pdb') +
                  ' -c '+chaindat                   +
                  ' -d '+initpdb                    +
                  ' -f no'
                  )
        if os.system(cmd):
            sys.stderr.write('ERROR(genChainDatFile) did not create chain.dat file for'+\
                             self.id+'\n')
        return True
    

    def checkChainDatFile(self):

        """check existence of chain.dat file. Try to produce one if non-existent"""

        cdf=joink(self.simdir,'chain.dat')
        if not isthere(cdf):
            sys.stderr.write('chain.dat missing for '+self.id+'... ')
            self.genChainDatFile()
        return False
            

    def secondary_type(self):

        """gather secondary structure content and estimate secondary class"""
        pt=self.secondary
        secdat=pt['mix']
        fa=float( len(re.compile('H').findall(secdat)) )/self.L
        fb=float( len(re.compile('E').findall(secdat)) )/self.L
        fc=1-fa-fb
        if fa>0.4 and fb<0.1:   type='a'
        elif fb>0.4 and fa<0.1: type='b'
        elif fa>0.2 and fb>0.2: type='c'
        else: type='d'
        pt['content']={'H':fa,'E':fb,'C':fc} #;print pt['content']
        pt['type']=type #;print type
        return True


    def addL(self):

        """add sequence length"""

        self.L=len(self.secondary['mix'])
    

    def submitGenericTasser(self,joblist='',protocol='initial',ibatch=0):

        """submit generic TASSER+SPICKER+AMBERMIN jobs"""

        from jobs.job import job,genJobList
        from tasser.tasser_utilities import calc_wall_time

        if not joblist: joblist=genJobList()
        wt=calc_wall_time(L=self.L) #estimated TASSER wall time, in days
        wt='%5.2lf'%( float(wt)+0.50) #add 1/2 day for spicker plus amber minimizations

        if protocol=='initial':
            outdir=self.simdir
            args=' '+self.id+' casB_gene4 '+outdir
            jobname='g'+self.id
        elif protocol=='second neighbors':
            outdir=joink(self.simdir,'eh6b')
            args=' '+self.id+' casB_gene5 '+outdir
            jobname='eh6b_'+self.id
        elif protocol=='secneig arfavor:1':
            outdir=joink(self.simdir,'1.0/eh6b')
            args=' '+self.id+' casB_gene5 '+outdir
            jobname='eh6b_'+self.id
        elif protocol=='megamix':
            outdir=self.megamix['simdir']
            args=' '+self.id+' casB_gene5 '+outdir+' megamix'
            jobname='meg_'+self.id
        Job=job(name='generic tasser',exe='submitGenericTasser.sh',exed=pcd,args=args)
        Job.qsub(jobname,outdir,wallt=wt,mem_limit='990',ddisk='400',joblist=joblist)
        ibatch=Job.dormant(ibatch)
        return ibatch,False
    

    def initMininized(self):

        """init 'minimized' dictionary"""

        from utilities.small_utilities import unTARme,chomp

        self.minimized={'outfile':'',}
        p=self.minimized
        outfile=joink(self.simdir,'out.minimized.tbz2')
        if not isthere(outfile): return True #just initialize self.minimized to previous line
        p['outfile']=outfile
        wd=unTARme(p['outfile']) #untar PDB structures
        currd=os.getcwd()
        os.chdir(wd)
        p['list']=chomp(os.popen('ls -1').readlines())
        os.chdir(currd)
        os.system('/bin/rm -r '+wd)
        for str in p['list']: p[str]={}
        return True
    

    def createTasserSpickerStructuresDicts(self,switch=''):

        """initialize 'tasser', 'spicker', and 'structures' dictionaries"""

        from utilities.small_utilities import unTARme
        if switch=='megamix':
            ptb=self.megamix
            for dict in ('tasser','spicker'):
                ptb[dict]={}  ;  ptb[dict]['outfile']=''
                outfile=joink(ptb['simdir'],'out.'+dict+'.tbz2')
                if isthere(outfile):
                    ptb[dict]['outfile']=outfile
                else:
                    sys.stderr.write('ERROR: no outfile for '+self.id+' '+dict+' simulation\n')
        elif switch=='secneig arfavor:1':
            self.arfavor1={}  ;  ptb=self.arfavor1
            ptb['simdir']=joink(self.simdir,'1.0/eh6b')
            for dict in ('tasser','spicker'):
                ptb[dict]={}
                ptb[dict]['outfile']=''
                outfile=joink(ptb['simdir'],'out.'+dict+'.tbz2')
                if isthere(outfile):
                    ptb[dict]['outfile']=outfile
                else:
                    sys.stderr.write('ERROR: no outfile for '+self.id+' '+dict+' simulation\n')
        elif switch=='second neighbors':
            self.secneig={}
            ptb=self.secneig
            ptb['simdir']=joink(self.simdir,'eh6b')
            for dict in ('tasser','spicker'):
                ptb[dict]={}
                ptb[dict]['outfile']=''
                outfile=joink(ptb['simdir'],'out.'+dict+'.tbz2')
                if isthere(outfile):
                    ptb[dict]['outfile']=outfile
                else:
                    sys.stderr.write('ERROR: no outfile for '+self.id+' '+dict+' simulation\n')
        return True
    

    def gatherStructures(self,switch=''):

        """gather structures from the TASSER simulations"""

        from utilities.small_utilities import unTARme,chomp

        if switch=='second neighbors':
            if 'secneig' not in self.__dict__.keys():
                self.secneig={'simdir':joink(self.simdir,'eh6b'),'outfile':'','list':[]}
            p0=self.secneig
            outfile=joink(p0['spicker']['outfile'])
            if isthere(outfile):
                wd=unTARme(outfile,fileL=['closc1.pdb',]) #untar structure to temporary directory
                if not wd: return false
                f0=joink(wd,'closc1.pdb')
                if not isthere(f0):
                    sys.stderr.write('ERROR: no closc1.pdb in out.spicker.tbz2\n')
                    pastry('/bin/rm -r '+wd)
                    return False
                f1=joink(scd,'structures_CA_secneig',self.id+'.fake')
                cmd='grep " CA " '+f0+'>'+f1 #;Bye(cmd)
                pastry(cmd)
                pastry('/bin/rm -r '+wd)
        elif switch=='secneig arfavor:1':
            if 'arfavor1' not in self.__dict__.keys(): return False
            p0=self.arfavor1
            outfile=joink(p0['spicker']['outfile'])
            if isthere(outfile):
                wd=unTARme(outfile,fileL=['closc1.pdb',]) #untar structure to temporary directory
                if not wd: return False
                f0=joink(wd,'closc1.pdb')
                if not isthere(f0):
                    sys.stderr.write('ERROR: no closc1.pdb in out.spicker.tbz2\n')
                    pastry('/bin/rm -r '+wd)
                    return False
                f1=joink(scd,'structures_CA_arfavor1',self.id+'.fake')
                cmd='grep " CA " '+f0+'>'+f1 #;Bye(cmd)
                pastry(cmd)
                pastry('/bin/rm -r '+wd)
        elif switch=='megamix':
            if 'megamix' not in self.__dict__.keys(): return False
            p0=self.megamix
            outfile=joink(p0['spicker']['outfile'])
            if isthere(outfile):
                wd=unTARme(outfile,fileL=['closc1.pdb',]) #untar structure to temporary directory
                if not wd: return False
                f0=joink(wd,'closc1.pdb')
                if not isthere(f0):
                    sys.stderr.write('ERROR: no closc1.pdb in out.spicker.tbz2\n')
                    pastry('/bin/rm -r '+wd)
                    return False
                f1=joink(scd,'structures_CA_megamix',self.id+'.fake')
                pastry('grep " CA " '+f0+'>'+f1)
                pastry('/bin/rm -r '+wd)                        
        else:
            strd=joink(scd,'structures_CA') #structure repository directory
            if  'minimized' not in self.__dict__: return False
            p0=self.minimized
            if not isthere(p0['outfile']): return False #no minimized structures
            wd=unTARme( p0['outfile'] ) #untar structures to temporary directory
            f0=joink(wd,p0['list'][0])   #full atom structure file, origin
            f1=joink(strd,p0['list'][0]) #CA-only structure file, destination
            cmd='grep " CA " '+f0+' > '+f1 #;Bye(cmd)
            os.system(cmd)
            os.system('/bin/rm -r '+wd)
        return False #there is nothing to update the pickle file in this method
    

    def nneighbor(self,switch=''):

        """enter nneighbor info"""

        if switch=='original to original':
            if 'original' not in self.__dict__.keys(): return False
            pt=self.original
            self.original['nneighbor']['to original']={'outfile':'','n':0}
            pt1=self.original['nneighbor']['to original']
            outfile=joink(scd,'tmX2allgene2',self.id+'.1.pdb.dat')
            if isthere(outfile):
                pt1['outfile']=outfile
                n=int(open(outfile,'r').readline().strip())  #number of neighbors
                pt1['n']=n
                return True
        elif switch=='secneig to secneig':
            if 'nneighbor' not in self.secneig.keys(): return False
            self.secneig['nneighbor']['to secneig']={'outfile':'','n':0}
            pt1=self.secneig['nneighbor']['to secneig']
            outfile=joink(scd,'tmX2all_secneig',self.id+'.dat')
            if isthere(outfile):
                pt1['outfile']=outfile
                n=int(open(outfile,'r').readline().strip())  #number of neighbors
                pt1['n']=n
            return True
        elif switch=='arfavor1 to arfavor1':
            if 'nneighbor' not in self.arfavor1.keys(): self.arfavor1['nneighbor']={}
            self.arfavor1['nneighbor']['to arfavor1']={'outfile':'','n':0}
            pt1=self.arfavor1['nneighbor']['to arfavor1']
            outfile=joink(scd,'tmX2all_arfavor1',self.id+'.fake.dat')
            if isthere(outfile):
                pt1['outfile']=outfile
                n=int(open(outfile,'r').readline().strip())  #number of neighbors
                pt1['n']=n
            return True
        elif switch=='megamix to megamix':
            if 'nneighbor' not in self.megamix.keys(): self.megamix['nneighbor']={}
            self.megamix['nneighbor']['to megamix']={'outfile':'','n':0}
            pt1=self.megamix['nneighbor']['to megamix']
            outfile=joink(scd,'tmX2all_megamix',self.id+'.fake.dat')
            if isthere(outfile):
                pt1['outfile']=outfile
                n=int(open(outfile,'r').readline().strip())  #number of neighbors
                pt1['n']=n
            return True
        elif switch=='snsd to snsd': #second neighbor single domain to itself
            if 'pdp' in self.secneig.keys() and self.secneig['pdp']['n']==1:
                if 'snsd' not in self.__dict__: self.snsd={}
                if 'nneighbor' not in self.snsd.keys(): self.snsd['nneighbor']={}
                self.snsd['nneighbor']['to snsd']={'outfile':'','n':0}
                pt1=self.snsd['nneighbor']['to snsd']
                outfile=joink(scd,'tmX2all_snsd_snsd',self.id+'.fake.dat')
                if isthere(outfile):
                    pt1['outfile']=outfile
                    n=int(open(outfile,'r').readline().strip())  #number of neighbors
                    pt1['n']=n
                return True
            return False
        elif switch=='secneig to pdb_jul07':
            if 'nneighbor' not in self.secneig.keys(): return False
            self.secneig['nneighbor']['to pdb_jul07']={'outfile':'','n':0}
            pt1=self.secneig['nneighbor']['to pdb_jul07']
            outfile=joink(scd,'tmX2all_snpdb_jul07',self.id+'.fake.dat')
            if isthere(outfile):
                pt1['outfile']=outfile
                n=int(open(outfile,'r').readline().strip())  #number of neighbors
                pt1['n']=n
                return True
        elif switch=='pdb300 to pdb300':
            if 'nneighbor' not in self.pdb300.keys(): self.pdb300={ 'nneighbor':{} }
            self.pdb300['nneighbor']['to pdb300']={'neighbors':[],'n':0}
            ptp=self.pdb300['nneighbor']['to pdb300']
            pt=open(joink(ocd,self.id+'.list'),'r')  ;  line=pt.readline()
            nneig=0 #because itself is allways included in the .list file
            while line:
                id,n,tm=line.split()
                if n!='1': break
                if id in cl.canlist and id!=self.id:
                    nneig+=1 #only keep first neighbors below 300 residues
                    ptp['neighbors'].append(id)
                line=pt.readline()
            ptp['n']=nneig
            return True
        elif switch=='second neighbors to PDB':
            pt=self.secneig
            if 'nneighbor' not in pt.keys():pt['nneighbor']={}
            pt['nneighbor']['to PDB']={'outfile':'','n':0}
            pt1=self.secneig['nneighbor']['to PDB']
            outfile=joink(scd,'tmX2all_snPDB',self.id+'.fake.dat')
            if isthere(outfile):
                pt1['outfile']=outfile
                n=int(open(outfile,'r').readline().strip())  #number of neighbors
                pt1['n']=n
                return True
        elif switch=='PDB300 to second neighbors':
            if 'pdb300' not in self.__dict__: self.pdb300={}
            pt=self.pdb300
            if 'nneighbor' not in pt.keys():pt['nneighbor']={}
            pt['nneighbor']['to secneig']={'outfile':'','n':0}
            pt1=self.pdb300['nneighbor']['to secneig']
            outfile=joink(scd,'tmX2all_snPDB',self.id+'.dat')
            if isthere(outfile):
                pt1['outfile']=outfile
                n=int(open(outfile,'r').readline().strip())  #number of neighbors
                pt1['n']=n
                return True
        elif switch=='original to PDB':
            if 'original' not in self.__dict__.keys(): self.__dict__['original']={}
            pt=self.original
            if 'nneighbor' not in pt.keys():pt['nneighbor']={}
            pt['nneighbor']['to PDB']={'outfile':'','n':0}
            pt1=self.original['nneighbor']['to PDB']
            outfile=joink(scd,'tmX2allgene2pdb',self.id+'.1.pdb.dat')
            if isthere(outfile):
                pt1['outfile']=outfile
                n=int(open(outfile,'r').readline().strip())  #number of neighbors
                pt1['n']=n
                if 'original' in self.secneig['nneighbor'].keys():
                    del self.secneig['nneighbor']['original']
                return True            
        return False


    def addPDBnneigh(self):
    
        """ add number of neighbors for real PDB structure to other real PDB structures"""

        pt=open(joink(ocd,self.id+'.list'),'r')  ;  line=pt.readline()
        nneig=-1 #because itself is allways included in the .list file
        while line:
            id,n,tm=line.split()
            if n!='1': break
            if id in cl.canlist: nneig+=1
            line=pt.readline()
        self.realpdbnneig=nneig
        return True


    def addNativeFile(self,switch=''):

        """add native file path"""

        if switch not in self.__dict__.keys(): return False
        if switch=='pdb300':
            self.__dict__['pdb300']['CA']=joink(cad,self.id)
        elif switch=='original':
            self.__dict__['original']['CA']=''
            CA=joink(scd,'structures_CA_pdbgene2',self.id+'.1.pdb')
            if isthere(CA): self.__dict__['original']['CA']=CA
        elif switch=='secneig':
            self.__dict__['secneig']['CA']=''
            CA=joink(scd,'structures_CA_secneig',self.id+'.fake')
            if isthere(CA): self.__dict__['secneig']['CA']=CA
        elif switch=='megamix':
            self.__dict__['megamix']['CA']=''
            CA=joink(scd,'structures_CA_megamix',self.id+'.fake')
            if isthere(CA): self.__dict__['megamix']['CA']=CA
        elif switch=='arfavor1':
            self.__dict__['arfavor1']['CA']=''
            CA=joink(scd,'structures_CA_arfavor1',self.id+'.fake')
            if isthere(CA): self.__dict__['arfavor1']['CA']=CA
        return True


    def addSeqDat(self):

        """ add xxxx.SEQ secondary structure file"""
        from secondary.ss2seqdat import genSeqdat
        
        SEQ=joink(scratchdir,'structural_neighbors/seq',self.id+'.SEQ')
        buffer=genSeqdat(self.pdb300['seq'],self.secondary['mix'])
        open(SEQ,'w').write(buffer)
        self.pdb300['SEQ']=SEQ
        if 'seqdat' in self.pdb300.keys(): del self.pdb300['seqdat']
        return True
    

    def addhbnet(self,switch=''):

        """add hydrogen bond contact map"""

        from utilities.small_utilities import junkName

        junkf=junkName()
        hbnetx=joink(codedir,'f77/tertiary/hbnet.x')
        if switch not in self.__dict__.keys(): return False
        pt0=self.__dict__[switch]
        pt0['hbmap']={'list':[]}
        if pt0['CA']:
            cmd='%s -pdbf %s -secf %s -outf %s'%\
                 (hbnetx,pt0['CA'],self.pdb300['SEQ'],junkf)
            os.system(cmd)
            pt=pt0['hbmap']  ;  pt['nhelical']=0  ;  pt['nbeta']=0
            for line in open(junkf,'r').readlines()[1:]:
                a,b=line.split()[0:2] ; a=int(a) ; b=int(b)
                if abs(a-b)==3: #remember TASSER uses (i,i+3) to define a hydrogen bond
                    pt['nhelical']+=1
                    if b<self.L: b+=1
                    else: a-=1
                else: pt['nbeta']+=1
                pt['list'].append( (a,b) )
            pt['N']=len(pt['list'])
            os.system('/bin/rm '+junkf)
        return True


    def findHBPromiscous(self,switch=''):

        """find residues that participate in more than two hydrogen bonds"""
        if switch not in self.__dict__.keys(): return False
        n=[]
        for i in range(self.L+1): n.append([])
        for (a,b) in self.__dict__[switch]['hbmap']['list']:
            n[a].append(b)
            n[b].append(a)
        for i in range(self.L+1): n[i].sort()
        promiscuous=[]
        for i in range(self.L+1):
            if len(n[i])>2: promiscuous.append(i)
            elif len(n[i])==2 and abs(n[i][0]-n[i][1])==1: promiscuous.append(i)
        self.__dict__[switch]['hbmap']['promiscuous']=promiscuous
        return True
            

    def contactOrder(self,switch=''):
        
        """add contact order of the hydrogen bond map"""

        try:
            pt=self.__dict__[switch]['hbmap']
            hblist=pt['list']
            pt['co']={'all':0.0,'beta':0.0}
            for (a,b) in hblist:
                d=abs(a-b)
                pt['co']['all']+=d
                if d>4: pt['co']['beta']+=d #HB bond in helix has 4 as sequence separation
            pt['co']['all']=pt['co']['all']/(pt['N']*self.L)
            if pt['nbeta']:pt['co']['beta']=pt['co']['beta']/(pt['nbeta']*self.L)
        except:
            pass
        return True
            

    def addRg(self,switch=''):

        """add radius of gyration"""
        if switch not in self.__dict__.keys(): return False
        try:
            a,b=os.popen('Rg.x -a %s -b'%(self.__dict__[switch]['CA'])).readline().split()
            self.__dict__[switch]['Rg']=float(a)
            if 'Rg_expected' not in self.__dict__.keys(): self.Rg_expected=float(b)
        except:
            self.__dict__[switch]['Rg']=0.0        
        return True
        

    def parseDomains(self,switch=''):

        """find domains within structure"""
        from utilities.small_utilities import junkName
        try:
            pt=self.__dict__[switch]
            junkf=junkName()
            lines=open(pt['CA'],'r').readlines()
            for i in range(len(lines)):
                lines[i]=lines[i][0:21]+'X'+lines[i][22:]
            open(junkf,'w').write(''.join(lines))
            outline=os.popen('pdp '+junkf+' X').readline().split()[1]
            os.system('/bin/rm '+junkf)
            domains=outline.split('/')
            n=len(domains)
            pt['pdp']={'outline':outline,'n':n}
            Lis=[]
            for idomain in range(len(domains)):
                pt['pdp'][idomain]={'chunks':[],'nchunks':0}
                chunks=domains[idomain].split(',')
                L=0
                for chunk in chunks:
                    a,b=chunk.split('-')  ;  a=int(a)  ; b=int(b)  ;  L=L+b-a+1
                    pt['pdp'][idomain]['chunks'].append((a,b))
                    pt['pdp'][idomain]['nchunks']+=1
                pt['pdp'][idomain]['L']=L
                Lis.append(L)
            Lis.sort()  ;  Lis.reverse()  ; Lmax=Lis[0]  ; neff=0.0
            for Li in Lis: neff+=float(Li)/Lmax
            pt['pdp']['neff']=neff
        except:
            pass
        return True


    def TMscore2nat(self,id='',switch=''):
        """calculate TMscore to native"""
        from tmScore.tm_score import gettmCA
        pdb1=self.pdb300['CA']  ;  pdb2=''
        try:
            if switch=='secneig':
                if self.secneig['CA']:
                    pdb2=self.secneig['CA']
                    tm=gettmCA(pdb1,pdb2)  ;  self.secneig['tm2nat']=tm
                    print '%s %5.3f'%(self.id,tm)
        except:
            pass
        return False


    def sidechainCO(self,cutoff=2,switch=''):
        """calculate CO from sidechain contacts"""
        from tertiary.contactMap import genCM
        from tertiary.contactOrder import contactOrder
        if switch not in self.__dict__.keys(): return False
        pt=self.__dict__[switch]
        pt['sideco']=None
        if pt['CA']: pt['sideco']=contactOrder(genCM(pt['CA']),L=self.L,cutoff=cutoff)
        return True


    def prospector_neighbors(self,task='submit',joblist=[]):
        """from prospector results on pdb300 structure,
        find evolutionary related pdb300 structures"""

        if task=='submit': #submit prospector jobs
            from jobs.job import job,genJobList
            if not joblist: joblist=genJobList()
            outdir=self.simdir  ;  jobname='p.'+self.id
            args='-a '+self.id+' -b /tmp/jose/pdb_jul07 -c /tmp/jose/pdb_jul07'\
                  +' -e '+self.simdir+' -f '+joink(cad,'LIST.reduced')\
                  +' -j pdb'+self.id+'.out.prosp.tbz2'
            Job=job(name='prospector',exe='prospector.py',
                    exed=joink(codedir,'python/prospector/prospector_genomesmay06/modified'),
                    shared='/tmp/jose/pdb_jul07:'+joink(scd,'../pdb_jul07.tbz2'),
                    args=args)
            Job.qsub(jobname,outdir,wallt='0.249',mem_limit='990',file_size='250',ddisk='2500',
                     joblist=joblist)
            return False
        elif task=='filter':
            from utilities.small_utilities import unTARme
            from prospector.prospManager import templ_file
            self.pdb300['prospneig']={}  ;  pt=self.pdb300['prospneig']
            prospout=joink(self.simdir,'pdb'+self.id+'.out.prosp.tbz2')
            if isthere(prospout):
                neighbors=self.pdb300['nneighbor']['to pdb300']['neighbors']
                pt['outfile']=prospout  ;  pt['neighbors']=[]  ;  wd=unTARme(prospout)
                tmpls=templ_file( joink(wd,'pdbbpdborienrev',self.id+'rap3orienrev5s.pdb') )
                tmpls.loadTemplates()  ;  props_ids=tmpls.tpls.keys()
                for id in neighbors:
                    if id not in props_ids: pt['neighbors'].append(id)
                pt['n']=len(pt['neighbors'])
                pastry( '/bin/rm -rf '+wd )
            return True
        elif task=='filtered tmX2all':
            unfiltered={}
            for line in open(joink(ocd,self.id+'.list'),'r').readlines():
                (id,order,tm)=line.split()
                if order!='1': break
                unfiltered[id]=tm
            if self.pdb300['prospneig']:
                buf=`len(self.pdb300['prospneig']['neighbors'])`+'\n'
                for id in self.pdb300['prospneig']['neighbors']: buf+=id+'  '+unfiltered[id]+'\n'
                open(joink(cad,'tmX2all.filtered',self.id+'.dat'),'w').write(buf)
            return False
        
########################################################################

    
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    pdbid ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'simdir': directory to store simulation results
    |_'L'; sequence length
    |_'Rg_expected': expected radius of gyration for a single-domain sequence of lenth L
    |_'realpdbnneig': number of neigh of a real PDB structure to rest of real PDB str
    |_'secondary': directory holding info on secondary structure
    | |_'HE-': sec. assignments from running linker42sum
    | |_'HEC': sec. assignments from running pdb2ss.py
    | |_'mix': sec. assignments from mixing 'HE-' and 'HEC' methods
    | |_'content':dictionary with secondary content. Keys are H,E,C
    | |_'type': secondary class. Types are a(alpha),b(beta),c(alphabeta),d(other)
    |_'minimized': directory holding info on structures
    |  |_'outfile': tbz2 file containing structure PDB files
    |  |_'list': list of structure basenames
    |  |_basename:directory holding info on a particular structure
    |_'pdb300':
    | |_'seq': amino acid sequence of native structure
    | |_'CA': native file name
    | |_'Rg': radius of gyration
    | |_'sideco': contact order from sidechains
    | |_'pdp': info on protein domain parses (pdp) on str selected to represent the gen str
    | | |_'outline': output line from pdp
    | | |_'n': number of domains
    | | |_'neff': effective number of domains
    | | |_i:info on domain i. i runs from 0 to n-1. i is an integer
    | |   |_'chunks': list of domain chunks. Example ['1-45','80-156']
    | |   |-'nchunks': number of chunks in which the domain is divided
    | |   |_'L':domain size
    | |_'SEQ': xxxxx.SEQ file containing native secondary structure
    | |_'hbmap': hydrogen bond contact map info
    | | |_'list': list of bonds
    | | |_'N': number of bonds
    | | |_'nhelical': number of helical bonds
    | | |_'nbeta': number of beta bonds
    | | |_'promiscuous': list of residues that could participate in more than two bonds
    | | |_'co': contact order info for the hydrogen bond network
    | |   |_'all': consider all contacts when calculating the contact order
    | |   |_'beta': consider only contacts with sequence separation > 3
    | |_'nneighbor':
    |   |_'to secneig': comparison to generic library 'secneig'
    |   | |_'outfile': file containing generic structures similar to this model
    |   | |_'n': number of generic structures similar to this model
    |   |_'to pdb300':
    |   | |_'neighbors': list of structural neighbors
    |   | |_'n':
    |   |_'to pdb300.sd': compare to single-domain pdb300 set
    |     |_'n': number of neighbors to the single-domain pdb300 set
    |_'original': simulations with original generic TASSER
    | |_'nneighbor':
    |   |_'to original': comparison to itself
    |   | |_'outfile': file containing PDB structures that are similar to this model
    |   | |_'n': number of PDB structures that are similar to this model
    |   |_'to PDB': comparison to pdb300 library
    |     |_'outfile': file containing PDB structures that are similar to this model
    |     |_'n': number of PDB structures that are similar to this model
    |_'secneig': simulations with second neighbor interaction added to TASSER
    | |_'CA': structure selected to represent the generic structure
    | |_'Rg': radius of gyration of structure selected to represent the generic structure
    | |_'sideco': contact order from sidechains    
    | |_'tm2nat': TMscore to native
    | |_'pdp': info on protein domain parse (pdp) on str selected to represent the gen str
    | | |_'outline': output line from pdp
    | | |_'n': number of domains
    | |_'hbmap': hydrogen bond contact map info
    | | |_'list': list of bonds
    | | |_'N': number of bonds
    | | |_'nhelical': number of helical bonds
    | | |_'nbeta': number of beta bonds
    | | |_'promiscuous': list of residues that could participate in more than two bonds
    | | |_'co': contact order info for the hydrogen bond network
    | |   |_'all': consider all contacts when calculating the contact order
    | |   |_'beta': consider only contacts with sequence separation > 3
    | |_'tasser': info to hold tasser output
    | | |_'outfile'
    | |_'spicker': info to hold spicker output
    | | |_'outfile'
    | |_'nneighbor':
    |   |_'to PDB': comparison to PDB library
    |   | |_'outfile': file containing PDB structures that are similar to this model
    |   | |_'n': number of PDB structures that are similar to this model
    |   |_'to pdb_jul07': comparison to PDB library
    |   | |_'outfile': file containing PDB structures that are similar to this model
    |   | |_'n': number of PDB structures that are similar to this model
    |   |_'to secneig': comparison to itself
    |     |_'outfile': file containing PDB structures that are similar to this model
    |     |_'n': number of PDB structures that are similar to this model
    |_'megamix': secondary structure assignments from a mixture of dssp,linker42,Yang methods
    | |_'CA': structure selected to represent the generic structure
    | |_'simdir': directory
    | |_'Rg': radius of gyration of structure selected to represent the generic structure
    | |_'sideco': contact order from sidechains    
    | |_'seqdat': name of associated seq.dat file
    | |_'chain.dat': name of associated chain.dat file
    | |_'tasser':
    | | |_'outfile':
    | |_'spicker':
    | | |_'outfile':
    | |_'nneighbor':
    |   |_'to megamix': comparison to generic library 'megamix'
    |     |_'outfile': file containing generic structures similar to this model
    |     |_'n': number of generic structures similar to this model    
    |_'arfavor1': secondary structure assignments from a mixture of dssp,linker42,Yang methods
    | |_'simdir': directory
    | |_'tasser':
    | | |_'outfile':
    | |_'spicker':
    | | |_'outfile':
    | |_'nneighbor':
    |   |_'to arfavor1': comparison to generic library 'arfavor1'
    |     |_'outfile': file containing generic structures similar to this model
    |     |_'n': number of generic structures similar to this model    
    |_'snsd': secneig results but considering the structure is a single-domain, as given by pdp
      |_'nneighbor':
        |_'to snsd':
          |_'outfile': file containing generic single-dom structures that are similar to this model
          |_'n': number of generic single-dom structures that are similar to this model
    """
    print globals()['datatree'].__doc__


def doitall(id):

    """a recount of everything we had done on each varset

    Useful if we screw-up a particular dump file, because we can redo
    it all"""

    c=PEPSI(id=id)
    c.genDumpf(cl.repod)

    methods=('createSimdir','extractSS','genrmsinp','genChainDatFile','addL',
             'submitGenericTasser','initMininized','gatherStructures',
             'secondary_type')
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


def genPDBstructd():

    """create directory with CA atoms for every protein with associated PDB code"""

    os.system('/bin/mkdir -p '+cad) #directory to store structures
    for header in canlist:
        for pdbd in pdbds:
            pdbfin=joink(pdbd,header+'.pdb')
            pdbfout=joink(cad,header)
            if isthere(pdbfin):
                os.system('cat '+pdbfin+' |grep ATOM| grep -e " CA " > '+pdbfout)
                break
        if not isthere(pdbfout):
            sys.stderr.write('ERROR: no '+header+' structure found\n')
    #create file LIST.structures
    listf=joink(cad,'LIST.structures')
    open( listf,'w').write(`len(canlist)`+'\n')
    os.system('cat '+canlistf+' >> '+listf)
    #create list file
    os.system('/bin/cp '+canlistf+' '+joink(cad,'list'))

        
def extractSS(id=''):

        """access protein associated with each PDB id, then extract its
        secondary structure"""

        from utilities.small_utilities import junkName

        ss={} #store secondary structures
        if not isthere(secmaster): #gerenate secmaster file
            currd=os.getcwd() #current directory
            wd=joink(currd,junkName()) #temporary working directory
            os.system('/bin/mkdir -p '+wd)
            os.chdir(wd)
            os.system('/bin/cp /gpfs1/active/jose/code/f77/secondary/helix .')
            os.system(joink(pcd,'linker42sum')) #calculates secondary assignments
            os.system('/bin/mv secmaster '+cad) #move secondary structure assignments
            os.chdir(currd)
            os.system('/bin/rm -r '+wd) #remove temporary working directory
        pin=open(secmaster,'r')
        header=pin.readline()[1:6]
        while header:
            ss[header]=pin.readline().strip()
            header=pin.readline()[1:6]
        pin.close()

        if id:
            c=cl.loadCAN(id)
            if c.extractSS(ss=ss[c.id])!=False: c.pickleDump()
        else:
            remaining=len(canlist)
            for c in cl.iterate():
                print remaining,c.id                
                if c.extractSS(ss=ss[c.id])!=False: c.pickleDump()
                remaining-=1


def extracSSfromDSSP(id=''):
    """extract secondary structure from DSSP asignments"""
    indexl={}  ;  sspt=open( joink(scratchdir,'db/ss.txt'), 'r' )
    n=sspt.tell()  ;  line=sspt.readline()
    while line:
        if 'secstr' in line: indexl[ line[1:7] ]=n
        n=sspt.tell()  ;  line=sspt.readline()
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.extracSSfromDSSP(sspt,indexl)!=False: c.pickleDump()
        remaining-=1


def megamixSS(id=''):
    """We prepare a megamix assignment for secondary structure, whereby we use helical
    assignments from dssp and strand assignments from mix"""
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.megamixSS()!=False: c.pickleDump()
        remaining-=1


def genSeqDatFile(id='',switch=''):   
    """generate seq.dat file"""
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.genSeqDatFile(switch=switch)!=False: c.pickleDump()
        remaining-=1


def genrmsinp(id=''):

    """generate rmsinp file and sequence length"""
    
    if id:
        c=cl.loadCAN(id)
        if c.genrmsinp()!=False: c.pickleDump()
    else:
        cl.cycleOverList('genrmsinp',silent=True)


def genChainDatFile(id='',switch=''):
    """generate chain.dat file"""
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.genChainDatFile(switch=switch)!=False: c.pickleDump()
        remaining-=1


def submitGenericTasser(id='',protocol='initial',listf=''):

    """submit generic TASSER+SPICKER+AMBERMIN jobs"""

    from jobs.job import genJobList
    
    joblist=genJobList() #;Bye(joblist)
    if id:
        c=cl.loadCAN(id)
        if c.submitGenericTasser(joblist=joblist,protocol=protocol)!=False: c.pickleDump()
    else:
        ibatch=0
        remaining=len(canlist)
        if listf:
            for id in chomp(open(listf,'r').readlines()):
                c=cl.loadCAN(id)
                ibatch,record=c.submitGenericTasser(joblist=joblist,protocol=protocol,
                                                    ibatch=ibatch)
                if record!=False: c.pickleDump()
        else:
            for c in cl.iterate():
                print remaining,c.id                
                if c.submitGenericTasser(joblist=joblist,protocol=protocol)!=False: c.pickleDump()
                remaining-=1


def addL(id=''):

    """add sequence length"""

    if id:
        c=cl.loadCAN(id)
        if c.addL()!=False: c.pickleDump()
    else:
        cl.cycleOverList('addL',silent=False)
    

def checkChainDatFile(id=''):

    """check existence of chain.dat file. Try to produce one if non-existent"""

    if id:
        c=cl.loadCAN(id)
        if c.checkChainDatFile()!=False: c.pickleDump()
    else:
        cl.cycleOverList('checkChainDatFile',silent=False)
    

def initMininized(id=''):

    """init 'minimized' dictionary"""

    if id:
        c=cl.loadCAN(id)
        if c.initMininized()!=False: c.pickleDump()
    else:
        cl.cycleOverList('initMininized',silent=False)

    
def gatherStructures(id='',switch=''):

    """gather structures from the TASSER simulations"""

    #prepare structures directory
    if switch=='second neighbors':
        strd=joink(scd,'structures_CA_secneig') #structure repository directory
    elif switch=='secneig arfavor:1':
        strd=joink(scd,'structures_CA_arfavor1')
    elif switch=='megamix':
        strd=joink(scd,'structures_CA_megamix')
    else:
        strd=joink(scd,'structures_CA')
    if not isthere(strd): os.system('/bin/mkdir -p '+strd)
    #iterate over CAN objects
    if id:
        c=cl.loadCAN(id)
        if c.gatherStructures(switch=switch)!=False: c.pickleDump()
    else:
        remaining=len(canlist)
        for c in cl.iterate():
            print remaining,c.id 
            if c.gatherStructures(switch=switch)!=False: c.pickleDump()
            remaining-=1

    from utilities.small_utilities import junkName
    junk=junkName()
    pastry('ls -1 '+strd+'>'+junk+' && /bin/mv '+junk+' '+strd+'/list') #list of structures
    pastry('/bin/tar jcf '+strd+'.tbz2 '+strd)                  #tar directory


def secondary_type(id=''):

    """gather secondary structure content and estimate secondary class"""

    if id:
        c=cl.loadCAN(id)
        if c.secondary_type()!=False: c.pickleDump()
    else:
        cl.cycleOverList('secondary_type',silent=False)

def secondary_type_file():

    """create secondary.all.dat file"""

    buf=''
    for c in cl.iterate():
        if 'type' not in c.secondary.keys(): Bye(c.secondary)
        buf+=c.id+' '+c.secondary['type']+'\n'
    open( joink(scd,'secondary.all.dat'),'w').write(buf)

    
def submitTMalign_lists(exe_type='yang02',switch='',sizechunk=None,wallt=None):

    """Do TMalign of all against all"""

    from TMalign.TMalign_lists import TMalign_lists,TMalign_listsII
    if sizechunk: sizechunk=int(sizechunk)
    if wallt: wallt=float(wallt)
    if switch=='second neighbors':
        strd=joink(scd,'structures_CA_secneig') #directory of structures
        lib='/tmp/jose/structures_CA_secneig:'+strd+'.tbz2'
        listf=joink(strd,'list')
        outd=joink(scd,'tmAll2All_secneig') ; pastry('/bin/mkdir -p '+outd)
        TMalign_lists(lib,listf,outd,exe_type='shashi',sizechunk=sizechunk,wallt=wallt)
    if switch=='secneig arfavor:1':
        strd=joink(scd,'structures_CA_arfavor1') #directory of structures
        lib='/tmp/jose/structures_CA_arfavor1:'+strd+'.tbz2'
        listf=joink(strd,'list')
        outd=joink(scd,'tmAll2All_arfavor1') ; pastry('/bin/mkdir -p '+outd)
        TMalign_lists(lib,listf,outd,exe_type='shashi',sizechunk=sizechunk,wallt=wallt)
    if switch=='megamix':
        strd=joink(scd,'structures_CA_megamix') #directory of structures
        lib='/tmp/jose/structures_CA_megamix:'+strd+'.tbz2'
        listf=joink(strd,'list')
        outd=joink(scd,'tmAll2All_megamix') ; pastry('/bin/mkdir -p '+outd)
        TMalign_lists(lib,listf,outd,exe_type='shashi',sizechunk=sizechunk,wallt=wallt,prefix='me')
    elif switch=='second neighbors to PDB':
        strd=joink(scd,'structures_CA_secneig') #directory of structures
        lib='/tmp/jose/structures_CA_secneig:'+strd+'.tbz2'
        listf=joink(strd,'list')
        strd2=cad #directory of structures
        lib2='/tmp/jose/pdb:'+strd2+'.tbz2'
        listf2=joink(strd2,'list')
        outd=joink(scd,'tmAll2All_snPDB') ; pastry('/bin/mkdir -p '+outd)
        TMalign_listsII(lib,listf,lib2,listf2,outd,prefix='sn2p',exe_type='shashi',wallt=wallt)
    elif switch=='second neighbors to pdb_jul07':
        strd=joink(scd,'structures_CA_secneig') #directory of structures
        lib='/tmp/jose/structures_CA_secneig:'+strd+'.tbz2'
        listf=joink(strd,'list')
        strd2=pdb_jul07 #directory of structures
        lib2='/tmp/jose/pdb_jul07:'+strd2+'.tbz2'
        listf2=joink(strd2,'list')
        outd=joink(scd,'tmAll2All_snpdb_jul07') ; pastry('/bin/mkdir -p '+outd)
        TMalign_listsII(lib,listf,lib2,listf2,outd,prefix='sn2pj',exe_type='shashi',wallt=wallt)
    elif switch=='original to pdb_jul07':
        strd=joink(scd,'structures_CA_pdbgene2') #directory of structures
        lib='/tmp/jose/structures_CA_pdbgene2:'+strd+'.tbz2'
        listf=joink(strd,'list')
        strd2=pdb_jul07 #directory of structures
        lib2='/tmp/jose/pdb_jul07:'+strd2+'.tbz2'
        listf2=joink(strd2,'list')
        outd=joink(scd,'tmAll2All_gene2pdb_jul07') ; pastry('/bin/mkdir -p '+outd)
        TMalign_listsII(lib,listf,lib2,listf2,outd,prefix='o2pj',exe_type='shashi',wallt=wallt)
    else:
        strd=joink(scd,'structures_CA') #directory of structures
        lib='/tmp/jose/structures_CA'':'+strd+'.tbz2'
        listf=joink(strd,'list')
        outd=joink(scd,'tmAll2All')
        os.system('/bin/mkdir -p '+outd)
        TMalign_lists(lib,listf,outd,exe_type=exe_type,wallt=wallt) #qsub jobs
        

def submitTMalign_listsII():

    """Do TMalign of all against all pdb library"""

    from TMalign.TMalign_lists import TMalign_listsII
    
    strd=joink(scd,'structures_CA_pdb') #directory of structures
    lib='/tmp/jose/structures_CA_pdb:'+strd+'.tbz2'
    listf=joink(strd,'list')
    strd=joink(scratchdir,'structural_neighbors','pdb') #directory of PDB structures
    lib2='/tmp/jose/pdb:'+strd+'.tbz2'
    listf2=joink(strd,'list')
    outd=joink(scd,'tmAll2AllPDB')
    os.system('/bin/mkdir -p '+outd)
    TMalign_listsII(lib,listf,lib2,listf2,outd,prefix='tm2pdb') #qsub jobs
        

def submitTMalign_listsPDB(exe_type='yang02'):

    """Do TMalign of all real PDB structures against themselves"""

    from TMalign.TMalign_lists import TMalign_lists
    
    lib='/tmp/jose/pdb'':'+cad+'.tbz2'
    listf=joink(cad,'list')
    outd=joink(cad,'tmAll2All')
    os.system('/bin/mkdir -p '+outd)
    TMalign_lists(lib,listf,outd,exe_type=exe_type) #qsub jobs
        

def gathertmAll2Alldat():

    """create tmAll2All.dat files"""

    from TMalign.TMalign_lists import gathertmAll2Alldat
    
    tmAll2Alldir=joink(scd,'tmAll2All')
    tmAll2Allfile=joink(scd,'tmAll2All.dat')
    listofstruct=joink(scd,'structures_CA','list')
    gathertmAll2Alldat(tmAll2Alldir,tmAll2Allfile,listofstruct)
        

def gentmX2Allfiles(switch=''):

    """create tmX2All directory"""

    from TMalign.TMalign_lists import gentmX2Allfiles

    if switch=='secneig arfavor:1':
        tmAll2Allfile=joink(scd,'tmAll2All_arfavor1.dat')
        tmX2all=joink(scd,'tmX2all_arfavor1')
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strs=chomp(open(joink(scd,'structures_CA_arfavor1','list') ).readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
        return False
    elif switch=='megamix':
        tmAll2Allfile=joink(scd,'tmAll2All_megamix.dat')
        tmX2all=joink(scd,'tmX2all_megamix')
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strs=chomp(open(joink(scd,'structures_CA_megamix','list') ).readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
        return False
    elif switch=='second neighbors':
        tmAll2Allfile=joink(scd,'tmAll2All_secneig.dat')
        tmX2all=joink(scd,'tmX2all_secneig')
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strs=chomp(open(joink(scd,'structures_CA_secneig','list') ).readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
        return False
    elif switch=='second neighbors to PDB':
        tmAll2Allfile=joink(scd,'tmAll2All_snPDB.dat')
        tmX2all=joink(scd,'tmX2all_snPDB')
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strs=chomp(open(joink(scd,'structures_CA_secneig','list') ).readlines())
        strs+=chomp(open(joink(cad,'list')).readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
        return False
    elif switch=='second neighbors to pdb_jul07':
        tmAll2Allfile=joink(scd,'tmAll2All_snpdb_jul07.dat')
        tmX2all=joink(scd,'tmX2all_snpdb_jul07')
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strs=chomp(open(joink(scd,'structures_CA_secneig','list') ).readlines())        
        strs+=chomp(open(joink(pdb_jul07,'list')).readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
        return False
    elif switch=='original to pdb_jul07':
        tmAll2Allfile=joink(scd,'tmAll2All_gene2pdb_jul07.dat')
        tmX2all=joink(scd,'tmX2all_gene2pdb_jul07')
        if not isthere(tmX2all): os.system('/bin/mkdir -p '+tmX2all)
        strs=chomp(open(joink(scd,'structures_CA_pdbgene2','list') ).readlines())        
        strs+=chomp(open(joink(pdb_jul07,'list')).readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
        return False
    else:
        tmAll2Allfile=joink(scd,'tmAll2All.dat')
        tmX2all=joink(scd,'tmX2all')
        strs=chomp(open(joink(scd,'structures_CA','list') ).readlines())
        gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
        return False
    

def gentmX2AllfilesII():

    """create tmX2allpdb directory"""

    from TMalign.TMalign_lists import gentmX2Allfiles
    
    tmAll2Allfile=joink(scd,'tmAll2AllPDB.dat')
    tmX2all=joink(scd,'tmX2allpdb')
    strs=chomp(open(joink(scd,'structures_CA_pdb','list') ).readlines())
    strs+=chomp(open(joink(cad,'list') ).readlines())
    gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)


def gentmX2AllfilesPDB():

    """create tmX2All directory"""

    from TMalign.TMalign_lists import gentmX2Allfiles
    
    tmAll2Allfile=joink(cad,'tmAll2All.dat')
    tmX2all=joink(cad,'tmX2all')
    strs=chomp(open(joink(cad,'list') ).readlines())
    gentmX2Allfiles(tmAll2Allfile,tmX2all,strs)
        

def nneighbors():

    """gather some statistics on number of neighbors, dump into nneighbors.dat"""

    from secondary.seqdat2ss import translateSec
    #gather who is helical, who is beta, and so on
    buff='#id    L nneighbors\n'
    tmX2alld=joink(scd,'tmX2all')
    tmX2allL=joink(tmX2alld,'list')
    list=[ l[0:5] for l in open(tmX2allL,'r').readlines()]
    for id in list:
        #estimate secondary structure content
        secdat=translateSec(joink(scd,id,'seq.dat'))
        n=len(secdat)
        fa=float( len(re.compile('H').findall(secdat)) )/n
        fb=float( len(re.compile('E').findall(secdat)) )/n
        fc=1-fa-fb
        if fa>0.4 and fb<0.1:   type='a'
        elif fb>0.4 and fa<0.1: type='b'
        elif fa>0.2 and fb>0.2: type='c'
        else: type='d'
        L=cl.loadCAN(id).L
        n=int( os.popen('head -1 '+joink(tmX2alld,id+'.1.pdb.dat')).readline() )
        buff+=id+' %3d %4d'%(L,n)+' '+type+'\n'
        print id+' %3d %4d'%(L,n)+' '+type
    open( joink(tmX2alld,'nneighbors.dat'),'w').write(buff)


def view(id):

    """viewver of object properties"""

    cl.loadCAN(id).view()
    

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


def nneighbor(id='',switch=''):

    """enter nneighbor info"""

    if id:
        c=cl.loadCAN(id)
        if c.nneighbor(switch=switch)!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            print remaining,c.id
            if c.nneighbor(switch=switch)!=False: c.pickleDump()
            remaining-=1
       
    
def output_coverage(switch=''):

    """output coverage of structures of one library to structures of other library"""

    version='' ;  flavor=''
    if switch=='original to original':
        version='original' ;  flavor='to original'
    elif switch=='original to PDB':
        version='original' ;  flavor='to PDB'
    elif switch=='secneig to secneig':
        version='secneig'  ;  flavor='to secneig'
    elif switch=='arfavor1 to arfavor1':
        version='arfavor1'  ;  flavor='to arfavor1'
    elif switch=='megamix to megamix':
        version='megamix'  ;  flavor='to megamix'
    elif switch=='secneig to pdb_jul07':
        version='secneig'  ;  flavor='to pdb_jul07'
    elif switch=='second neighbors to PDB':
        version='secneig'  ;  flavor='to PDB'
    elif switch=='pdb300 to secneig':
        version='pdb300'  ;  flavor='to secneig'
    elif switch=='pdb300 to pdb300':
        version='pdb300'  ;  flavor='to pdb300'
    elif switch=='snsd to snsd':
        version='snsd'  ;  flavor='to snsd'

    Ls=['0-25','25-50','50-75','75-100','100-125','125-150','150-175','175-200',
        '200-225','225-250','250-275','275-300']
    hist={'N':{'all':[0]*12,'abc':[0]*12,'a':[0]*12,'b':[0]*12,'c':[0]*12,'d':[0]*12},
          'cov':{'all':[0.]*12,'abc':[0.]*12,'a':[0.]*12,'b':[0.]*12,'c':[0.]*12,'d':[0.]*12} }

    for c in cl.iterate():
        try:
            if switch=='real PDB to real PDB':
                n=c.realpdbnneig
            else:
                n=c.__dict__[version]['nneighbor'][flavor]['n']
            i=c.L/25
            hist['N']['all'][i]+=1
            sec=c.secondary['type']
            hist['N'][sec][i]+=1
            if sec!='d': hist['N']['abc'][i]+=1
            if n>0: #covered by the PDB
                hist['cov']['all'][i]+=1
                hist['cov'][sec][i]+=1
                if sec!='d': hist['cov']['abc'][i]+=1
        except: pass
    buf= '#  L      N  cov Nabc covabc Na  cov.a  Nb cov.b  Nc cov.c  Nd cov.d\n'
    buf+='#  1      2   3   4     5    6     7    8    9    10  11    12  13\n'
    for i in range(0,12):
        buf+='%7s'%(Ls[i])
        for sectype in ('all','abc','a','b','c','d'):
            if hist['N'][sectype][i]:
                hist['cov'][sectype][i]=hist['cov'][sectype][i]/hist['N'][sectype][i]
            buf+=' %4d %4.2lf'%(hist['N'][sectype][i],hist['cov'][sectype][i])
        buf+='\n'
    print buf


def addPDBnneigh(id=''):
    
    """ add number of neighbors for real PDB structure to other real PDB structures"""

    if id:
        c=cl.loadCAN(id)
        if c.addPDBnneigh()!=False: c.pickleDump()
    else:
        remaining=len(cl.canlist)
        for c in cl.iterate():
            print remaining,c.id
            if c.addPDBnneigh()!=False: c.pickleDump()
            remaining-=1
            

def outputRealnneig(id=''):

    """output number of neighbors for the real PDB's to themselves!"""

    idlist=cl.canlist
    if id: idlist=[id,]
    for id in idlist:
        c=cl.loadCAN(id)
        print '%s %4d %s'%(c.id,c.realpdbnneig,c.secondary['type'])
        

def addSeq(id=''):

    """add sequence of native PDB protein"""

    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.addSeq()!=False: c.pickleDump()
        remaining-=1
            

def addhbnet(id='',switch=''):

    """add hydrogen bond contact map"""
    
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.addhbnet(switch=switch)!=False: c.pickleDump()
        remaining-=1


def findHBPromiscous(id='',switch=''):

    """find residues that participate in more than two hydrogen bonds"""
    
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.findHBPromiscous(switch=switch)!=False: c.pickleDump()
        remaining-=1


def contactOrder(id='',switch=''):
        
    """add contact order"""

    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.contactOrder(switch=switch)!=False: c.pickleDump()
        remaining-=1


def addNativeFile(id='',switch=''):

    """add native file path"""
    
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.addNativeFile(switch=switch)!=False: c.pickleDump()
        remaining-=1


def addSeqDat(id=''):

    """ add xxxx.SEQ secondary structure file"""

    os.system('/bin/mkdir -p '+joink(scratchdir,'structural_neighbors/seq'))
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(cl.canlist)
    for id in idlist:
        c=cl.loadCAN(id)
        print remaining,c.id
        if c.addSeqDat()!=False: c.pickleDump()
        remaining-=1


def hbnetFeatures(id='',switch='secneig'):

    """output feature of hydrogen bond map"""

    print '#T: secondary structure type'
    print '#L: sequence length'
    print '#N: ratio number of hydrogen bonds of generic to native structure'
    print '#co_all: ratio contact order of all HB contacts of generic to native structure'
    print '#co_beta: ratio contact order of HB beta contacts of generic to native structure'
    print '#Rge Rgn Rgs: expected, native, and generic radius of gyration'
    print '#nfn nfg: effective number of domains for native and generic structures'
    print '#pp: number of promiscuous residues in the native structure'
    print '#sp: number of promiscuous residues in the generic structure'
    print '#npp: number of neighbors of native structure to native PDB300'
    print '#nps: number of neighbors of native structure to generic PDB300'
    print '#nss: number of neighbors of generic structure to generic PDB300'
    print '#nsp: number of neighbors of generic structure to native PDB300'
    print '#id   T  L     N   co_all  Nbeta co_beta Rge  Rgn  Rgs  nfn  nfg  pp  sp  npp  nps  nss  nsp'
    print '# 1   2  3     4     5       6     7      8    9    10   11   12  13  14   15  16    17   18'
    idlist=cl.canlist
    if id: idlist=[id,]
    for id in idlist:
        c=cl.loadCAN(id)
        try:
            ptp=c.pdb300['hbmap']
            pts=c.__dict__['secneig']['hbmap']
            buf='%s %s %3d %6.3lf %6.3lf '%(c.id,c.secondary['type'],c.L,float(pts['N'])/ptp['N'],
                                         pts['co']['all']/ptp['co']['all'])
            if ptp['nbeta']:
                buf+='%6.3lf %6.3lf '%(float(pts['nbeta'])/ptp['nbeta'],
                                       pts['co']['beta']/ptp['co']['beta'])
            else: buf+='  N/A     N/A '
            buf+='%4.1lf %4.1lf %4.1lf '%(c.Rg_expected,c.pdb300['Rg'],c.secneig['Rg'])
            buf+='%4.2lf %4.2lf '%(c.pdb300['pdp']['neff'],c.secneig['pdp']['neff'])
            buf+='%3d %3d '%(len(ptp['promiscuous']),len(pts['promiscuous']))
            buf+='%4d %4d %4d %4d'%(c.pdb300['nneighbor']['to pdb300']['n'],
                            c.pdb300['nneighbor']['to secneig']['n'],
                            c.secneig['nneighbor']['to secneig']['n'],
                            c.secneig['nneighbor']['to PDB']['n']
                            )
            print buf
        except:
            pass
            

def addRg(id='',switch=''):
            
    """add radius of gyration"""
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    for id in idlist:
        print remaining,id
        c=cl.loadCAN(id)
        if c.addRg(switch=switch)!=False: c.pickleDump()
        remaining-=1

    
def parseDomains(id='',switch=''):

    """find domains within structure"""

    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    for id in idlist:
        print remaining,id
        c=cl.loadCAN(id)
        if c.parseDomains(switch=switch)!=False: c.pickleDump()
        remaining-=1


def filterBy_neff(switch='pdb300',neff=1.01):

    """return a list of targets that have efective number of domains below neff argument"""

    neff=float(neff)
    filtdom=[]
    idlist=cl.canlist
    for id in idlist:
        c=cl.loadCAN(id)
        try:
            if c.__dict__[switch]['pdp']['neff']<neff: filtdom.append(id)
        except:
            pass
    print '\n'.join(filtdom)


def neighFilt_neff(id='',switch='pdb300',switchfilter='pdb300'):

    """print number of neighbors for each query if we remove multiple domains"""

    filtdom=[]
    if switchfilter=='pdb300':
        filtdom=chomp(open(joink(scd,'analysis/pdb300.singledomain.list'),'r').readlines())
    elif switchfilter=='secneig':
        filtdom=chomp(open(joink(scd,'analysis/secneig.singledomain.list'),'r').readlines())

    idlist=cl.canlist
    if id: idlist=[id,]
    for id in idlist:
        c=cl.loadCAN(id)  ;  neigfile=''  ;  pt=None
        try:
            pt=c.__dict__[switch]['nneighbor']
        except:
            continue
        pt[switchfilter+'.sd']={'n':0}  ;  extract=None
        if switchfilter=='pdb300':
            if switch=='pdb300':
                neigfile=joink(scd,'../pdb/tmX2all/'+c.id+'.dat')
                extract=lambda line: line[0:5]
            elif switch=='secneig':
                neigfile=joink(scd,'tmX2all_secneig/'+c.id+'.fake.dat')
                extract=lambda line: line[0:5]
        elif switchfilter=='secneig':
            if switch=='secneig':
                extract=lambda line: line[0:5]
                neigfile=joink(scd,'tmX2all_secneig/'+c.id+'.fake.dat')
        try:
            pin=open(neigfile,'r')  ;  n=int(pin.readline())
            neighborsfiltered=[]
            if n>0:
                line=pin.readline()
                while line:
                    neighbor=extract(line) #; print neighbor
                    if neighbor in filtdom: neighborsfiltered.append(neighbor)
                    line=pin.readline()
                pt[switchfilter+'.sd']['n']=len(neighborsfiltered)
            print '%s %4d %4d'%(id,len(neighborsfiltered),n)
        except:
            pass
        

def create_tmX2all_snsd_snsd():
    """create the tmX2all_snsd_snsd directory

    The directory contains neighbor files for each secneig structure whose secondar
    structure originates in a single-domain PDB structure. All neighbors are also secneig
    structures whose secondary structures originates in single-domain PDB structures."""          

    xxx='/gpfs1/scratch/jose/structural_neighbors/out/tmX2all_secneig/'
        
    ids=chomp( open(scd+'/analysis/pdb300.singledomain.list','r').readlines() )
    remaining=len(ids)
    outdx=joink(scd,'tmX2all_snsd_snsd')  ;  os.system('/bin/mkdir -p '+outdx)
    list=[]
    for id in ids:
        print id,remaining
        c=cl.loadCAN(id)
        try:
            neigf=c.secneig['nneighbor']['to secneig']['outfile']
            neigfb=os.path.basename(neigf)
            ptin=open(neigf,'r')  ;  n=int(ptin.readline())  ;  buffer=[]
            for i in range(0,n):
                line=ptin.readline()
                if line[0:5] in ids: buffer.append(line)
            buffer=`len(buffer)`+'\n'+''.join(buffer)
            open(joink(outdx,neigfb),'w').write(buffer)
            list.append(neigfb)
        except:
            pass
        remaining-=1
    list=`len(list)`+'\n'+'\n'.join(list)
    open(joink(outdx,'LIST.files'),'w').write(list)
        

def TMscore2nat(id='',switch=''):
    """calculate TMscore to native"""
    
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    for id in idlist:
        #print remaining,id
        c=cl.loadCAN(id)
        if switch=='secneig':
            if c.TMscore2nat(switch=switch)!=False: c.pickleDump()
        remaining-=1
        

def sidechainCO(id='',cutoff=2,switch='pdb300'):
    """calculate CO from sidechain contacts"""
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    for id in idlist:
        print remaining,id
        c=cl.loadCAN(id)
        if c.sidechainCO(cutoff=cutoff,switch=switch)!=False: c.pickleDump()
        remaining-=1


def compareSec(id=''):
    """compare secondary structure assignments"""
    from secondary.seqdat2ss import translateSec
    print '#          HE-        mix         seqdat       dssp      megamix'
    print '# id    H   E  -    H   E   C   H   E   C   H   E   C   H   E   C sectype'
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    for id in idlist:
        c=cl.loadCAN(id)  ;  pt=c.secondary  ;  buf=c.id+' '
        for x in ('H','E','-'): buf+='%3d '%(int(len(re.compile(x).findall(pt['HE-']))))
        for x in ('H','E','C'): buf+='%3d '%(int(len(re.compile(x).findall(pt['mix']))))
        sec=translateSec(joink(c.simdir,'seq.dat'))
        for x in ('H','E','C'): buf+='%3d '%(int(len(re.compile(x).findall(sec))))
        for x in ('H','E','C'): buf+='%3d '%(int(len(re.compile(x).findall(pt['dssp']))))
        for x in ('H','E','C'): buf+='%3d '%(int(len(re.compile(x).findall(pt['megamix']))))
        print buf+pt['type']
        

def checkSomething(id='',listfile=''):

    """do something!"""
    idlist=cl.canlist
    if listfile: idlist=chomp(open(listfile,'r').readlines())
    if id: idlist=[id,]
    for id in idlist:
        c=cl.loadCAN(id)
        if c.L > 250:
            if c.pdb300['nneighbor']['to pdb300']['n']>0:
                pt=open(joink(ocd,c.id+'.list'),'r')
                line=pt.readline()
                id2=''
                while line:
                    id2,n,tm=line.split()
                    if n=='1' and id2!=c.id: break
                    if n!='1': break
                    line=pt.readline()
                if id2:
                    c2=cl.loadCAN(id2)
                    cmd='/gpfs1/u/jose/bin/align %s %s 3'%(c.pdb300['seq'],c2.pdb300['seq'])
                    lines=os.popen(cmd).readlines()
                    identical=float(lines[4].split()[2])
                    print '%s %3d %s %3d %s %4.2lf'%\
                          (c.id,c.L,id2,c2.L,tm,identical/min(c.L,c2.L))


def contactOrderInfo(id='',switch='',filterRg=1.0,filterNHB=1.0):
    filterRg=float(filterRg)  ;  filterNHB=float(filterNHB)
    idlist=cl.canlist
    if id: idlist=[id,]
    remaining=len(idlist)
    if not switch: switch='secneig'
    print '#id    L   nat nat   N_HB nat   nat  fake  fake   N_HB fake  fake  sec-type'
    print '#          Rg  side       HB    beta  Rg   side         HB   beta'
    for id in idlist:
        c=cl.loadCAN(id)
        if switch not in c.__dict__.keys(): continue
        pt=c.__dict__[switch]
        try:
            if c.pdb300['Rg']>(1+filterRg)*c.Rg_expected or \
               c.pdb300['Rg']<(1-filterRg)*c.Rg_expected: continue #filter by Rg
            if c.pdb300['hbmap']['N']>(1+filterRg)*pt['hbmap']['N'] or \
               c.pdb300['hbmap']['N']<(1-filterRg)*pt['hbmap']['N']: continue #filter by Rg
            buf='%s %3d %4.1lf %5.3lf %3d %5.3lf %5.3lf'%\
                 (id,c.L,c.pdb300['Rg'],c.pdb300['sideco'],c.pdb300['hbmap']['N'],
                  c.pdb300['hbmap']['co']['all'],c.pdb300['hbmap']['co']['beta'])
            if 'Rg' in pt.keys():
                buf+=' %4.1lf '%(pt['Rg'])
            if 'sideco' in pt.keys():
                buf+=' %5.3lf '%(pt['sideco'])
            if 'hbmap' in pt.keys():
                buf+=' %3d %5.3lf '%(pt['hbmap']['N'],pt['hbmap']['co']['all'])
            if 'beta' in pt['hbmap']['co'].keys():
                buf+=' %5.3lf'%(pt['hbmap']['co']['beta'])
            buf+=' %s'%(c.secondary['type'])
            print buf
        except:
            pass


def prospector_neighbors(*kargs,**kwargs):
    """from prospector results on pdb300 structure, find evolutionary related pdb300 structures
    example: ... -c id=101m_"""
    if ('task','submit') in kwargs.items(): kwargs['pause']=0.6 #pause 0.6sec between qsub's
    elif ('task','filtered tmX2all') in kwargs.items():
        os.system('/bin/mkdir -p '+joink(cad,'tmX2all.filtered'))
    cl.cycleOverList('prospector_neighbors',*kargs,**kwargs)
            
    
def help(*kargs):

    """list available services"""
    
    services=['createList','datatree','doitall','updateCollection','initCollection',
              'createSimdir','genPDBstructd','extractSS','genrmsinp','genChainDatFile',
              'addL','submitGenericTasser','checkChainDatFile','initMininized',
              'gatherStructures','gathertmAll2Alldat','nneighbors','gentmX2AllfilesII',
              'secondary_type','secondary_type_file','gentmX2AllfilesPDB',
              'createTasserSpickerStructuresDicts','addPDBnneigh','addSeq',
              'addNativeFile','addSeqDat','findHBPromiscous','addRg','parseDomains',
              'filterBy_neff','TMscore2nat','create_tmX2all_snsd_snsd',
              'sidechainCO']
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
