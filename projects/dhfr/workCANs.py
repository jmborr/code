#!/usr/bin/python

import os,sys,re,pdb
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir,projectsdir
from utilities.small_utilities import chomp,Bye,junkName
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
joink=os.path.join; jk=joink
isthere=os.path.exists
basename=os.path.basename

my_project='dhfr'                   #substitute 'my_project'
pcd=joink(codedir,'projects',my_project)  #code directory
pjd=joink(projectsdir,my_project)         #projects directory
pjdout=joink(projectsdir,my_project,'out')
scd=joink(scratchdir,my_project)          #scratach directory
scd=joink(scratchdir,my_project,'out')
collf=joink(pjd,'varsets.can')            #you may substitute scd for pjd
canlistf=joink(pjd,'out','list')          #you may substitute scd for pjd
repod=joink(pjd,'repod')                  #you may substitute scd for pjd

def distrE(E,type='BASIS'):
    """gaussian distributions for energy resolution
    energy units are 1E-6eV"""
    from math import exp,pi
    #fwhm={'original':32,'BASIS':3}
    z=1.0/2.354820045
    sigma={'original':z*32,'BASIS':z*3}
    ss=sigma[type] ; ss=ss*ss
    return exp(-E*E/(2*ss))/sqrt(2*pi*ss)

def convT(t,type='BASIS'):
    """fourier transform of distrE
    Time units are in picoseconds"""
    from math import exp,pi
    sigma={'original':0.049,'BASIS':0.521}
    ss=sigma[type] ; ss=ss*ss
    return exp(-t*t/(2*ss))/sqrt(2*pi*ss)
    
    
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


    def foo(self):
        print self.id
        return False #True will pickledump


    def gensimdir(self):
        """simulation directory"""
        self.simdir=joink(pjdout,self.id)
        if not isthere(self.simdir):
            os.system('/bin/mkdir -p '+self.simdir)
        return True
    
    def getinitPDB(self,recover=''):
        """retrieve PDB file"""
        from Bio.PDB import PDBList
        if recover:
            self.pdbfile=joink(self.simdir,'pdb'+self.id[0:4]+'.ent')
            return True
        self.pdbfile=PDBList().retrieve_pdb_file(self.id[0:4].upper(),pdir=self.simdir)
        return True

    def storeL(self):
        """store protein length"""
        self.L=53
        return True
    
    def parsePDB(self,recover=''):
        """parse PDB file for tleap use

        (1)remove non ATOM lines
        (2)remove other chains than ' ' or 'A'
        (3)remove hydrogen and 'D' atoms"""
        if recover:
            self.pdbparsed=joink(self.simdir,'parsed_'+basename(self.pdbfile))
            return True
        buf=''
        for line in open(self.pdbfile,'r'):
            p=line            
            if p[16] not in (' ','A'): continue
            if p[0:4]=='LINK': continue
            if p[0:4]!='ATOM' and p[17:20]!=' FE': continue
            if p[13] in ('H','D'): continue
            if p[14:16]=='DZ': continue
            if self.id=='1vcx_Fe' or self.id=='1vcx_Fe_D':
                if 'CYS' in p: p=p.replace('CYS','CYM') #CYS-->CYM
                if 'FE    FE' in p: p=p.replace('FE    FE','FE   FE3')
            buf+=p
        buf+='TER\nEND\n'
        if self.id=='1vcx_Fe' or self.id=='1vcx_Fe_D':
            buf=buf.split('\n')
            buf='\n'.join( buf[0:-4]+['TER',]+buf[-4:] )
        self.pdbparsed=joink(self.simdir,'parsed_'+basename(self.pdbfile))
        open(self.pdbparsed,'w').write(buf)
        return True

    def getcharge(self):
        """find net charge of parsed pdb file"""
        buf='source leaprc.ff98\n'
        if self.isDeutd: buf='source leaprc.ff98Dsimple\n'
        if self.id in ('1vcx_Fe','1vcx_Fe_D'): buf+='loadOff FE3.lib\n'
        buf+=self.id+'=loadPdb '+self.pdbparsed+'\n'
        buf+='charge '+self.id+'\n'
        buf+='quit\n'
        script=junkName() ; open(script,'w').write(buf)
        tlout=junkName()
        os.system('tleap -f '+script+' > '+tlout)
        pattern=re.compile('unperturbed charge:\s*(-*\d+\.\d+)')
        match=pattern.search(''.join(open(tlout,'r').readlines()))
        if match:
            self.Q=float(match.group(1)) ;print self.Q
        else:
            sys.stderr.write('ERROR getcharge failed\n')
        os.system('/bin/rm '+script+' '+tlout)
        return True

    def get_bounding_box(self):
        """extract the bounding box from the solvated system"""
        prmtop=open(self.parmfile,'r') ; line=prmtop.readline()
        self.box=None
        while line:
            if 'FLAG BOX_DIMENSIONS' in line:
                line=prmtop.readline()
                self.box=' '.join(prmtop.readline().split()[-3:])
                break
            line=prmtop.readline()
        print self.box
        return True

    def gen_tleap_in(self,mode='solvated',recover=''):
        """generate input tleap
        (1)add missing hydrogens
        (2)solvate box with 10 Angstroms at least around the protein
        (3)add counterions"""
        self.intleap=''
        buf='source leaprc.ff98\n'
        if self.isDeutd: buf='source leaprc.ff98Dsimple\n'
        buf+=self.id+'=loadPdb '+self.pdbparsed+'\n'
        if mode=='solvated':
            buf+='solvateBox '+self.id+' SPCBOX 10\n'
            if int(self.Q)!=0:
                ion='Na+'
                if self.Q>0: ion='Cl-'
                buf+='addIons2 '+self.id+' '+ion+' 0\n'
            self.parmfile=joink(self.simdir,self.id+'.top')
            self.crdfile=joink(self.simdir,self.id+'.crd')
            buf+='saveAmberParm '+self.id+' '+self.parmfile+' '+\
                  self.crdfile+'\n'
            buf+='quit\n'
            self.intleap=joink(self.simdir,'in.tleap')
            if recover:return True
            open(self.intleap,'w').write(buf)
        elif mode=='unsolvated':            
            if 'unsolvated' not in self.__dict__: self.unsolv={}
            self.unsolv['parmfile']=joink(self.simdir,self.id+'.uns.top')
            self.unsolv['crdfile']=joink(self.simdir,self.id+'.uns.crd')
            buf+='set '+self.id+' box {'+self.box+'}\n'
            buf+='saveAmberParm '+self.id+' '+self.unsolv['parmfile']+' '+\
                  self.unsolv['crdfile']+'\n'
            buf+='quit\n'
            self.unsolv['intleap']=joink(self.simdir,'in.uns.tleap')
            if recover:return True
            open(self.unsolv['intleap'],'w').write(buf)
        return True

    def genparmtop(self,mode='solvated',recover=''):
        """generate paramenter, topology, and preminPDB files"""
        if mode=='solvated':
            if not self.intleap: sys.stderr.write('no in.tleap file')
            self.preminPDB=joink(self.simdir,self.id+'.premin.pdb')
            self.outtleap=joink(self.simdir,'out.tleap')
            if recover: return True
            os.system('tleap -f '+self.intleap+' &> '+self.outtleap)
            #create pre-minimization PDB file
            os.system('ambpdb -p '+self.parmfile+' < '+self.crdfile+\
                      ' > '+self.preminPDB)
        elif mode=='unsolvated':
            if 'unsolv' not in self.__dict__.keys(): self.unsolv={}
            self.unsolv['parmfile']=joink(self.simdir,self.id+'.uns.top')
            self.unsolv['crdfile']=joink(self.simdir,self.id+'.uns.crd')
            self.unsolv['preminPDB']=joink(self.simdir,self.id+\
                                           '.uns.premin.pdb')
            if recover: return True
            os.system('ambpdb -p '+self.unsolv['parmfile']+\
                      ' < '+self.unsolv['crdfile']+\
                      ' > '+self.unsolv['preminPDB'])
        return True

    def deuteratePDB(self):
        """Insert hydrogens and change to deuterated amino acid names"""
        if self.id!='1vcx_Fe_D':
            sys.stderr.write('not 1vcx_Fe_D\n')
            return False
        from amber.amber9 import leap,tleap
        from seq.letters import three2one
        tmpdb=junkName() #; print tmpdb #will include hydrogens
        buf="""source leaprc.ff98
loadOff FE3.lib
loadamberparams /home/software/amber9/dat/contrib/FE.frcmod
x=loadPdb """+self.pdbparsed+"""
savepdb x """+tmpdb
        tleap(seq=buf).go()
        lp=leap()
        buf=''
        for l in open(tmpdb,'r').readlines():
            if l[0:5]=='ATOM ':
                AAA=l[17:20]
                if AAA in lp.isstandarAAAname.keys():
                    DDD='D'+three2one[ lp.basenameAAA[AAA] ]
                    if lp.isstandarAAAname[AAA]: DDD+='O'
                    else: DDD+=AAA[-1]
                    l=l[:17]+DDD+l[20:]
            buf+=l
        open(tmpdb,'w').write(buf)
        os.system('/bin/mv '+tmpdb+' '+self.pdbparsed)
        
    def genPackageForMin(self):
        """package input for minimization"""
        name=self.id+'.ForMin.tbz2'
        parentdir=joink(os.getcwd(),junkName())
        os.system('/bin/mkdir -p '+parentdir+'/'+self.id)        
        os.chdir( joink(parentdir,self.id) )
        os.system('/bin/cp -r /home/software/amber9/Equil .')
        buf=''.join(open('./Equil/submit','r').readlines())
        buf=buf.replace('somejobname',self.id)
        buf=buf.replace('projectdir','rubredoxin')
        buf=buf.replace('subdir',self.id)
        buf=buf.replace('pid',self.id)
        buf=buf.replace('theIon','Na')
        open('./Equil/submit','w').write(buf)
        os.system('/bin/cp '+self.preminPDB+' Equil/srt/'+self.id+'.pdb')
        os.system('/bin/cp '+self.parmfile+' Equil/srt/'+self.id+'.top')
        os.system('/bin/cp '+self.crdfile+' Equil/srt/'+self.id+'.crd.srt')
        os.chdir(parentdir)
        cmd='tar jcf '+name+' '+self.id+' && mv '+name+' '+self.simdir
        os.system(cmd)
        self.packForMin=joink(self.simdir,name)
        return False

    def genPDBid(self):
        self.pdbid=self.id[0:4]
        return True

    def createNMRdist(self,recover=''):
        """create NRM distance file according to Pratul"""
        templ=""" &rst
 iat=_AT1_,_AT2_,
    iresid=0,
         r1=0.0, r2=2.2, r3=2.3, r4=999.0,
         rk2=25.0, rk3=25.0,
 &end
"""
        if self.id=='1vcx_Fe' or self.id=='1vcx_Fe_D':
            self.NMR={} ; p=self.NMR
            p['amberdistf']=joink(self.simdir,self.id+'.amber.dist')
            if recover: return True
            at1s=('82','133','560','603') ; at2='789'
            buf=''
            for at1 in at1s:
                buff=templ.replace('_AT1_',at1)
                buf+=buff.replace('_AT2_',at2)
            open(p['amberdistf'],'w').write(buf)
        return True
    
    def createNMRdist2(self):
        """create NRM distance file according to AMBER tutorial"""
        if self.id=='1vcx_Fe':
            self.NMR={} ; p=self.NMR
            #create simple distance input
            p['distf']=joink(self.simdir,self.id+'.dist')
            open(p['distf'],'w').write(\
            """  5 CYM  SG   54   FE   FE  2.300
  8 CYM  SG   54   FE   FE  2.300
  38 CYM  SG   54   FE   FE  2.300
  41 CYM  SG   54   FE   FE  2.300""")
            #create custom mapping file
            p['mapf']=joink(self.simdir,'map.DG-AMBER')
            os.system('/bin/cp $AMBERHOME/dat/map.DG-AMBER '+p['mapf'])
            open(p['mapf'],'a').write(\
                """
RESIDUE FE
MAPPING         FE      =       FE""")
            #generate amber distance file
            p['amberdistf']=joink(self.simdir,self.id+'.amber.dist')
            os.system('makeDIST_RST -upb '+p['distf']+\
                      ' -pdb '+self.pdbparsed+' -rst '+p['amberdistf']+\
                      ' -map '+p['mapf'])
            buf=''.join(open(p['amberdistf'],'r').readlines())
            subs={'r1= 1.30':'r1= 1.70', 'r2= 1.80':'r2= 2.20',
                  'rk2=20.0':'rk2=25.0', 'rk3=20.0':'rk3=25.0'}
            for (a,b) in subs.items(): buf=buf.replace(a,b)
            open(p['amberdistf'],'w').write(buf)
            return True

    def setEquilInput(self):
        """import certain file names neccessary for equilibrium"""
        self.preminPDB=joink(self.simdir,self.id+'.premin.pdb')
        self.crdfile=joink(self.simdir,self.id+'.crd')
        self.parmfile=joink(self.simdir,self.id+'.top')
        return True
        
    def relax(self,T='300'):
        """minimize and relax the premin pdb structure"""
        #copy template Equil directory and a few input files
        equil0=os.environ['AMBERHOME']+'/minprotoc/Equil'
        
        equil0=os.environ['AMBERHOME']+'/minprotoc/Equil'
        Td=joink(self.simdir,'T_'+T) ; os.system('/bin/mkdir -p '+Td)
        equil=joink(Td,basename(equil0))
        pdbf=joink(equil,'srt',self.id+'.pdb')
        crdf=joink(equil,'srt',self.id+'.crd.srt')
        parf=joink(equil,'srt',self.id+'.top')
        cmd='/bin/cp -rf '+equil0+' '+Td+'/ && '+\
             '/bin/cp '+self.preminPDB+' '+pdbf+' && '+\
             '/bin/cp '+self.crdfile+' '+crdf+' && '+\
             '/bin/cp '+self.parmfile+' '+parf
        os.system(cmd)
        #fill template submit script
        subf=joink(equil,'submit') #;Bye(subf)
        equilJaguar='/lustre/scr144/jmborr/dhfr/'+self.id+'/T_'+T+'/Equil'
        buf=''.join(open(subf,'r').readlines())
        buf=buf.replace('_JOBNAME_','equil_'+T)
        buf=buf.replace('_LOGDIR_',equilJaguar)
        buf=buf.replace('_WALLTIME_','1:00:00')
        buf=buf.replace('_NCORES_','128')
        buf=buf.replace('_EQUIL_',equilJaguar)
        buf=buf.replace('_PID_',self.id)
        open(subf,'w').write(buf) #replace subf with updated values
        #final preparations
        cmd='cd '+joink(equil,'files')+' && '+\
             './generate '+self.preminPDB+' Na && '+\
             './temperature.py '+T
        print '********** Manually execute the following commands ********'
        print cmd
        return False
        
                   
    def quasiharmonic(self):
        """obtain quasiharmonic models"""
        from amber.Prod import Prod
        self.Prod=Prod(self.id,joink(self.simdir,'Prod'))
        from amber.amber9 import ptraj
        traj=ptraj(self.unsolv['parmfile'],self.Prod.prev['crdlist'])
        traj.rms(outrms=rms.dat)['self']
        traj.correl()['self']
        traj.quasiharmonic()['self']
        traj.go()
        return False

    def genIsDeutd(self):
        """flag as deuterated"""
        self.isDeutd=False
        if self.id in ('dhp__D',): self.isDeutd=True
        return True

    def genInit(self,T='300'):
        """put all needed input files for the production run in jaguar
        in a single init directory"""
        Td=joink(self.simdir,'T_'+T)
        equil=joink(Td,'Equil')
        initd=joink(equil,'init') 
        os.system('/bin/mkdir -p '+initd)
        pairs={self.parmfile:joink(initd,self.id+'.top'),
               #self.unsolv['parmfile']:joink(initd,self.id+'.uns.top'),
               self.crdfile :joink(initd,self.id+'.ref'),
               #self.unsolv['crdfile']:joink(initd,self.id+'.uns.ref'),
               joink(equil,'md16',self.id+'.crd.md16'):joink(initd,self.id+'.rst')
               }
        for (a,b) in pairs.items(): os.system('/bin/cp '+a+' '+b)
        #create input file for pmemd
        buf=''.join(open(joink(equil,'files/md16.in'),'r').readlines())
        resline=re.compile('(RES \d+ \d+)').search(buf).group(1)
        buf="""
&cntrl
imin=0, irest=1, ntx=5,
ntt=0, tempi=_TEMP_ ,temp0=_TEMP_, tautp=1.0,
ntb=1, ntc=2, ntf=2, nsnb=10,
nstlim=200000, dt=0.001,
ntwe=1000, ntpr=1000, ntwx=20, ntwv=20
cut=10.0, ntr=1, nmropt=0,
&end
"""
        buf=buf.replace('_TEMP_',T)
        inf=joink(initd,self.id+'.in')
        open(inf,'w').write(buf)
        buf="""#Apply restraints to the solute atoms
 1.0
"""+resline+"""
END
END"""
        open(inf,'a').write(buf)
        cmd='cd '+equil+' && tar cf init.tar init' ;  os.system(cmd)
        return False

    def genInit_rubredoxin(self,T='300'):
        """put all needed input files for the production run in jaguar
        in a single init directory"""
        Td=joink(self.simdir,'T_'+T)
        equil=joink(Td,'Equil')
        initd=joink(equil,'init') 
        os.system('/bin/mkdir -p '+initd)
        pairs={self.parmfile:joink(initd,self.id+'.top'),
               self.unsolv['parmfile']:joink(initd,self.id+'.uns.top'),
               self.crdfile :joink(initd,self.id+'.ref'),
               self.unsolv['crdfile']:joink(initd,self.id+'.uns.ref'),
               joink(equil,'md16',self.id+'.crd.md16'):joink(initd,self.id+'.rst')
               }
        for (a,b) in pairs.items(): os.system('/bin/cp '+a+' '+b)
        #create input file for pmemd
        buf=''.join(open(joink(equil,'files/md16.in'),'r').readlines())
        resline=re.compile('(RES \d+ \d+)').search(buf).group(1)
        buf="""
&cntrl
imin=0, irest=1, ntx=5,
ntt=0, tempi=300.0 ,temp0=300.0, tautp=1.0,
ntb=1, ntc=2, ntf=2, nsnb=10,
nstlim=500000, dt=0.001,
ntwe=1000, ntpr=1000, ntwx=200, ntwv=200
cut=10.0, ntr=1, nmropt=1,
&end
"""
        inf=joink(initd,self.id+'.in')
        open(inf,'w').write(buf)
        os.system('/bin/cat '+self.NMR['amberdistf']+' >> '+inf)
        buf="""&wt TYPE='END' &end
#Apply restraints to the solute atoms
 1.0
"""+resline+"""
END
END"""
        open(inf,'a').write(buf)
        cmd='cd '+equil+' && tar cf init.tar init' ;  os.system(cmd)
        return False

    def checkInit(self,T='300'):
        """check init directory generated from Equil with a short
        simulation"""
        Td=joink(self.simdir,'T_'+T)
        equil=joink(Td,'Equil')
        initd=joink(equil,'init') 
        currd=os.getcwd()
        wd=joink(currd,junkName())
        print ' *** checking init in '+wd
        os.system('/bin/mkdir -p '+joink(wd,'curr'))
        cmd='/bin/cp -r '+initd+' '+wd ; os.system(cmd)
        os.system('cd '+wd+' && /bin/cp init/'+self.id+'.rst curr/')
        buf=''.join(open(joink(wd,'init',self.id+'.in'),'r').readlines())
        buf=re.compile( 'nstlim=\d+').sub( 'nstlim=2000', buf)
        open(joink(wd,'init',self.id+'.in'),'w').write(buf)
        cmd='cd '+wd+' && '+\
            'time sander -O -i init/'+self.id+'.in -p init/'+self.id+'.top -c curr/'+self.id+'.rst -o curr/'+self.id+'.out -x curr/'+self.id+'.crd -v curr/'+self.id+'.vel -ref init/'+self.id+'.ref'
        print ' *** submitting sander\n'+cmd
        if os.system(cmd): return False
        #check rms of unsolvated system
        buf= 'trajin curr/'+self.id+'.crd\n'+\
             'strip :WAT\nstrip :Na+\n'+\
             'trajout curr/'+self.id+'.uns.crd\n'+\
             'go\n'
        scriptf=joink(wd,junkName()) ; open(scriptf,'w').write(buf)
        os.system('cd '+wd+' && ptraj init/'+self.id+'.top < '+scriptf)
        rms=joink(wd,'rms')
        buf= 'trajin curr/'+self.id+'.uns.crd\n'+\
             'reference init/'+self.id+'.uns.ref\n'+\
             'rms reference out '+rms+'\n'+\
             'go\n'
        scriptf=joink(wd,junkName()) ; open(scriptf,'w').write(buf)
        os.system('cd '+wd+' && ptraj init/'+self.id+'.uns.top < '+scriptf)
        rmax=0.0
        for line in open(rms,'r').readlines():
            t,r=line.split() ; r=float(r)
            if r>rmax: rmax=r
        print 'max RMS='+str(rmax)+' for unsolvated, unsalted system'
        print 'Please REMOVE working directory '+wd
        return False

    def count_frames(self,T='300',simLength=''):
        """check number of frames correct"""
        from amber.Prod import Prod
        TT='T_'+T
        if TT not in self.solution.keys(): self.solution[TT]={}
        pr=Prod(self.id,joink(self.simdir,TT,simLength,'Prod'),
                module_name='amber.Prod')
        pr.count_frames(mode='unsolvated') #check number of frames
        return False
        
    def importProd(self,T='300',simLength=''):
        """import Prod structure directory"""
        from amber.Prod import Prod
        if 'solution' not in self.__dict__: self.solution={}
        TT='T_'+T
        if TT not in self.solution.keys(): self.solution[TT]={}
        pr=Prod(self.id,joink(self.simdir,TT,simLength,'Prod'),
                module_name='amber.Prod')
        pr.importStripTopCrd(self.unsolv['parmfile'],
                             self.unsolv['crdfile'])
        pr.count_frames(mode='unsolvated') #check number of frames
        #several RMS calculations to check consistency of simulation
        pr.getAvg()          #averages
        pr.rmsJump(outcrd=True)#consec rms and output fitted trajectory
        pr.rmsRef()        #rms with respect to reference
        pr.centr_and_avg() #find centroid
        pr.rmsCentr()      #rms with respect to centroid
        pr.rmsFirst()      #rms with respect to first frame
        pr.fit2centr()     #fit trajectory to centroid
        self.solution[TT]['Prod'+simLength]=pr #store
        return True


    def createMMTKtraj(self,T='300',recover=''):
        """create an MMTK trajectory from the consecutive trajectory"""
        import Numeric,pdb
        from MMTK import Units
        from amber.amber_to_nc import convertAmberToMMTK
        prmpdb=self.unsolv['preminPDB']
        Prod=self.solution['T_'+T]['Prod']
        if 'nMoldyn' not in self.solution['T_'+T].keys():
            self.solution['T_'+T]['nMoldyn']={}
        pM=self.solution['T_'+T]['nMoldyn']
        pM['simdir']=joink(self.simdir,'T_'+T,'nMoldyn')
        nc_file=basename(Prod.strip['fit2prev']).replace('.crd.gz','.nc')
        pM['nc']=joink(pM['simdir'],nc_file)
        if not isthere(pM['simdir']):
            os.system('mkdir -p '+pM['simdir'])
            
        prmcrd='junk347' ; prmvel='junk304'
            
        #prmcrd=junkName() ; prmvel=junkName() #temporary files
        #prmcrd0=joink(Prod.simdir,Prod.strip['fit2prev'])
        #print 'unzipping %s to %s'%(basename(prmcrd0),prmcrd)
        #cmd='gunzip -c %s > %s'%(prmcrd0,prmcrd)
        #os.system(cmd)
        #prmvel0=joink(Prod.simdir,Prod.strip['velfile'])
        #print 'unzipping %s to %s'%(basename(prmvel0),prmvel)
        #cmd='gunzip -c %s > %s'%(prmvel0,prmvel)
        #os.system(cmd)
        pM['pdbfixed']=joink(pM['simdir'],
                             basename(prmpdb).replace('.pdb','.fixed.pdb'))
        pM['pdb']=joink(pM['simdir'],
                        basename(prmpdb).replace('.pdb','.mmkt.pdb'))
        box=[float(x) for x in self.box.split(' ')]
        print 'converting to MMTK trajectory...'
        convertAmberToMMTK(prmpdb,prmcrd,pM['nc'],prmvel=prmvel,
                           format='crd',dtime=Prod.times['twx'],
                           box=box,
                           pdbfixed=pM['pdbfixed'],pdbtest=pM['pdb'])
        os.system('/bin/rm %s %s'%(prmcrd,prmvel))
        print '...converting finished'
        return True
        
    def get_mass(self):
        """retrieve mass(a.m.u) of the protein alone"""
        buf=''
        for line in open(self.preminPDB,'r').readlines():
            if line[0:3]=='TER': break
            buf+=line
        pdbf=junkName() ; open(pdbf,'w').write(buf)
        massf=junkName()
        cmd='pdb_to_mass '+pdbf+' '+massf ; os.system(cmd)
        mass=0.0
        for line in open(massf,'r').readlines():
            mass+=sum([ float(m) for m in line.split()])
        self.mass=mass
        return True

    def get_masses(self):
        """load atomic masses"""
        import numpy
        massf=junkName()
        os.system('pdb_to_mass '+self.preminPDB+' '+massf)
        m=[float(x) for x in open(massf).read().split()]        
        if self.isDeutd:
            for i in range(len(m)):
                if m[i]<1.02: m[i]=2.014
        self.masses=numpy.array(m)
        os.system('/bin/rm '+massf)
        return True
        
    def neededWaters(self,r=0.18):
        """print number of needed waters molecules for
        a given mass(H20)-to-mass(protein) ratio"""
        r=float(r)
        watermass=18.0 ; Nwat=int(self.mass/watermass)
        Nrwat=int(r*self.mass/watermass)
        if not 'powder' in self.__dict__.keys(): self.powder={}
        self.powder['nwat']=Nrwat
        return True

    def extractHfromPDB(self,pdb='init'):
        """extract hydrogen lines from preminPDB"""
        pdbf=self.pdbfile
        if pdb=='premin': pdbf=self.preminPDB
        buf=''
        for line in open(pdbf,'r').readlines():
            if line[0:3]=='TER': break
            if line[0:5]=='ATOM ' and line[13]=='H': buf+=line
        if pdb=='init':
            self.onlyHpdb=self.pdbfile.replace('.ent','.onlyH.ent')
            open(self.onlyHpdb,'w').write(buf)
        elif pdb=='premin':
            self.onlyHpreminPDB=self.preminPDB.replace('.pdb','.onlyH.pdb')
            open(self.onlyHpreminPDB,'w').write(buf)
        return True
        
    def getHinRes(self):
        """make a list of hydrogen numbers in every residue"""
        from copy import copy
        HinRes=[]
        l=[] ; currresi=1 ; N=0
        for line in open(self.onlyHpreminPDB,'r').readlines():
            resi=int(line[22:26])
            if resi==currresi:
                l.append(N)
            else: #new residue
                HinRes.append( copy(l) )
                l=[N,]
                currresi=resi
            N+=1
        HinRes.append( copy(l) )
        self.HinRes=HinRes
        return True

    def keepH(self,T='300',simLength=''):
        """generate trajectory containing only the hydrogens"""
        print 'generate trajectory containing only the hydrogens'
        from amber.amber10 import top,crd,ptraj
        self.solution['T_'+T]['onlyH'+simLength]={} #store stuff for H only
        pH=self.solution['T_'+T]['onlyH'+simLength]
        pH['simdir']=joink(self.simdir,'T_'+T,simLength,'onlyH')
        if not isthere(pH['simdir']): os.system('mkdir -p '+pH['simdir'])
        #check number of hydrogens read
        Nhydr=top(self.unsolv['parmfile']).numberAtoms(mask='@H*')
        pH['Nhydr']=Nhydr #store also here
        if 'Nhydr' in self.__dict__ and self.Nhydr!=Nhydr:
            Bye('different number of hydrogens than stored!')
        else:
            self.Nhydr=Nhydr
        Prod=self.solution['T_'+T]['Prod'+simLength] #Production object
        if 'crdfile' in pH.keys():
            pH['crd.fit2centr']=pH['crdfile']
            del pH['crdfile']
        #for keyword in ('fit2centr','fit2prev'):
        for keyword in ('fit2prev',):
            crdf=jk(Prod.simdir,Prod.strip[keyword])
            outf=basename(crdf.replace('.crd.gz','.H.crd'))
            outf=joink(pH['simdir'],outf)
            ptr=ptraj(self.unsolv['parmfile'],crdf)
            ptr.strip(entities=['!@H*',]) #keep only hydrogens
            ptr.go(outcrd=outf)
            print 'gzip H-only trajectory'
            os.system('gzip -f '+outf)
            pH['crd.'+keyword]=outf+'.gz'
        return True

    def keepHvel(self,T='300',simLength=''):
        """same as keepH but for velocities"""
        print 'generate velocity trajectory containing only the hydrogens'
        from amber.amber10 import top,crd,ptraj
        pH=self.solution['T_'+T]['onlyH'+simLength]
        #check number of hydrogens read
        Prod=self.solution['T_'+T]['Prod'+simLength] #Production object
        crdf=Prod.strip['velfile']
        outf=basename(crdf.replace('.vel.gz','.vel.H'))
        outf=joink(pH['simdir'],outf)
        ptr=ptraj(self.unsolv['parmfile'],crdf)
        ptr.strip(entities=['!@H*',]) #keep only hydrogens
        ptr.go(outcrd=outf)
        print 'gzip H-only velocity trajectory'
        os.system('gzip -f '+outf)
        pH['vel']=outf+'.gz'
        return True
        
    def Bfactors(self,T='300',mask=':H*',simLength=''):
        """calculate Bfactors using trajectory fit to centroid.
        Centroid is calculated using only CA atoms"""
        pr=self.solution['T_'+T]['Prod'+simLength]
        pr.Bfactors()
        pr.BfactorsByRes()
        pr.BfactorsByMask(mask='@H*')
        return True

    def av_msd2(self,T='300',simLength='',task='onlyH'):
        """average mean displacements"""
        from amber.amber10 import crd
        if task=='onlyH':
            onlyH_p=self.solution['T_'+T]['onlyH'+simLength]
            crdf=onlyH_p['crd.fit2prev']
            crdobj=crd(crdf,self.Nhydr)
            msd2s=crdobj.msd2()
            buf='#average msd^2 from trajectory with crd::msd2()\n'
            ix=1
            for x in msd2s:
                buf+='%4d %6.3f\n'%(ix,x); ix+=1
            outf=joink(onlyH_p['simdir'],'av_msd2.dat')
            onlyH_p['av_msd2']=outf
            open(outf,'w').write(buf)
            return True
        
    def av_H_msd(self,T='300',simLength='',pdb=''):
        """average mean square displacement of H or D from trajectory
        fitted to centroid"""    
        from math import sqrt,pi
        av=0.0 ; N=0
        if pdb=='init':
            for line in open(self.onlyHpdb,'r').readlines():
                av+=float(line[60:66]) ; N+=1
            print av/N
            msd=sqrt(av*3/(8*N))/pi ;print 'T=300 msd=',msd
            self.av_H_msd=msd
            return False
        Prod=self.solution['T_'+T]['Prod'+simLength]
        for line in open(Prod.strip['bfactmask'],'r').readlines():
            (i,b)=line.split()
            av+=float(b) ; N+=1
        msd=sqrt(av*3/(8*N))/pi ;print 'T='+T+' msd=',msd
        self.solution['T_'+T]['onlyH'+simLength]['av_H_msd']=msd
        return True

    def msdRatio(self,T='300',simLength=''):
        """ratio between computed and calculated m.s.d"""
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        ratio=[1/2.2]*self.Nhydr #default ratio
        Prod.strip['bfactmask']
        #pH[]
        return False
        
    def avPerResHmsd(self,T='300',simLength='',pdb=''):
        """average mean square displacement of H or D from trajectory
        fitted to centroid, per residue"""        
        from math import sqrt,pi
        if pdb=='init': #from experiment
            for l in open(self.onlyHpdb,'r').readlines():
                ires=int(l[23:26])-1
                bf[ires]+=float(l[60:66])
                nbf[ires]+=1
            buf=''
            for ires in range(self.L):
                bf[ires]=sqrt(3*bf[ires]/(8*nbf[ires]))/pi
                buf+=str(ires)+' '+str(bf[ires])+'\n'
            self.av_perRes_H_msd=joink(self.simdir,'av_perRes_H_msd.dat')
            open(self.av_perRes_H_msd,'w').write(buf)
            return True
        else:
            Prod=self.solution['T_'+T]['Prod'+simLength]
            bf=[]
            for line in open(Prod.strip['bfactmask'],'r').readlines():
                (i,b)=line.split()  ;  bf.append(float(b))
            msd=[0.0]*self.L ; buf=''
            for ires in range(self.L):
                Hlist=self.HinRes[ires]
                for iH in Hlist: msd[ires]+=bf[iH]
                msd[ires]=sqrt(3*msd[ires]/(8*len(Hlist)))/pi
                buf+=str(ires)+' '+str(msd[ires])+'\n'
            pH=self.solution['T_'+T]['onlyH'+simLength]
            pH['msd_perRes']=joink(pH['simdir'],'msd_perRes.dat')
            open(pH['msd_perRes'],'w').write(buf)
            return True

    def getSIIScF(self,nt=200,T='300',scale='log10',simLength=''):
        """calculate Self-Intermediate Incoherent Scattering Function
        nt:number of time points for which to calculate the function"""
        from scatter.scatter import SIIScF
        from amber.amber10 import crd
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        crdobj.loadframes(skip=skipframes)
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']
        tl=(crdobj.nframe-10)*Prod.times['twx']
        qq=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] #Angstroms^(-1)
        scattering=SIIScF(qq,crdobj.frames,dt,tl,nt,scale=scale)
        buf=''
        for iq in range(len(qq)):
            q=scattering[iq][0]
            buf+='q='+str(q)+'\n'
            for (t,s) in scattering[iq][1]: buf+=str(t)+' '+str(s)+'\n'
        Bye(buf)
        pH['SIIScF'+scale]=crdf.replace('.crd.gz','.SIIScF.'+scale+'.dat')
        open(pH['SIIScF'+scale],'w').write(buf)
        return True

    def getSIIScF_sts(self,nt=200,T='300',scale='linear',simLength=''):
        """calculate Self-Intermediate Incoherent Scattering Function
        in the short time scale (0.2psec to 50 psec)
        nt:number of time points for which to calculate the function"""
        from scatter.scatter import SIIScF
        from amber.amber10 import crd
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']
        tl=100 #100 picoseconds
        nt=int(tl/dt)
        qq=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] #Angstroms^(-1)
        crdobj.loadframes(skip=skipframes)
        scattering=SIIScF(qq,crdobj.frames,dt,tl,nt,scale=scale)
        buf=''
        for iq in range(len(qq)):
            q=scattering[iq][0]
            buf+='q='+str(q)+'\n'
            for (t,s) in scattering[iq][1]: buf+=str(t)+' '+str(s)+'\n'
        pH['SIIScF_sts']=crdf.replace('.crd.gz','.SIIScF_sts.dat')
        open(pH['SIIScF_sts'],'w').write(buf)
        return True

    def getVelCorr_sts(self,nt=500,T='300',scale='linear',simLength=''):
        """calculate velocity autocorrelation function
        in the short time scale (0.2psec to 50 psec)
        nt:number of time points for which to calculate the function"""
        from scatter.scatter import autoCorr
        from amber.amber10 import crd
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['vel'] #velocities of hydrogen
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']
        tl=100 #100 picoseconds
        nt=int(tl/dt)
        crdobj.loadframes(skip=skipframes)
        scattering=autoCorr(crdobj.frames,dt,tl,nt,scale=scale)
        buf=''
        for (t,s) in scattering: buf+=str(t)+' '+str(s)+'\n'
        pH['velACorr_sts']=crdf.replace('.vel.H.gz','.velACorr_sts.dat')
        open(pH['velACorr_sts'],'w').write(buf)
        return True

    def getDensStates(self,T='300',simLength=''):
        """fourier transform of velACorr_sts"""
        from scatter.scatter import fourier_abs
        Prod=self.solution['T_'+T]['Prod'+simLength]
        dt=Prod.times['twx'] #space in between conformations
        factor=dt/4.14 #from picoseconds to 10E-6eV)
        pH=self.solution['T_'+T]['onlyH'+simLength]
        pt=open(pH['velACorr_sts'],'r')
        tt=[] ; ss=[]
        line=pt.readline()
        while line:
            (t,s)=line.split()
            tt.append(float(t)) ; ss.append(float(s))
            line=pt.readline()
        ss=fourier_abs(ss) ; N=1 ; buf=''
        for x in ss:
            buf+=str(N*factor)+' '+str(x)+'\n'
            N+=1
        pH['DenSt_sts']=pH['velACorr_sts'].replace('velACorr','DenSt')
        open(pH['DenSt_sts'],'w').write(buf) ; print pH['DenSt_sts']
        return True
        
    def getSq_sts(self,T='300',conv=None,simLength=''):
        """fourier transform SIIScF_sts

        conv: a Gaussian to multiply I(q,t)"""
        from scatter.scatter import fourier_abs
        Prod=self.solution['T_'+T]['Prod'+simLength]
        dt=Prod.times['twx'] #space in between conformations
        factor=dt/4.14 #from picoseconds to 10E-6eV)
        pH=self.solution['T_'+T]['onlyH'+simLength]
        pt=open(pH['SIIScF_sts'],'r')
        pH['Sq_sts']=pH['SIIScF_sts'].replace('SIIScF','Sq')
        buf=''
        qline=pt.readline() ; line=pt.readline()
        while line:
            buf+=qline
            tt=[] ; ss=[]
            while line and 'q' not in line:
                (t,s)=line.split()
                tt.append(float(t)) ; ss.append(float(s))
                line=pt.readline()
            qline=line
            if conv:
                for it in range(len(tt)):
                    t=tt[it]
                    ss[it]=ss[it]*conv(t)
            ss=fourier_abs(ss)
            N=2
            for x in ss[1:]: #avoid first entry, which is delta function
                buf+=str(N*factor)+' '+str(x)+'\n'
                N+=1
            line=pt.readline()
        pH['Sq_sts']=pH['SIIScF_sts'].replace('SIIScF','Sq')
        open(pH['Sq_sts'],'w').write(buf)
        return True
        
    def get_msd(self,nt=200,T='300',scale='log10',simLength=''):
        """calculate time-dependent average atomic fluctuation
        nt:number of time points for which to calculate the function"""
        from scatter.scatter import msd
        from amber.amber10 import crd
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        crdobj.loadframes(skip=skipframes)
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']
        tl=(crdobj.nframe-10)*Prod.times['twx']
        scattering=msd(crdobj.frames,dt,tl,nt,scale=scale)
        pH['msd'+scale]=crdf.replace('.crd.gz','.msd.'+scale+'.dat')
        pt=open(pH['msd'+scale],'w')
        for (t,s) in scattering: pt.write(str(t)+' '+str(s)+'\n')
        pt.close()
        return True


    def getEISF(self,T='300',simLength=''):
        """calculate Self-Intermediate Incoherent Scattering Function
        nt:number of time points for which to calculate the function"""
        from scatter.scatter import EISF
        from amber.amber10 import crd
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        centNframe=Prod.strip['centNframe']
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        crdobj.loadframes(skip=skipframes)
        qq=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] #Angstroms^(-1)
        scattering=EISF(qq,crdobj.frames)
        buf=''
        for (q,s) in scattering: buf+=str(q)+' '+str(s)+'\n'
        Bye(buf)
        pH['EISF']=crdf.replace('.crd.gz','.EISF.dat') #; print pH['EISF']
        open(pH['EISF'],'w').write(buf)
        return True

    def getEISFperRes(self,T='300',simLength=''):
        """calculate Self-Intermediate Incoherent Scattering Function
        per residue (average of the hydrogens within a single residue)
        for different scattering vectors"""
        from scatter.scatter import EISF
        from amber.amber10 import crd
        from numpy import rank
        Prod=self.solution['T_'+T]['Prod'+simLength]
        centNframe=Prod.strip['centNframe'] - 1      
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        crdobj.loadframes(skip=skipframes)
        qq=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
        scattering=[]
        for q in qq: scattering.append('#q='+str(q)+'\n')
        for resi in range(self.L):
            #extract the trajectory of the hydrogens in residue "resi"
            a=self.HinRes[resi][0] ; b=1+self.HinRes[resi][-1]
            frs=(crdobj.frames).transpose((1,0,2))[a:b].transpose((1,0,2))
            z=EISF(qq,frs,centNframe)
            for iq in range(len(qq)):
                (q,s)=z[iq]
                scattering[iq]+=str(resi+1)+' '+str(s)+'\n'
        pH['EISFperRes']=crdf.replace('.crd.gz','.EISFperRes.dat')
        print pH['EISFperRes']
        open(pH['EISFperRes'],'w').write('\n'.join(scattering))
        return True

    def rmsCentr(self,T='300',simLength=''):
        from amber.Prod import Prod
        TT='T_'+T
        Prod=self.solution['T_'+T]['Prod'+simLength]
        Prod.rmsCentr()
        return True

    def Ktraj(self,T='270',subdir=None):
        """create trajectory of kinetic energy"""
        from amber.amber10 import crd
        from pdbfiles.smallThings import insertBfact
        from utilities.readingWritingFiles import write_from_numpy
        Prod=self.solution['T_'+T]['Prod']
        if subdir: Prod=self.__dict__['T_'+T][subdir]['Prod']
        Kcrd=joink(Prod.simdir,'strip','K.crd');
        if isthere(Kcrd): os.system('/bin/rm '+Kcrd)
        Kpt=open(Kcrd,'a')
        velf=joink(Prod.simdir,Prod.strip['velfile'])
        nat=Prod.strip['nat']; crdobj=crd(velf,nat)
        Kav,Kdev=crdobj.kinetic_energy(self.masses[0:nat],outptf=Kpt)
        #gzip trajectory
        os.system('gzip -f '+Kcrd)
        Prod.strip['Kcrd']=joink('strip','K.crd.gz')
        #save average Kinetic energy as .pdb
        buf=insertBfact(self.unsolv['preminPDB'],Kav)
        Prod.strip['Kavpdb']=joink('strip','Kav.pdb')
        open(jk(Prod.simdir,Prod.strip['Kavpdb']),'w').write(buf)
        #save average Kinetic energy as .dat
        buf='\n'.join([str(x) for x in Kav])
        Prod.strip['Kavdat']=joink('strip','Kav.dat')
        open(jk(Prod.simdir,Prod.strip['Kavdat']),'w').write(buf)
        Prod.dump()
        return True


#################### end of PEPSICAN  ##################
    
    
def datatree():
    """print the datatree of one varset (a PEPSI(CAN) object)

    ojbect ( a PEPSI(CAN) object )
    |_'id': string identifying the CAN object
    |_'dumpf': picke dump file
    |_'dumpbf': basename of pickle dump file
    |_'dumpdf': dirname of pickle dump file
    |_'simdir': directory to store simulation results
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


def initCollection(*kargs,**kwargs):

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
    idlist=canlist
    if 'id' in kwargs.keys(): idlist=[kwargs['id'],]
    remaining=len(idlist)
    for id in idlist:
        c=PEPSI(id=id)
        c.genDumpf(cl.repod)
        c.pickleDump()
        print remaining,id
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
    
def gensimdir(*kargs,**kwargs):
    """simulation directory"""
    cl.cycleOverList('gensimdir',*kargs,**kwargs)

def getinitPDB(*kargs,**kwargs):
    """retrieve PDB file"""
    cl.cycleOverList('getinitPDB',*kargs,**kwargs)

def parsePDB(*kargs,**kwargs):
    """parse PDB file for tleap use

    (1)remove non ATOM lines
    (2)remove other chains than ' ' or 'A'
    (3)remove hydrogen and 'D' atoms"""
    cl.cycleOverList('parsePDB',*kargs,**kwargs)

def storeL(*kargs,**kwargs):
    """store protein length"""
    cl.cycleOverList('storeL',*kargs,**kwargs)

def get_bounding_box(*kargs,**kwargs):
    """extract the bounding box from the solvated system"""
    cl.cycleOverList('get_bounding_box',*kargs,**kwargs)
    
def gen_tleap_in(*kargs,**kwargs):
    """generate input tleap"""
    cl.cycleOverList('gen_tleap_in',*kargs,**kwargs)

def getcharge(*kargs,**kwargs):
    """find net charge of parsed pdb file"""
    cl.cycleOverList('getcharge',*kargs,**kwargs)

def genparmtop(*kargs,**kwargs):
    """generate paramenter and topology files"""
    cl.cycleOverList('genparmtop',*kargs,**kwargs)

def genPackageMin(*kargs,**kwargs):
    """package input for minimization"""
    cl.cycleOverList('genPackageMin',*kargs,**kwargs)

def genPackageForMin(*kargs,**kwargs):
    """package input for minimization"""
    cl.cycleOverList('genPackageForMin',*kargs,**kwargs)

def genPDBid(*kargs,**kwargs):
    """initialize the PDB id entry"""
    cl.cycleOverList('genPDBid',*kargs,**kwargs)

def createNMRdist(*kargs,**kwargs):
     cl.cycleOverList('createNMRdist',*kargs,**kwargs)

def createFeLib():
    """create library for Fe+3 atom"""   
    tleap_in=junkName()
    open(tleap_in,'w').write("""FE3=createUnit FE3
res_fe3=createResidue \"FE3\"
atom_fe3=createAtom FE FE 3
set atom_fe3 element Fe
set atom_fe3 type FE
add res_fe3 atom_fe3
add FE3 res_fe3
saveOff FE3 FE3.lib
quit""")
    os.system('tleap -sf '+tleap_in+' && /bin/rm '+tleap_in)
    print 'DO: sudo mv FE3.lib $AMBERHOME/dat/leap/lib'
    frcmod='FE.frcmod'
    open(frcmod,'w').write("""#parameters for FE
MASS
FE 55.845

NONBON
FE   1.25000   0.40000""")
    print 'DO: sudo mv '+frcmod+' $AMBERHOME/dat/contrib'

def setEquilInput(*kargs,**kwargs):
    """import certain file names neccessary for equilibrium"""
    cl.cycleOverList('setEquilInput',*kargs,**kwargs)

def relax(*kargs,**kwargs):
    """minimize and relax the premin pdb structure"""
    cl.cycleOverList('relax',*kargs,**kwargs)

def quasiharmonic(*kargs,**kwargs):
    """perform quasiharmonic analysis"""
    cl.cycleOverList('quasiharmonic',*kargs,**kwargs)

def deuteratePDB(*kargs,**kwargs):
    """modify preminPDB by hydrogen-deuterium exchange and rewrite
    to parsedPDb"""
    cl.cycleOverList('deuteratePDB',*kargs,**kwargs)

def genIsDeutd(*kargs,**kwargs):
    """flag as deuterated"""
    cl.cycleOverList('genIsDeutd',*kargs,**kwargs)

def genInit(*kargs,**kwargs):
    """put all needed input files for the production run in jaguar
    in a single init directory"""
    cl.cycleOverList('genInit',*kargs,**kwargs)

def checkInit(*kargs,**kwargs):
    """check init directory generated from Equil with a short
    simulation"""
    cl.cycleOverList('checkInit',*kargs,**kwargs)

def get_mass(*kargs,**kwargs):
    """retrieve mass of the protein alone"""
    cl.cycleOverList('get_mass',*kargs,**kwargs)

def get_masses(*kargs,**kwargs):
    """retrieve mass of the protein alone"""
    cl.cycleOverList('get_masses',*kargs,**kwargs)

def neededWaters(*kargs,**kwargs):
    """print number of needed waters molecules for
    a given mass(H20)-to-mass(protein) ratio"""
    cl.cycleOverList('neededWaters',*kargs,**kwargs)

def avPerResHmsd(*kargs,**kwargs):
    """average mean square displacement of H or D from trajectory
    fitted to centroid, per residue"""        
    cl.cycleOverList('avPerResHmsd',*kargs,**kwargs)
    
def getSIIScF(*kargs,**kwargs):
    """calculate different scattering-related functions"""
    cl.cycleOverList('getSIIScF',*kargs,**kwargs)

def getSIIScF_sts(*kargs,**kwargs):
    """calculate different scattering-related functions"""
    cl.cycleOverList('getSIIScF_sts',*kargs,**kwargs)

def getVelCorr_sts(*kargs,**kwargs):
    """calculate velocity autocorrelation function"""
    cl.cycleOverList('getVelCorr_sts',*kargs,**kwargs)

def getDensStates(*kargs,**kwargs):
    """fourier transform of velACorr_sts"""
    cl.cycleOverList('getDensStates',*kargs,**kwargs)

def getSq_sts(*kargs,**kwargs):
    """fourier transform SIIScF_sts"""
    cl.cycleOverList('getSq_sts',*kargs,**kwargs)

def count_frames(*kargs,**kwargs):
    """import Prod structure directory"""
    cl.cycleOverList('count_frames',*kargs,**kwargs)

def importProd(*kargs,**kwargs):
    """import Prod structure directory"""
    cl.cycleOverList('importProd',*kargs,**kwargs)

def extractHfromPDB(*kargs,**kwargs):
    """extract hydrogen lines from preminPDB"""
    cl.cycleOverList('extractHfromPDB',*kargs,**kwargs)

def getHinRes(*kargs,**kwargs):
    """make a list of hydrogen numbers in every residue"""
    cl.cycleOverList('getHinRes',*kargs,**kwargs)

def keepH(*kargs,**kwargs):
    """generate one trajectory containing only the hydrogens"""
    cl.cycleOverList('keepH',*kargs,**kwargs)

def keepHvel(*kargs,**kwargs):
    """same as keepH but for velocities"""
    cl.cycleOverList('keepHvel',*kargs,**kwargs)

def Bfactors(*kargs,**kwargs):
    """calculate Bfactors using trajectory fit to centroid.
    Centroid is calculated using only CA atoms"""
    cl.cycleOverList('Bfactors',*kargs,**kwargs)

def av_msd2(*kargs,**kwargs):
    cl.cycleOverList('av_msd2',*kargs,**kwargs)

def getResolution(*kargs,**kwargs):
    """create gaussian for resolution
    The original data was taken at 32 10-6, the BASIS data was
    collected at 3 10-6, these are fwhm resolutions"""
    fwhm={'original':32E-6,'BASIS':3E-6}
    z=2.354820045
    sigma={'original':z*32E-6,'BASIS':z*3E-6}

def get_msd(*kargs,**kwargs):
    """calculate Self-Intermediate Incoherent Scattering Function
    nt:number of time points for which to calculate the function"""
    cl.cycleOverList('get_msd',*kargs,**kwargs)

def av_H_msd(*kargs,**kwargs):
    """average mean square displacement of H or D from trajectory
    fitted to centroid"""
    cl.cycleOverList('av_H_msd',*kargs,**kwargs)

def Ktraj(*kargs,**kwargs):
    """create trajectory of kinetic energy"""
    cl.cycleOverList('Ktraj',*kargs,**kwargs)
    
def rmsCentr(*kargs,**kwargs):
    cl.cycleOverList('rmsCentr',*kargs,**kwargs)

def getEISF(*kargs,**kwargs):
    """calculate Self-Intermediate Incoherent Scattering Function"""
    cl.cycleOverList('getEISF',*kargs,**kwargs)

def getEISFperRes(*kargs,**kwargs):
    """calculate Self-Intermediate Incoherent Scattering Function per
    residue (average of the hydrogens within a single residue) for
    different scattering vectors"""
    cl.cycleOverList('getEISFperRes',*kargs,**kwargs)

def createMMTKtraj(*kargs,**kwargs):
    """create an MMTK trajectory from the consecutive trajectory"""
    cl.cycleOverList('createMMTKtraj',*kargs,**kwargs)

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
