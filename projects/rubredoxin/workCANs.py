#!/usr/bin/python

import pdb
import os,sys,re,numpy
from tempfile import mkstemp
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir,projectsdir
from utilities.small_utilities import chomp,Bye,junkName
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
pwrite=sys.stdout.write
joink=os.path.join; jk=os.path.join
isthere=os.path.exists
basename=os.path.basename

my_project='rubredoxin'                   #substitute 'my_project'
pcd=joink(codedir,'projects',my_project)  #code directory
pjd=joink(projectsdir,my_project)         #projects directory
pjdout=joink(projectsdir,my_project,'out')
scd=joink(scratchdir,my_project,'out')
collf=joink(pjd,'varsets.can')            #you may substitute scd for pjd
canlistf=joink(pjd,'out','list')          #you may substitute scd for pjd
repod=joink(pjd,'repod')                  #you may substitute scd for pjd

def distrE(E,type='BASIS'):
    """gaussian distributions for energy resolution
    energy units are 1E-6eV"""
    from math import exp,pi
    #fwhm={'corrected':24,'original':32,'BASIS':3}
    z=1.0/2.354820045
    sigma={'original':z*32,'BASIS':z*3}
    ss=sigma[type] ; ss=ss*ss
    return exp(-E*E/(2*ss))/sqrt(2*pi*ss)

def convT(t,type='BASIS',ss=None):
    """fourier transform of distrE
    Time units are in picoseconds"""
    from math import exp,pi,sqrt
     #time resolution windows, in picoseconds
    sigma={'dcs-nist':108.8,'original':49.0,'BASIS':521.0,
           'corrected':65}
    if not ss: ss=sigma[type]
    ss=2*ss*ss ; f=1/sqrt(2*pi*ss)
    return f*exp(-t*t/ss)
            
def convT_BASIS(t):
    return convT(t,type='BASIS')

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
        Prod=self.solution['T_300']['Prod']
        Prod.strip['fit2centr.evecsall.project']='strip/evecs.all.project.dat.gz'
        return True#False


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
            if self.id in ('ref15-shel','ref15-shel_D'):
                self.pdbfile=joink(self.simdir,'pdbref15-shel.ent')
            print self.pdbfile
            return True
        
        self.pdbfile=PDBList().retrieve_pdb_file(self.id[0:4].upper(),pdir=self.simdir)
        return True

    def storeL(self):
        """store protein length"""
        self.L=53
        if self.id=='awq': self.L=171
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
            if p[0:4]!='ATOM' and 'FE' not in p[17:20]: continue
            if p[13] in ('H','D'): continue
            if p[14:16]=='DZ': continue
            if self.id in ('1vcx_Fe','1vcx_Fe_D','1fhm_D','1fhm',
                           'ref15-shel','ref15-shel_D'):
                if 'CYS' in p: p=p.replace('CYS','CYM') #CYS-->CYM
                if ' FE   FE   ' in p:
                    p=p.replace(' FE   FE   ',' FE   FE3  ')
                elif 'FE    FE' in p:
                    p=p.replace('FE    FE','FE   FE3')
                elif 'FE   FE2 A  55' in p:
                    p=p.replace('FE   FE2 A  55','FE   FE3 A  54')
            buf+=p
        buf+='TER\nEND\n'
        if self.id in ('1vcx_Fe','1vcx_Fe_D','1fhm','1fhm_D',
                       'ref15-shel','ref15-shel_D'):
            buf=buf.split('\n')
            buf='\n'.join( buf[0:-4]+['TER',]+buf[-4:] )
        self.pdbparsed=joink(self.simdir,'parsed_'+basename(self.pdbfile))
        open(self.pdbparsed,'w').write(buf)
        return True

    def genIsDeutd(self):
        """flag as deuterated"""
        self.isDeutd=False
        if self.id in ('1vcx_Fe_D','1fhm_D','ref15-shel_D'):
            self.isDeutd=True
        return True

    def getcharge(self):
        """find net charge of parsed pdb file"""
        buf='source leaprc.ff98\n'
        if self.isDeutd: buf='source leaprc.ff98Dsimple\n'
        if self.id in ('1vcx_Fe','1vcx_Fe_D','1fhm','1fhm_D',
                       'ref15-shel','ref15-shel_D'):
            buf+='loadOff FE3.lib\n'
        buf+=self.id+'=loadPdb '+self.pdbparsed+'\n'
        buf+='charge '+self.id+'\n'
        buf+='quit\n'
        script=junkName() ; open(script,'w').write(buf)
        tlout=junkName()
        os.system('tleap9 -f '+script+' > '+tlout)
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
        (1)add missing hydrogens"""        
        self.intleap=''
        buf='source leaprc.ff98\n'
        if self.isDeutd: buf='source leaprc.ff98Dsimple\n'
        if self.id in ('1vcx_Fe','1vcx_Fe_D','1fhm','1fhm_D',
                       'ref15-shel','ref15-shel_D'):
            buf+='loadOff FE3.lib\n'
            buf+='loadamberparams FE.frcmod\n'
        buf+=self.id+'=loadPdb '+self.pdbparsed+'\n'
        if mode=='solvated':
            """(2)solvate box with 10 Angstroms at least around the protein
            (3)add counterions"""
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
            pdb.set_trace()
        return True

    def genparmtop(self,mode='solvated',recover=''):
        """generate paramenter, topology, and preminPDB files"""
        if mode=='solvated':
            if not self.intleap: sys.stderr.write('no in.tleap file')
            self.preminPDB=joink(self.simdir,self.id+'.premin.pdb')
            self.outtleap=joink(self.simdir,'out.tleap')
            if recover: return True
            os.system('tleap9 -f '+self.intleap+' &> '+self.outtleap)
            #create pre-minimization PDB file
            os.system('ambpdb -p '+self.parmfile+' < '+self.crdfile+\
                      ' > '+self.preminPDB)
        elif mode=='unsolvated':
            pU=self.unsolv
            intleap=pU['intleap']
            pU['outtleap']=joink(self.simdir,'out.uns.tleap')
            outtleap=pU['outtleap']
            if not intleap: sys.stderr.write('no in.tleap file')
            outtleap=joink(self.simdir,'out.tleap')
            if recover: return True
            os.system('tleap9 -f '+intleap+' &> '+outtleap)
            #create pre-minimization PDB file
            pU['preminPDB']=joink(self.simdir,self.id+'.uns.premin.pdb')
            os.system('ambpdb -p '+pU['parmfile']+' < '+pU['crdfile']+\
                      ' > '+pU['preminPDB'])
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
        self.pdbid=self.id[0:-2]
        return True

    def workDSSP(self,task='',T='300',simLength='10ns'):
        """appy DSSP and related"""
        if task=='read and HB map':
            from secondary.dssp import dssp
            from tempfile import mkstemp
            pdbfile=self.unsolv['preminPDB']
            (a,dsspfile)=mkstemp()
            cmd='dsspcmbi.x %s %s'%(pdbfile,dsspfile)
            os.system(cmd)
            dobj=dssp(dsspfile)
            dobj.HBmap()
            self.dssp={'dobj':dobj}
            os.system('/bin/rm '+dsspfile)
            return True
        elif task=='average HBmap':
            from amber.Prod import Prod
            pr=self.solution['T_'+T]['Prod'+simLength]

            #produce 1000 PDB's from trajectory
            #from amber.amber10 import ptraj
            #from tempfile import mkdtemp
            #npdb=1000; nframe=int(pr.prev['last'])
            #delta=int(nframe/npdb)
            #crdf=joink(pr.simdir,pr.strip['fit2centr'])
            #ptr=ptraj(self.unsolv['parmfile'],crdf)
            #ptr.trajin(start=1,end=nframe,offset=delta)
            #wd=mkdtemp(dir='/projects/tmp')
            #ptr.trajout(joink(wd,'snap'),args='pdb')
            #pdb.set_trace()
            #ptr.go()

            #apply DSSP to each PDB

            wd='/projects/tmp/tmpt8H6XK'; nframe=5000; delta=50
            
            from secondary.dssp import dssp
            from tempfile import mkstemp
            (a,dsspfile)=mkstemp()
            pdb.set_trace()
            Nu=self.dssp['dobj'].nres
            avhbm=numpy.zeros(N*N).reshape(N,N)
            for isnap in range(1,nframe,delta):
                pdbfile=joink(wd,'snap.'+str(isnap))
                os.system( 'dsspcmbi.x %s %s'%(pdbfile,dsspfile) )
                dobj=dssp(dsspfile)

                dobj.HBmap()
            os.system('/bin/rm '+dsspfile)
            #average the DSSP files

            #remove the PDB and DSS files

            #create average HBmap

            #clean-up
            os.system('/bin/rm -rf wd')
            return False #True
        return False
            
    def createNMRdist(self,recover=''):
        """create NRM distance file according to Pratul"""
        templ=""" &rst
 iat=_AT1_,_AT2_,
    iresid=0,
         r1=0.0, r2=2.2, r3=2.3, r4=999.0,
         rk2=25.0, rk3=25.0,
 &end
"""
        if self.id in ('1vcx_Fe','1vcx_Fe_D','1fhm','1fhm_D',
                       'ref15-shel','ref15-shel_D'):
            self.NMR={} ; p=self.NMR
            p['amberdistf']=joink(self.simdir,self.id+'.amber.dist')
            if recover: return True
            at1s=('82','133','560','603') ; at2='789'
            if self.id in ('1fhm','1fhm_D',):
                at1s=('106','146','573','616') ; at2='782'
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
       
    def relax(self,T='300'):
        """minimize and relax the premin pdb structure"""
        #copy template Equil directory and a few input files
        equil0=os.environ['AMBERHOME']+'/minprotoc/Equil'
        if self.id in ('1vcx_Fe','1vcx_Fe_D','1fhm','1fhm_D','ref15-shel',
                       'ref15-shel_D'):
            equil0=os.environ['AMBERHOME']+'/minprotoc/nmr/Equil'
        Td=joink(self.simdir,'T_'+T) ; os.system('/bin/mkdir -p '+Td)
        equil=joink(Td,basename(equil0)) #; Bye(equil)
        pdbf=joink(equil,'srt',self.id+'.pdb')
        crdf=joink(equil,'srt',self.id+'.crd.srt')
        parf=joink(equil,'srt',self.id+'.top')
        RST=open(self.NMR['amberdistf']).read()
        cmd='/bin/cp -rf '+equil0+' '+Td+'/ && '+\
             '/bin/cp '+self.preminPDB+' '+pdbf+' && '+\
             '/bin/cp '+self.crdfile+' '+crdf+' && '+\
             '/bin/cp '+self.parmfile+' '+parf
        os.system(cmd)
        #fill template submit script
        subf=joink(equil,'submit') #;Bye(subf)
        equilJaguar='/lustre/scr144/jmborr/rubredoxin/'+\
                     self.id+'/T_'+T+'/Equil'
        buf=''.join(open(subf,'r').readlines())
        buf=buf.replace('_JOBNAME_','equil_'+T)
        buf=buf.replace('_LOGDIR_',equilJaguar)
        buf=buf.replace('_WALLTIME_','2:00:00')
        buf=buf.replace('_NCORES_','64')
        buf=buf.replace('_EQUIL_',equilJaguar)
        buf=buf.replace('_PID_',self.id)
        open(subf,'w').write(buf) #replace subf with updated values
        #fill md14,md15,md16 template files with restrains info
        for run in 'md14','md15','md16':
            templatef=joink(equil,'files',run+'.in.template')
            buf=open(templatef).read()
            buf=buf.replace('_APPEND_RST_FILE_HERE_\n',RST)
            open(templatef,'w').write(buf)
        #final preparations
        cmd='cd '+joink(equil,'files')+' && '+\
             './generate '+self.preminPDB+' Na && '+\
             './temperature.py '+T
        print "****** MANUALLY execute the following commands *******"
        print cmd
        return False

    def are_orthogonal(self,T='300'):
        """check if quasiharmonic modes are orthogonal"""
        from amber.Prod import Prod
        from amber.amber10 import quasiharmonic_mode
        pr=self.solution['T_'+T]['Prod']
        pt=open(joink(pr.simdir,pr.strip['fit2centr.evecsall']))
        q1=quasiharmonic_mode(pt)
        q2=quasiharmonic_mode(pt)
        print q1.massmul(q2,self.masses)
        print q1*q2#/(q1.norm()*q2.norm())
        print q1.massmul(q1,self.masses)
        print q1*q1
        return False
    
    def quasiharmonic(self,T='300',nvecs=1000,simLength='',task=''):
        from amber.Prod import Prod
        from amber.amber10 import ptraj
        pr=self.solution['T_'+T]['Prod'+simLength]
        if task=='project modes':
            """use Amber tools to project modes onto trajectory"""
            sc="""trajin _TRAJ_
projection modes evecs.all.dat out _PROJF_ beg 1 end _NVECS_ start 1 stop 50000 offset 50"""
            dict={'_TRAJ_':self.id+'.00001_50000.fit2centr.crd',
                  '_PROJF_': self.id+'.00001_50000_50_modes.fit2centr.crd',
                  '_NVECS_' : str(nvecs)}
            for k,v in dict.items(): sc=sc.replace(k,v)
            a,scf=mkstemp(dir='/projects/tmp')
            open(scf,'w').write(sc)
            cmd='cd %s && gunzip %s && ptraj %s < %s && '\
                 'gzip %s && /bin/rm %s'%\
                 (jk(pr.simdir,'strip'),
                  self.id+'.00001_50000.fit2centr.crd.gz',
                  self.unsolv['parmfile'],
                  scf,
                  self.id+'.00001_50000.fit2centr.crd',
                  scf
                  )
            os.system(cmd)
            return False
        elif task=='find msd^2(t)':
            """find msd^2(t) for each mode"""
            from amber.amber10 import evecs_projection
            from scatter.scatterModes import msd2
            """Calculate rmsd^2(t) from file projection.all.crd.gz"""
            projf=jk(pr.simdir,'strip/projection.all.crd.gz')
            proj=evecs_projection(projf)
            rr=numpy.array(proj[:])
            dt=0.02        #time in between frames (in picoseconds)
            tl=40.0        #limit time, in picoseconds
            nt=int(tl/dt)  #number of msd2 points to compute
            scattering=msd2(rr,dt,tl,nt)
            #store msd2 for each mode
            dir=jk(pr.simdir,'strip/evecs_msd2')
            os.system('/bin/mkdir -p '+dir)
            for i in range(proj.nmodes):
                id='%03d'%(i+1,);
                buf='#msd2 for mode '+id+'\n# t msd2\n'
                t=0
                for x in scattering[i]:
                    buf+='%4.2f %6.3f\n'%(t,x); t+=dt
                open(jk(dir,id+'.dat'),'w').write(buf)
            #store msd2 from average for all modes
            buf= '#msd2 for modes %3d to %3d\n'%(1,proj.nmodes)
            buf+='#time(psec) msd2\n'
            time=0
            for x in scattering.sum(axis=0):
                buf+='%4.2f %6.3f\n'%(time,x/proj.nmodes); time+=dt
            outf=jk(pr.simdir,'strip/evecs_msd2_average.dat')
            open(outf,'w').write(buf)
            return False
        elif task=='find MSD2-excess':
            """find MSD2-excess for each quasiharmonic mode from ptraj"""
            from amber.amber10 import evecs_projection
            from stats.statistics import linreg
            projf=jk(pr.simdir,'strip/projection.all.crd.gz')
            proj=evecs_projection(projf)
            msd2_dir=jk(pr.simdir,'strip/evecs_msd2')
            excessf=jk(msd2_dir,'msd2-excess.dat')
            buf='#i-hydr slope intercept corr-coeff  tmax msd2max excess\n'
            excesss=[]
            xl= 1.0 #1picosec
            xL=10.0 #10picosecs
            for i in range(proj.nmodes):
                infile=joink(msd2_dir,'%03d.dat'%(i+1,))
                t=[] ; t0=0.0 ; msd2=[] ; m0=0.0
                for line in open(infile).readlines()[2:]:
                    x,y=[float(z) for z in line.split()]
                    if x<xl:
                        if y>m0:
                            m0=y ; t0=x
                    elif x<xL:
                        t.append(float(x)); msd2.append(float(y))
                    else:
                        break #do not store values for t>xL
                #linear regression in the [1.0,10.0]psec time window
                regr=linreg(t,msd2); a=regr['slope']; b=regr['intercept']
                excess=m0-(b+a*t0); excesss.append(excess)
                buf+='%3d %f %f %f %f %f %f\n'%(i+1,a,b,regr['corrcoef'],
                                                t0,m0,excess)
            open(excessf,'w').write(buf)
            return False 
        elif task=='fit mode projections to Normal distribution':
            from amber.amber10 import evecs_projection
            projf=jk(pr.simdir,'strip/projection.all.crd.gz')
            proj=evecs_projection(projf)
            f=proj.fit2Normal()
            outdir=jk(pr.simdir,'strip/projection.all.hists')
            os.system('/bin/mkdir -p '+outdir)
            buf='#imode  mean      std       chi\n'
            for i in range( len(f['mean']) ):
                l='%4d %9.6f %9.6f %9.6f\n'%\
                   (i+1,f['mean'][i],f['std'][i],f['chi'][i])
                buf+=l
                buf2='#histogram of projections for mode %03d\n'%(i+1,)
                buf2+='#proj hist\n'
                for j in range( len(f['xvals'][i]) ):
                    buf2+='%8.3f %5d\n'%(f['xvals'][i][j],f['histogram'][i][j])
                open(jk(outdir,'%03d.dat'%(i+1,)),'w').write(buf2)
            outf=jk(pr.simdir,'strip/projection.all.fit2norm.dat')
            open(outf,'w').write(buf)
            return False
        elif task=='fract_boson_peak':

            """find fraction of modes important for the Boson peak in
            the structure factor.  We create structure factor plots by
            consecutively removing the contribution of one mode to the
            trajectory, then extract hydrogens into a trajectory, then
            compute the structure factor. We start with the mode with
            highest MSD2-excess, then continue with modes having
            decreasing MSD2-excess values."""
            from amber.amber10 import ptraj,evecs,top,crd
            from scatter.scatter import SIIScF,storeSIIScF,FourierSIIScf,store_sQW
            #store the MSD2-excess values, will serve for ranking of modes
            msd2f=jk(pr.simdir,'strip/evecs_msd2/msd2-excess.dat')
            msd2s={}
            for line in open(msd2f).readlines()[1:]:
                l=line.split() ; imode=int(l[0])
                msd2=float(l[-1]) ; msd2s[imode]=msd2
            #sort in drecreasing MSD2-excess values
            pairs=[ (msd2,imode) for (imode,msd2) in msd2s.items()]
            pairs.sort(); pairs.reverse()
            #instantiate topology object
            topf=jk(pr.simdir,'strip/ref15-shel.top')
            topobj=top(topf)
            nat=topobj.ITITL['NATOM']
            #instantiate evecs object and initialize some stuff
            evecsf=jk(pr.simdir,'strip/evecs.all.dat')
            evecsobj=evecs(evecsf)
            evecsobj.loadModes( nq=len(pairs) )  #don't load all modes
            evecsobj.massProcess( topobj.AMASS ) #pass atomic mass
            sM=evecsobj.sM
            #initialize scattering parameters
            qq=[2.0,] #Angstroms^(-1)
            dt=pr.times['twx']
            tl=40.0            #40picosec
            nt=int(tl/dt)      #number of time points
            tt=[]; ww=[]; dw=4.136/(nt*dt) #dw in meV when dt in psec
            for i in range(nt): tt.append(i*dt); ww.append(i*dw)
            #filter only Hhydrogens
            nH=self.Nhydr; fH=[]
            for iatom in range(1,nat+1):
                if iatom in self.unsolv['hydr']['atomnos']:
                    fH+=[1,1,1]
                else:
                    fH+=[0,0,0]
            fH=numpy.array(fH)

            print 'initialize trajectories'
            (a,prevcrdf)=mkstemp(dir='/projects/tmp')
            crdf='strip/ref15-shel.00001_50000.fit2centr.crd.gz'
            crdf=jk(pr.simdir,crdf)
            os.system('gunzip -c %s > %s'%(crdf,prevcrdf))
            print 'remove average conformation from trajectory'
            crdobj=crd(prevcrdf,nat)
            crdobj.loadframes()
            for i in range(crdobj.nframe):
                crdobj.frames[i]-=evecsobj.avgX
            #instantiate Hydrogen only frames
            Hframes=numpy.zeros(crdobj.nframe*nH*3).reshape(crdobj.nframe,nH,3)
            read=0
            outdir=jk(pr.simdir,'strip/evecs_fract_Boson_peak')
            os.system('/bin/mkdir -p '+outdir)
            #copy msd2-excess to outdir
            os.system('/bin/cp %s %s/'%(msd2f,outdir))
            #cycle through the modes
            for (msd2,imode) in pairs:
                read+=1
                print 'Processing mode %d %6.3f (%d)'%(imode,msd2,read)
                print '  remove mode projection from trajectory'
                Q=(evecsobj.qs[imode-1]).m.ravel() #quasiharmonic mode
                for i in range(crdobj.nframe):
                    fr=crdobj.frames[i].ravel()
                    q=(fr*sM*Q).sum() #mode projection
                    fr-=(sM*Q*q)      #remove projection
                    #store update frame
                    crdobj.frames[i]=fr.reshape(nat,3)
                    #keep only Hydrogens
                    Hframes[i]=(fr.compress((fr*fH)!=0.0)).reshape(nH,3)
                print '  intermediate factor scattering[nxyz][nq][nt]'
                scattering=SIIScF(qq,Hframes,dt,tl,nt)
                #average over all hydrogen atoms
                scattering=scattering.sum(axis=0)/len(scattering)
                #output to file
                outf=jk(outdir,'%03d.isf'%read)
                storeSIIScF(outf,scattering,qq,tt)
                print '  structure factor'
                sQW=FourierSIIScf(scattering,tt=tt,resolFun=convT_BASIS)
                #ouput to file
                outf=outf.replace('.isf','.sqw')
                store_sQW(outf,sQW,qq,ww)
            return False
        elif task=='ratio of areas':
            
            """Find how many modes we need to reduce the area under
            the Boson Peak by 50% with respect to the first area. We
            calculate the area as the integral [1meV,20mev] range. We
            fit the baselines previous to computing the ration of the
            areas. We calculate the baseline as the average of the
            structure factor in the [50mev,100mev] range."""
            
            buf= '#relative ratio area under Boson Peak\n'
            buf+='#subs ratio baseline\n'

            nmodes=200
            inpdir=jk(pr.simdir,'strip/evecs_fract_Boson_peak')
            areas=[]; baselines=[]
            for i in range(nmodes):
                area=0.0; nsampled=0; msampled=0; baseline=0.0
                inpf=jk(inpdir,'%03d.sqw'%(i+1,))
                for l in open(inpf).readlines()[1:]:
                    x,s=[ float(z) for z in l.split() ]
                    if x>=1.0 and x<=15.0:
                        nsampled+=1; area+=s
                    elif x>50.0 and x<100.0:
                        msampled+=1
                        baseline+=s
                baselines.append(baseline/msampled)
                areas.append(area*msampled/(nsampled*baseline))
            for i in range(nmodes):
                baseline=baselines[i]
                ratio=(areas[i]-areas[-1])/(areas[0]-areas[-1])
                buf+='%3d   %4.2f %6.4f\n'%(i,ratio,baseline)
            open(jk(inpdir,'ratio_of_areas.dat'),'w').write(buf)
            return False
        elif task=='excess areas':
            """Find the excess of S(Q,w) over a baseline given by a
            power-law. See Boukenter98 for justification of the
            baseline."""
            from math import exp,log  #natual logarithm
            buf= '#Excess area of S(Q,w) over a power-law baseline\n'
            buf+='#mode baseline S(Q,w) excess\n'
            nmodes=200;
            lnwm=0.00 #ln(w)_min=0.00 -> w_min= 1.00meV
            lnwM=2.91 #ln(w)_min=2.91 -> w_min=18.35meV
            inpdir=jk(pr.simdir,'strip/evecs_fract_Boson_peak')
            for i in range(nmodes):
                #find the two points for the power law
                inpf=jk(inpdir,'%03d.sqw'%(i+1,))                
                pw=0.0; ps=0.0; lnw1=None; lns1=None; lnw2=None; lns2=None
                for l in open(inpf).readlines()[2:]:
                    w,s=[ float(z) for z in l.split() ]
                    if log(w)>lnwm and not lnw1:#find first point
                        nw,ns=[ float(z) for z in l.split() ]
                        lnw1=( log(pw)+log(w)+log(nw) )/3.0
                        lns1=( log(ps)+log(s)+log(ns) )/3.0
                    if log(w)>lnwM and not lnw2:#find second point
                        nw,ns=[ float(z) for z in l.split() ]
                        lnw2=( log(pw)+log(w)+log(nw) )/3.0
                        lns2=( log(ps)+log(s)+log(ns) )/3.0
                    if lnw1 and lnw2: break
                    pw=w; ps=s
                area=0.0; blin=0.0; n=0; a=(lns2-lns1)/(lnw2-lnw1) #slope
                buf2='#  w   S_(w)\n'
                for l in open(inpf).readlines()[2:]:
                    w,s=[ float(z) for z in l.split() ]
                    if log(w)>lnwm:
                        area+=s
                        blin+=exp( lns1+a*(log(w)-lnw1) )
                        n+=1
                    if log(w)>lnwM: break
                buf+='%3d %f %f %f\n'%(1+i,blin/n,area/n,(area-blin)/n)
                open(jk(inpdir,'excess_areas.dat'),'w').write(buf)
            return False
        elif task=='cross temperatures mode projections':
            """find projection of the the three modes at T=360K having
            highest electrostatic fluctuatins (the three amigos)
            against the first 100 modes obtained for each other
            temperature"""
            from amber.amber10 import top,evecs
            from math import sqrt
            pr=self.solution['T_360']['Prod'+simLength]
            #instantiate topology object
            topf=jk(pr.simdir,'strip/ref15-shel.top')
            topobj=top(topf)
            #instantiate evecs object
            evecsf=jk(pr.simdir,'strip/evecs.all.dat')
            evecsobjA=evecs(evecsf)
            evecsobjA.loadModes(nq=7)
            #load mass matrix
            evecsobjA.massProcess( topobj.AMASS ) #pass atomic mass
            sM=evecsobjA.sM
            #cycle over the modes
            amigosID=[2,6,7] #mode indexes for the three amigos
            nmodes=100;
            buf= '#projections between three modes at T=360K and modes at other temperatures\n'
            buf2='modes with top similarity\n'
            for T in (100,140,180,220,260,300):
                avsim=0.0
                pps={} #store projections
                for id in amigosID: pps[id]=numpy.zeros(nmodes)
                args=tuple([T,]+amigosID)
                buf+='#\n#T=%3d\n#\n# i p(i,%1d) p(i,%1d) p(i,%1d)\n'%args
                pr=self.solution['T_'+str(T)]['Prod'+simLength]
                evecsf=jk(pr.simdir,'strip/evecs.all.dat')
                evecsobjB=evecs(evecsf)
                evecsobjB.loadModes(nq=nmodes)
                for Bmode in range(nmodes):
                    J=(evecsobjB.qs[Bmode]).m.ravel() #quasiharmonic mode
                    ps=[]
                    for Amode in amigosID:
                        Q=(evecsobjA.qs[Amode-1]).m.ravel()
                        pps[Amode][Bmode]=numpy.dot(sM*J,sM*Q)
                        ps.append(pps[Amode][Bmode]) #projection
                    ps.insert(0,1+Bmode)
                    buf+='%3d %6.3f %6.3f %6.3f\n'%tuple(ps)
                #sort projections, output top three for every mode
                buf2+='T=%3d\n'%T
                for Amode in amigosID:
                    buf2+='mode%d'%Amode
                    topidx=numpy.abs(pps[Amode]).argsort()[::-1][0:3]
                    for idx in topidx: #top three
                        avsim+=abs(pps[Amode][idx])
                        buf2+=' %3d(%6.3f)'%(1+idx,pps[Amode][idx])
                    buf2+='\n'
                buf2+='average similarity: %f\n'%(avsim/9.,)
            print buf2
            pr=self.solution['T_360']['Prod'+simLength]
            outf=jk(pr.simdir,'strip/crossTmodeProj.dat')
            open(outf,'w').write(buf)
            return False
        elif task=='Z-score ILE40 ILE7 methyl hydrogens':
            """How to assess the contribution of ALA40 and ALA7 methyl
hydrogens to the excess-MSD2 of mode 7 at T=100 K:

since q(t)= E^t * M^0.5 * x(t), we look at matrix elements e_{7,j}
corresponding to eigenvector of mode 7. In particular, we look at
e_{7,j} with j corresponding to degrees of freedom of ALA40 and ALA7
methyl hydrogens. Let:

   S_1=\sum_{j} K_j e_{7,j}^2 / \sum_{j} K_j, where K_j=1 if j is one
   of the degrees of freedom of ALA40 and ALA7 methyl hydrogens, and
   zero otherwise.

   Let S_0 and S_2 be the average and standard deviation of the
   distribution e_{7,j}^2 values where j correspond to one degree of
   freedom of any hydrogen in rubredoxin, respectively.

We can define a Z-score for the set of e_{7,j}^2 elements with j
corresponding to degrees of freedom of ALA40 and ALA7 methyl
hydrogens:

  Z=(S_1-S_0)/S_2

If Z>>1, then we will know that ALA40 and ALA7 methyl hydrogens will
be overrepresented in mode 7 when compare to the \"average
hydrogen\". Then we can stablish that the hihg excess-MSD^2 of mode 7
is correlated with the high excess-MSD^2 of the ALA40 and ALA7 methyl
hydrogens at this temperature (T=100 K)."""
            from amber.amber10 import ptraj,evecs,top,crd
            Hns=[121,122,123,591,592,593] #at num ILE7,ILE40 methyl hydr
            #instantiate topology object
            topf=jk(pr.simdir,'strip/ref15-shel.top')
            topobj=top(topf)
            nat=topobj.ITITL['NATOM']
            #instantiate evecs object and initialize some stuff
            evecsf=jk(pr.simdir,'strip/evecs.all.dat')
            evecsobj=evecs(evecsf)
            evecsobj.loadModes( nq=7 )  #load first seven modes
            eigenmode=evecsobj.qs[6].m #mode 7
            q=[] #displacements of methyl hydrogens in mode 7
            for index in Hns: q.append(eigenmode[index-1])
            q=numpy.array(q)
            Q=[] #displacements of hydrogens in mode 7
            for index in self.unsolv['hydr']['atomnos']:
                Q.append(eigenmode[index-1])
            Q=numpy.array(Q)
            #Z-score of displacements
            Q=(Q*Q).sum(axis=1); q=(q*q).sum(axis=1)
            Z_score=( q.mean()-Q.mean() ) / Q.var()
            print 'Z score=',str(Z_score) #Z_score=7131
            for i in range(self.Nhydr):
                k=self.unsolv['hydr']['atomnos'][i]
                print '%3d %f'%(k,Q[i])

            return False
        """obtain quasiharmonic modes for CA, heavy-atoms,
        hydrogen-atoms, and all-atoms"""
        script=junkName() ; traj=junkName()
        os.system('gunzip -c %s > %s'%\
                  (joink(pr.simdir,pr.strip['fit2centr']),traj) )
        masks={'CA':'@CA','heavyat':'!@H*','H':'@H*','all':''}
        for (maskk,maskv) in masks.items():
            print 'analyzing quasiharmonic modes for '+maskk+' atoms'
            junkf=junkName()
            open(script,'w').write("""trajin %s 1 50000 5
matrix mwcovar name mwcovar %s out mwcovar.dat byatom
analyze matrix mwcovar out evecs.dat thermo vecs %d
analyez matrix mwcovar out %s thermo vecs 0"""%(traj,maskv,nvecs,junkf))
            cmd='ptraj10 %s < %s > eigen.dat '%\
                 (self.unsolv['parmfile'],script)
            os.system(cmd+' && /bin/rm  '+junkf)

            mw='strip/mwcovar.'+maskk+'.dat'
            eg='strip/evecs.'+maskk+'.dat'
            ev='strip/eigen.'+maskk+'.dat'
            os.system( 'mv mwcovar.dat '+joink(pr.simdir,mw) )
            os.system( 'mv evecs.dat '+joink(pr.simdir,eg)  )
            os.system( 'mv eigen.dat '+joink(pr.simdir,ev)  )
            pr.strip['fit2centr.mwcovar.'+maskk]=mw
            pr.strip['fit2centr.evecs'+maskk]=eg
            pr.strip['fit2centr.eigen'+maskk]=ev
            pr.dump() #store these dictionary items
        os.system('/bin/rm %s %s '%(traj,script))
        return True
    

    def quasiHforH(self,T='300',nvecs=1000,simLength='',
                   task='init onlyH',filter=None,outf=None):
        """obtain quasiharmonic modes with my own library"""
        from amber.Prod import Prod
        from quasiH.quasi import modes
        from amber.amber10 import crd
        from utilities.readingWritingFiles import write_from_numpy
        if 'onlyH' in task: #work for only hydrogens trajectory
            pH=self.solution['T_'+T]['onlyH'+simLength]
            nat=pH['Nhydr']; nd=nat*3; qH=modes(nd); qH_p=None
            sd=jk(pH['simdir'],'qH');
            if 'init' in task:
                 os.system('mkdir -p '+sd)
                 if 'qH' in pH.keys():
                     pwrite('already initialized. continue?(y/n): ')
                     yn=sys.stdin.readline().strip()
                     if yn[0].lower()=='n':
                         return False
                 pH['qH']={'sd':sd}; qH_p=pH['qH']
            else:
                qH_p=pH['qH']
                qH.pickleLoad( open(jk(sd,qH_p['obj'])) )
                if 'projectTrajectory' in task:
                    """find mode components of each trajectory snapshot
                    """
                    crdf=pH['crd.fit2prev']; crdobj=crd(crdf,nat)
                    qH.loadM(numpy.ones(nd))
                    qH.findVarCovar(crdobj.nextravel) #var-cov matrix
                    crdobj.rewind() #rewind to first frame
                    pout=open(jk(qH_p['sd'],'crd'),'w')
                    pout.write('#trajectory of quasiharmonic modes\n')
                    qH.projectTrajectory(crdobj.nextravel,pout)
                    qH_p['evals']=jk(qH_p['sd'],'evals.dat')
                    buf='\n'.join([str(x) for x in qH.evals])
                    open(qH_p['evals'],'w').write(buf)
                    os.system('gzip -f '+jk(sd,'crd'))
                    qH_p['crd']=jk(sd,'crd.gz')#compress trajectory
                elif 'extractTrajetory' in task:
                    #note that nat=Nhydr, even though qH_p['crd'] is a
                    #trajectory of modes (there are 3*nat of them). It
                    #is OK as long as we do frame.ravel() after
                    #reading each frame
                    maxSize=50E+6; crdobj=crd(qH_p['crd'],nat)
                    sys.stdout.write('enter mode number (1 for first): ')
                    J=int(sys.stdin.readline().strip()); I=J-1
                    sys.stdout.write('enter file name: ')
                    outf=sys.stdin.readline().strip(); pt=open(outf,'w')
                    buf='#trajectory of mode '+str(J)+'\n'; iframe=0;
                    while crdobj.loadframe():
                        f=crdobj.frame.ravel();iframe+=1
                        buf+=str(f[I])+'\n'; #print 'iframe=',iframe
                        if len(buf)>maxSize:
                            pt.write(buf); buf=''
                    pt.write(buf); pt.close()
                    return False #no need to update the PEPSI object
                elif 'filterTrajectory' in task:

                    """filter trajectory of hydrogens by keeping only
                    the components of selected modes. We pass array
                    'filter' with filter.shape=(3*N_hydrg,) array of
                    zeros and ones, and will keep only those modes for
                    which the entry in the array is one"""
                    
                    #note below nat=Nhydr, even though qH_p['crd'] is a
                    #trajectory of modes (there are 3*nat of them). It
                    #is OK as long as we do frame.ravel() after
                    #reading each frame
                    crdobj=crd(qH_p['crd'],nat) #traj in mode components
                    #enter trajectory file name
                    if not outf:
                        sys.stdout.write('enter file name: ')
                        outf=sys.stdin.readline().strip();
                    pt=open(outf,'w')
                    #filter trajectory
                    iframe=0;
                    while crdobj.loadframe():
                        q=crdobj.frame.ravel();iframe+=1
                        x=qH.avF.ravel()+qH.mode2cart(q*filter) #cart. XYZ
                        write_from_numpy(pt,x,ncols=7)
                    print str(iframe),' frames recorded'
                    os.system('gzip -f '+outf) #save space!
                    return False #we did not modified the qH object
                elif 'fix chiSquare' in task:
                    """chiSquare was not well calculated in the next
                    'elif'. We fix that here"""
                    from stats.histograms import chiSquare
                    buf='#Chi-square values of mode vale distributions\n'
                    for i in range(1,nd):
                        file=jk(qH_p['mode_hists'],'%04d.dat'%i)
                        X=[]; Y=[]
                        for l in open(file).readlines():
                            if l[0]=='#': continue
                            (x,y)=[float(z) for z in l.split()]
                            X.append(x); Y.append(y)
                        X=numpy.array(X); Y=numpy.array(Y)
                        chi=chiSquare( X, Y, lowcut=0.05 )
                        buf+='%04d %f\n'%(i,chi)
                    open(qH_p['chisq'],'w').write(buf)
                elif 'chiSquare' in task:
                    """calculate fit to of mode value distribution to
                    normal distribution"""
                    #note that nat=Nhydr, even though qH_p['crd'] is a
                    #trajectory of modes (there are 3*nat of them). It
                    #is OK as long as we do frame.ravel() after
                    #reading each frame
                    crdobj=crd(qH_p['crd'],nat);
                    chi2,hist,bins=qH.chisquare(crdobj.next,crdobj.rewind)
                    #store chi2
                    qH_p['chisq']=jk(sd,'chisquare.dat')
                    buf='#Chi-square values of mode vale distributions\n'+\
                         '\n'.join([str(x)for x in chi2])
                    open(qH_p['chisq'],'w').write(buf)
                    #store histograms
                    nbin=len(bins[0])
                    qH_p['mode_hists']=jk(sd,'mode_hists')
                    os.system('/bin/mkdir -p '+qH_p['mode_hists'])
                    for i in range(3*nat):
                        h=hist[i]; b=bins[i]
                        outf=jk(sd,'mode_hists','%04d.dat'%(i+1))
                        buf='#distribution of mode values for %04d\n'%(i+1)
                        for j in range(nbin): buf+='%f %f\n'%(b[j],h[j])
                        open(outf,'w').write(buf)
                elif 'displacements as bfactors' in task:
                    """output pdb files, one per mode. Bfactors are
                    displacements of the atoms in the particular
                    mode"""
                    pdbf=self.unsolv['preminPDB']#template pdb file
                    dir=jk(sd,'color_displ')     #place to store pdb files
                    os.system('/bin/mkdir -p '+dir)
                    indexes=[ x-1 for x in self.unsolv['hydr']['atomnos'] ]
                    qH.colorDispl(pdbf,dir,indexes=indexes)
                    qH_p['color_displ']=dir
        #save updated gn object to file 
        qH_p['obj']=jk(sd,'qH.dump');
        qH.pickleDump( open(qH_p['obj'],'w') )
        return True

    def getNvecs(self,T='300',simLength=''):
        """find number of quasiharmonic modes computed"""
        from amber.amber10 import evecs
        Prod=self.solution['T_'+T]['Prod'+simLength]
        for maskk in ['CA','heavyat','H','all']:
            evecsf=joink(Prod.simdir,Prod.strip['fit2centr.evecs'+maskk])
            ev=evecs(evecsf)
            Prod.strip['fit2centr.evecs'+maskk+'.nvecs']=ev.find_nq()
        Prod.dump()
        return True
    
    def quasiharmonic_movies(self,T='300',simLength='',nmodes=200):
        """create PDB movies for each quasiharmonic mode"""
        from amber.amber10 import evecs
        nmodes=int(nmodes)
        
        Prod=self.solution['T_'+T]['Prod'+simLength]
        topfile=self.unsolv['parmfile']
        key='fit2centr.evecsall'; outdir='strip/'+key+'d'
        ev=evecs( joink(Prod.simdir,Prod.strip[key]) )
        outf=joink(Prod.simdir,outdir)
        ev.create_movies( topfile,outdir=outf,nmodes=nmodes )
        #Prod.strip['fit2centr.evecsalld']=outdir
        #Prod.dump()
        return False
 
    def generate_movies(self,T='300',simLength=''):
        """create PDB movies for each quasiharmonic mode with gen_movies"""
        from amber.amber10 import evecs
        
        Prod=self.solution['T_'+T]['Prod'+simLength]
        templatepdb=joink(self.simdir,'ref15-shel.uns.premin.pdb')
        outdir=joink(Prod.simdir,'strip/gen_movies')
        ev=evecs( joink(Prod.simdir,Prod.strip['fit2centr.evecsall']) )
        ev.gen_movies(templatepdb,outdir=outdir)
 
        return False
    
    def projectCRDontoEVECS(self,T='300',simLength=''):
        """For every snapshot, find the projection of each trajectory
        snapshot onto the quasi-harmonic modes"""
        from amber.amber10 import crd,top,evecs
        Prod=self.solution['T_'+T]['Prod'+simLength]
        #instantiate object topology
        topf=joink(Prod.simdir,'strip/ref15-shel.top')
        topobj=top(topf)
        #instantiate object trajectory
        crdf=joink(Prod.simdir,Prod.strip['fit2centr'])
        crdobj=crd(crdf,topobj.ITITL['NATOM'])
        #instantiate object evecs
        evecf=joink(Prod.simdir,Prod.strip['fit2centr.evecsall'])
        evecobj=evecs(evecf)
        evecobj.massProcess(topobj.AMASS) #retrieve masses
        #project trajectory onto evecs
        outf=evecf.replace('.all.dat','.all.project.dat')#store projections
        evecobj.loadModes(nq=200)
        evecobj.projectCRD(crdobj,outf)
        #gzip and store to CAN object
        os.system('gzip -f '+outf)
        outf='strip/'+basename(outf)+'.gz' #update outf variable
        Prod.strip['fit2centr.evecsall.project']=outf
        Prod.dump() #record new (key,value) added to Prod object
        return True

    def corr_evecs_rmsd(self,T='300',simLength=''):
        """Load the projections of the quasiharmonic modes onto the
        trajectory, then find correlation of these projections to rmsd
        to reference structure values"""
        from utilities.readingWritingFiles import read_column
        from amber.amber10 import qhm_proj
        from stats.statistics import corrcoeff
        Prod=self.solution['T_'+T]['Prod'+simLength]
        rmsf=joink(Prod.simdir,Prod.strip['rmsRef'])
        rmss=numpy.array([float(x) for x in read_column(rmsf,2)])
        #correlations between RMSD and absolute value of the projections
        outf=joink(Prod.simdir,Prod.strip['fit2centr.evecsall.project'])
        qhmp=qhm_proj(outf) #initialize projections object
        corrs=[]
        for m in qhmp.loadProjections():            
            corrs.append(corrcoeff(numpy.abs(m),rmss))
        #store correlations
        buf='#correlation coefficient between RMSD values and\n'+\
             '#projections of the frames onto the quasi-harmonic modes\n'+\
             '#\n'+\
             '#mode-number correl-coeff\n'
        for imode in range(len(corrs)):
            buf+='%04d %6.3f\n'%(imode+1,corrs[imode]); imode+=1
        key='fit2centr.evecsProj_rms.corr'
        Prod.strip[key]='strip/'+key+'.dat'
        outf=joink(Prod.simdir,Prod.strip[key])
        open(outf,'w').write(buf)
        Prod.dump()
        return True

    def contactSurvivalTimes(self,T='300'):
        """given a set of native contacts, calculate set of survival
        times for a contact map trajectory"""
        import numpy
        #read native contacts, put them in a list
        cmnat=[]
        pcmnat=open(self.natmap,'r') ; Nnat=int(pcmnat.readline().strip())
        for i in range(Nnat):
            c=[int(x) for x in pcmnat.readline().split()]
            c.sort()
            cmnat.append(c)
        #read contact map trayectoy. Create a
        #11111101110101100100010000 string for every native contact,
        #representing whether it is present of absent in each frame of
        #the trayectory
        cm={}
        for c in cmnat: cm[c]=''
        pr=self.solution['T_'+T]['Prod']
        pcm=open(pr.strip['cm'],'r') ; l=pcm.readline()
        frame=0
        while l:
            frame+=1 ; Nc=int(l.strip())
            for ic in range(Nc):
                c=[int(x) for x in pcmnat.readline().split()]
                c.sort()
                if c in cmnat: 
                    cm[c]+='1'
                else:
                    cm[c]+='0'
            l=pcm.readline()
        #calculate survival time for every native contact
        buf=''
        for c in cmnat:
            tr=numpy.array( [int(x) for x in list(cm[c])] )
            emin=tr.sum() ; e=emin ; tmins=[1,]
            for t in range(len(tr)):
                e=e+1-2*tr[t]
                if e==emin:
                    tmins.append(t)
                elif e<emin:
                    emin=e ; tmins=[t,]
            tmin=float(tmins.sum())/len(tmins)
            buf+='%3d %3d %f\n'%(c[0],c[1],tmin)
        pr.strip['survT']='strip/survivalTimesMap.dat'
        open(joink(pr.simdir,pr.strip['survT'])).write(buf)
        pr.dump()
        return True
    
    def genInit(self,T='300'):
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
ntwe=1000, ntpr=1000, ntwx=20, ntwv=20
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
        buf=re.compile( 'nstlim=\d+').sub( 'nstlim=5000', buf)
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
        #pr.doStrip(self.unsolv['parmfile'],self.unsolv['crdfile'],
        #           maskl=[':WAT',':Na+']) #remove water and sodium ions
        if 'avg' not in pr.prev.keys(): pr.getAvg()#averages
        pr.rmsJump(outcrd=True)#consec rms and output fitted trajectory
        pr.rmsRef()        #rms with respect to reference
        pr.centr_and_avg() #find centroid
        pr.rmsCentr()      #rms with respect to centroid
        pr.rmsFirst()      #rms with respect to first frame
        pr.fit2centr()     #fit trajectory to centroid
        self.solution[TT]['Prod'+simLength]=pr #store
        return True

    def getZMatrix(self,T='300',simLength='',task='init all'):
        """calculate trajectory of Z-Matrix"""
        from amber.amber10 import crd
        from tertiary.mopin import mopin
        Prod=self.solution['T_'+T]['Prod'+simLength] #Prod object
        nat=Prod.strip['nat']
        MP=mopin()
        if task=='init all':
            #initialize the Z-Matrix 'factory'
            MP.loadPDB(self.unsolv['preminPDB'])
            crdobj=crd(jk(Prod.simdir,Prod.strip['fit2prev']),nat)
            outf='strip/zmatrix.crd'; absoutf=jk(Prod.simdir,outf)
            MP.zmTraj(crdobj.next,open(absoutf,'w'))#create Z-Matrix traj
            #store trajectory
            Prod.strip['zmatrix.crd']=outf+'.gz'
            os.system('gzip -f '+absoutf) #compress trajectory
            Prod.dump()
            #store mopin object into file
            self.mopinf=jk(self.simdir,'mopin.dump')
            #store connectivity map
            MP.topmap() #find connectivity map
            MP.pickleDump(open(self.mopinf,'w'))
        if task=='dihedrals only':
            outf='strip/dihedrals.crd'; absoutf=jk(Prod.simdir,outf)
            MP.pickleLoad(open(self.mopinf))
            crdobj=crd(jk(Prod.simdir,Prod.strip['fit2prev']),nat)
            hists=MP.dihedrals(crdobj.next,open(absoutf,'w'))
            outd=jk(Prod.simdir,Prod.strip['dihedrals_histograms'])
            os.system('/bin/mkdir -p '+outd)
            for i in range(nat):
                cmt='#histogr of dihedral angles for atom C-index %04d\n'%i
                buf=cmt+'\n'.join([ str(x) for x in hists[i] ])
                open('%s/%04d.dat'%(outd,i),'w').write(buf)
            MP.pickleDump(open(self.mopinf,'w'))
        if task=='torsions only':
            outf='strip/torsions.crd'; absoutf=jk(Prod.simdir,outf)
            MP.pickleLoad(open(self.mopinf))
            MP.topmap() #find connectivity map
            crdobj=crd(jk(Prod.simdir,Prod.strip['fit2prev']),nat)
            hists=MP.torsions(crdobj.next,open(absoutf,'w'))
            outd=jk(Prod.simdir,Prod.strip['torsion_histograms'])
            os.system('/bin/mkdir -p '+outd)
            for i in range(nat):
                cmt='#histogr of torsion angles for atom C-index %04d\n'%i
                buf=cmt+'\n'.join([ str(x) for x in hists[i] ])
                open('%s/%04d.dat'%(outd,i),'w').write(buf)
            MP.pickleDump(open(self.mopinf,'w'))

        return True

    def dihedralHist(self,T='300',simLength=''):
        """store the dihedral distribution for each atom"""
        from utilities.readingWritingFiles import write_from_numpy
        Prod=self.solution['T_'+T]['Prod'+simLength]
        nat=Prod.strip['nat']
        hists=numpy.zeros(nat*360).reshape(nat,360)

        Prod.strip['zmatrix.crd']='strip/zmatrix.crd.gz'
        
        crdf=jk(Prod.simdir,Prod.strip['zmatrix.crd'])
        tmpf=mkstemp()[-1]
        os.system('gunzip -c %s > %s'%(crdf,tmpf));
        #go frame by frame
        pt=open(tmpf); l=pt.readline(); iframe=0
        while l:
            frame=[float(x) for x in l.split()]; iframe+=1
            frame=numpy.array(frame).reshape(nat,3)
            dihedrals=(frame.transpose())[-1]
            for i in range(nat):
                angle=int(dihedrals[i])
                hists[i][angle]+=1
            l=pt.readline()
            if not iframe%1000: print iframe
            #only for debugging
            #print 'iframe=',iframe;
            #if iframe==10000: break
            #end of only for debugging
        pt.close()
        #store histograms
        Prod.strip['dihedrals_histograms']='strip/dihedrals_histograms'
        outd=jk(Prod.simdir,Prod.strip['dihedrals_histograms'])
        os.system('/bin/mkdir -p '+outd)
        for i in range(nat):
            cmt='#histogram of dihedral angles for atom C-index %04d\n'%i
            buf=cmt+'\n'.join([ str(x) for x in hists[i] ])
            open('%s/%04d.dat'%(outd,i),'w').write(buf)
        Prod.dump()
        os.system('/bin/rm '+tmpf)
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
        if self.isDeutd: mass+=self.Nhydr
        self.mass=mass; print 'mass=%fa.m.u'%(mass,)
        return True
        
    def needed_waters(self,r='0.18'):
        """print number of needed waters molecules for
        a given mass(H20)-to-mass(protein) ratio"""
        x=float(r)
        watermass=18.0 ; #Nwat=int(self.mass/watermass)
        if self.isDeutd: watermass=22
        Nrwat=int(x*self.mass/watermass)
        if not 'hydrated' in self.__dict__.keys(): self.hydrated={}
        self.hydrated[r]={}
        self.hydrated[r]['nwat']=Nrwat; print 'n-waters=%d'%(Nrwat,)
        return True

    def prepare_raw_pdb(self,r='0.18'):
        """Extract as many waters as needed from the solvated premin
        in order to set the hydration ratio to preselected value"""

        buf="""trajin _PREMINCRD_
strip :_STRIPMASK_
solvent byres :WAT
closestwater _NWAT_ :1-_L_
trajout _OUTPDB_ pdb"""
        pr=self.hydrated[r]
        if not 'simdir' in pr.keys():
            simdir=joink(self.simdir,'hydrated',r)
            os.system('/bin/mkdir -p '+simdir)
            pr['simdir']=simdir
        pr['rawpdb']=joink(pr['simdir'],self.id+'.raw.pdb')

        for key,val in {'_PREMINCRD_':self.crdfile,
                        '_STRIPMASK_':'Na+',
                        '_NWAT_'     :str(pr['nwat']),
                        '_L_'        :str(self.L),
                        '_OUTPDB_'   :pr['rawpdb']}.items():
            print key,val
            buf=buf.replace(key,val)
        infile=junkName(); open(infile,'w').write(buf)
        #run ptraj
        cmd='ptraj9 %s < %s && '%(self.parmfile,infile)+\
             '/bin/rm %s    && '%(infile)+\
             '/bin/mv %s %s    '%(pr['rawpdb']+'.1',pr['rawpdb'])
        os.system(cmd)
        return True

    def extractHfromPDB(self,task='init'):
        """extract hydrogen lines from preminPDB"""
        pdbf=self.pdbfile
        if task=='premin': pdbf=self.preminPDB
        buf=''
        for line in open(pdbf,'r').readlines():
            if line[17:20]=='WAT': break
            if line[0:5]=='ATOM ' and line[13]=='H': buf+=line
        if task=='init':
            self.onlyHpdb=self.pdbfile.replace('.ent','.onlyH.ent')
            open(self.onlyHpdb,'w').write(buf)
        elif task=='premin':
            self.onlyHpreminPDB=self.preminPDB.replace('.pdb','.onlyH.pdb')
            open(self.onlyHpreminPDB,'w').write(buf)
        elif task=='uns.premin':
            self.unsolv['onlyHpreminPDB']=self.unsolv['preminPDB'].replace('.pdb','.onlyH.pdb')
            open(self.unsolv['onlyHpreminPDB'],'w').write(buf)
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

    def hydrWork(self,task='init'):
        """create dictionary 'hydr' under self.unsolv and add properties"""
        if task=='init':
            self.unsolv['hydr']={}
        elif task=='list index':
            """create list of atom numbers for the hydrogens with PDB
            convention"""
            ls=[]
            for l in open(self.unsolv['preminPDB']).readlines():
                if l[0:5]=='ATOM ' and l[13]=='H': ls.append(l[6:11])
            self.unsolv['hydr']['atomnos']=numpy.array(ls,dtype=int)
        return True

    def Ktraj(self,T='270',simLength=''):
        """create trajectory of kinetic energy"""
        from amber.amber10 import crd
        from pdbfiles.smallThings import insertBfact
        from utilities.readingWritingFiles import write_from_numpy
        Prod=self.solution['T_'+T]['Prod'+simLength] #Production object
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

    def keepH(self,T='300',simLength='',recover=False):
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
            buf='different number of hydrogens than stored!\n%d vs %d'%\
                 (self.Nhydr,Nhydr)
            Bye(buf)
        else:
            self.Nhydr=Nhydr
        #pdb.set_trace()
        Prod=self.solution['T_'+T]['Prod'+simLength] #Production object
        if 'crdfile' in pH.keys():
            pH['crd.fit2centr']=pH['crdfile']
            del pH['crdfile']
        #for keyword in ('crdfile','fit2prev','fit2centr'):
        for keyword in ('fit2prev','fit2centr'):
            print 'processing '+keyword
            crdf=joink(Prod.simdir,Prod.strip[keyword])
            outf=basename(crdf.replace('.crd.gz','.H.crd'))
            outf=joink(pH['simdir'],outf)
            ptr=ptraj(self.unsolv['parmfile'],crdf)
            ptr.strip(entities=['!@H*',]) #keep only hydrogens
            if not recover:
                ptr.go(outcrd=outf,outfile='./junk.log')
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
        crdf=joink(Prod.simdir,Prod.strip['velfile'])
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
        #pr.BfactorsByMask(mask='@H*')
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

    def av_msd2(self,T='300',simLength='',task=''):
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
            outf=joink(onlyH_p['simdir'],'H.av.msd2.dat')
            onlyH_p['H.av.msd2']=outf
            open(outf,'w').write(buf)
            return True
        elif task=='onlyH gnmodes':
            g_p=self.solution['T_'+T]['onlyH'+simLength]['gn']
            crdf=g_p['crd'];
            crdobj=crd(crdf,self.Nhydr)
            msd2=crdobj.msd2()
            print '\n'.join([str(x) for x in msd2])
        
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


    def createMMTKtraj(self,T='300',recover=''):
        """create an MMTK trajectory from the consecutive trajectory"""
        import pdb
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

        prmcrd=junkName() ; prmvel=junkName() #temporary files
        prmcrd0=joink(Prod.simdir,Prod.strip['fit2prev'])
        print 'unzipping %s to %s'%(basename(prmcrd0),prmcrd)
        cmd='gunzip -c %s > %s'%(prmcrd0,prmcrd)
        os.system(cmd)
        prmvel0=joink(Prod.simdir,Prod.strip['velfile'])
        print 'unzipping %s to %s'%(basename(prmvel0),prmcrd)
        cmd='gunzip -c %s > %s'%(prmvel0,prmvel)
        os.system(cmd)
        pM['pdbfixed']=joink(pM['simdir'],
                             basename(prmpdb).replace('.pdb','.fixed.pdb'))
        pM['pdb']=joink(pM['simdir'],
                        basename(prmpdb).replace('.pdb','.mmkt.pdb'))
        print 'converting to MMTK trajectory...'
        convertAmberToMMTK(prmpdb,prmcrd,pM['nc'],prmvel=prmvel,
                           format='crd',dtime=Prod.times['twx'],
                           pdbfixed=pM['pdbfixed'],pdbtest=pM['pdb'])
        os.system('/bin/rm %s %s'%(prmcrd,prmvel))
        print '...converting finished'
        return True

    def averageISF(self,selection=None,out_isf=None,out_dsf=None,
                   T='300',keyword=None):
        from Scientific.IO.NetCDF import NetCDFFile
        from Numeric import zeros,Float
        from my_nMoldyn.misc import cloneNetCDF
        pM=self.solution['T_'+T]['nMoldyn']
        if not selection: selection='*'
        if not keyword: keyword='av'
        #average isf
        gooble=joink(pM['simdir'],'isf',selection+'.nc')
        files=chomp(os.popen('ls -1 '+gooble).readlines())

        isf=zeros( NetCDFFile(files[0],'r').variables['sf'].shape,Float)
        for file in files:
            isf+=NetCDFFile(file,'r').variables['sf'][:][:]
        isf=isf/len(files)
        #  save results
        if not out_isf: out_isf=joink(pM['simdir'],keyword+'_isf.nc')
        ncfout=cloneNetCDF(out_isf,files[0])
        ncfout.variables['sf'][:]=isf
        ncfout.close()       
        pM[keyword+'_isf']=out_isf
        #average dsf
        gooble=joink(pM['simdir'],'isf_spectr',selection+'.nc')
        files=chomp(os.popen('ls -1 '+gooble).readlines())
        dsf=zeros(NetCDFFile(files[0],'r').variables['dsf'].shape,Float)
        for file in files:
            dsf+=NetCDFFile(file,'r').variables['dsf'][:][:]
        dsf=dsf/len(files)
        #  save results
        if not out_dsf: out_dsf=joink(pM['simdir'],keyword+'_dsf.nc')
        ncfout=cloneNetCDF(out_dsf,files[0])
        ncfout.variables['dsf'][:]=dsf
        ncfout.close()       
        pM[keyword+'_dsf']=out_dsf
        return False
        
    def printISF(self,T='300',infile='/projects/rubredoxin/out/1fhm_D/T_180/nMoldyn/av_isf.nc',outfile='/projects/rubredoxin/out/1fhm_D/T_180/nMoldyn/junk.png',outtxt='/projects/rubredoxin/out/1fhm_D/T_180/nMoldyn/junk.txt'):
        from scipy.io.netcdf import netcdf_file
        from pylab import hold,plot,savefig,imshow,contour
        ncfile=netcdf_file(infile,'r')
        q=ncfile.variables['q'][:] ; L=len(q)
        t=ncfile.variables['time'][:]
        z=ncfile.variables['sf']
        hold(True)
        all=['']*len(t)
        for i in range(len(t)): all[i]+='%f '%(t[i])
        for iq in range(len(q)):
            plot(t[0:],z[iq][0:])
            for i in range(len(t)): all[i]+='%f '%z[iq][i]
        open(outtxt,'w').write('\n'.join(all))
        #pdb.set_trace()
        #matrixplot3d(ncfile.variables['sf'][:][:])
        #imshow(z,origin='lower', extent=[q[0],q[-1],t[0],t[-1]])
        #contour(z,origin='lower', extent=[q[0],q[-1],t[0],t[-1]])
        savefig(outfile)
        return False
    
    def printDSF(self,T='300',infile='/projects/rubredoxin/out/1fhm_D/T_180/nMoldyn/av_dsf.nc',outfile='/projects/rubredoxin/out/1fhm_D/T_180/nMoldyn/junk.png',outtxt='/projects/rubredoxin/out/1fhm_D/T_180/nMoldyn/junk.txt'):
        from scipy.io.netcdf import netcdf_file
        from pylab import hold,plot,savefig,imshow,contour
        ncfile=netcdf_file(infile,'r')
        q=ncfile.variables['q'][:] ; L=len(q)
        t=ncfile.variables['frequency'][:]
        z=ncfile.variables['dsf']
        hold(True)
        all=['']*len(t)
        for i in range(len(t)): all[i]+='%f '%(t[i])
        for iq in range(len(q)):
            plot(t[0:50],z[iq][0:50])
            for i in range(len(t)): all[i]+='%f '%z[iq][i]
        open(outtxt,'w').write('\n'.join(all))
            
        #pdb.set_trace()
        #matrixplot3d(ncfile.variables['sf'][:][:])
        #imshow(z,origin='lower', extent=[q[0],q[-1],t[0],t[-1]])
        #contour(z,origin='lower', extent=[q[0],q[-1],t[0],t[-1]])
        savefig(outfile)
        return False

    def findWeights(self,deut_sel=None,outfile=None):
        """find incoherent scattering weights
        deut_sel: deuterium selection"""
        from nMoldyn.misc import inputFileRead
        from tempfile import gettempdir, mktemp
        template="""from MMTK import * 
title = 'Incoherent Scattering Function' 
trajectory = ['*TRAJ*',] 
log_file = 'isf_*BTRAJ*' 
time_info = (0, 49999, 1) 
atoms = {'Protein.0': ['*'], 'FE3': ['*']} 
deuter = *DEUTER*
q_vector_set = ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0], 0.100000, 50, None) 
weights = 'incoherent' 
ft_window = 50 
units_q = 1/Units.nm 
time_steps = 2000 
frequency_points = 2000 
output_files = {'fft': 'ISF_SPECT_*ISFOUT*', 'isf': 'ISF_*ISFOUT*'}"""
        pM=self.solution['T_100']['nMoldyn']
        query=template.replace( '*TRAJ*',pM['nc'] )
        query=query.replace( '*BTRAJ*',basename(pM['nc']) )
        query=query.replace( '*ISFOUT*',basename(pM['nc']) )
        if self.isDeutd:
            query=query.replace("*DEUTER*","{'Protein.0': ['Hydrogen']}")
        else:
            query=query.replace("*DEUTER*","None")
        inputf=mktemp() ; open(inputf,'w').write(query)
        data = inputFileRead(inputf)
        os.system('/bin/rm '+inputf)
        w={}
        for atom in data.atoms.atomList():
            w[atom.index]='%12.10f %s'%(data.weights[atom],atom.name)
        indexes=w.keys()  ; indexes.sort()
        buf=''
        for index in indexes:
            buf+='%04d %s\n'%(index,w[index])
        outfile=joink(self.simdir,'incoherent_weights.dat')
        open(outfile,'w').write(buf)
        self.incoh_weights=outfile ;print outfile
        return True

    def checknMoldynISF(self,T='300'):
        """calculate intermediate and final structure factors"""
        from copy import deepcopy
        from tempfile import gettempdir, mktemp
        from MMTK import Collection
        from nMoldyn.misc import inputFileRead
        from my_nMoldyn.core import IncoherentScatteringFunction
        from nMoldyn.core import ScatteringFunctionFFT
        isf=IncoherentScatteringFunction
        template="""from MMTK import * 
title = 'Incoherent Scattering Function' 
trajectory = ['*TRAJ*',] 
log_file = 'isf_*BTRAJ*' 
time_info = (0, 49999, 1) 
atoms = {'Protein.0': ['Hydrogen']} 
deuter = *DEUTER*
q_vector_set = ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0], 0.100000, 50, None) 
weights = 'incoherent' 
ft_window = 50 
units_q = 1/Units.nm 
time_steps = 25000
frequency_points = 25000
output_files = {'fft': 'ISF_SPECT_*ISFOUT*', 'isf': 'ISF_*ISFOUT*'}"""
        pM=self.solution['T_'+T]['nMoldyn']
        isfdir=joink(pM['simdir'],'isf')
        os.system('/bin/mkdir -p %s'%(isfdir))
        isf_spectrdir=joink(pM['simdir'],'isf_spectr')
        os.system('/bin/mkdir -p %s'%(isf_spectrdir))
        query=template.replace( '*TRAJ*',pM['nc'] )
        query=query.replace( '*BTRAJ*',basename(pM['nc']) )
        query=query.replace( '*ISFOUT*',basename(pM['nc']) )
        if self.isDeutd:
            query=query.replace("*DEUTER*","{'Protein.0': ['Hydrogen']}")
        else:
            query=query.replace("*DEUTER*","None")
        inputf=mktemp() ; open(inputf,'w').write(query)
        data = inputFileRead(inputf)
        os.system('/bin/rm '+inputf)
        #calculate scattering for each atom, store in appropriate file
        remaining=len(data.atoms)

        span=0.02*data.time_steps #selected span, in picoseconds
        window=500                #width for energy resolution in picosec
        data.ft_window=100.0*min(1.0,window/span)
        
        for atom in data.atoms.atomList()[0:1]:
            print 'proccessing %s-%s (%d remaining)'%(self.id,atom.name,
                                                      remaining)
            file='/projects/rubredoxin/out/1fhm_D/T_'+T+\
                  '/nMoldyn/junk_isf.nc'
            atoms=Collection()
            atoms.addObject(atom)
            ISF=isf(data.trajectory,
                    qVect=deepcopy(data.q_vector_set),
                    ncDataFN=file,
                    atoms=atoms,
                    bincoh=1.0,
                    timeInfo=data.time_info,
                    nsteps=data.time_steps)
            fdle='/projects/rubredoxin/out/1fhm_D/T_'+T+\
                  '/nMoldyn/junk_dsf.nc'
            ScatteringFunctionFFT(file,
                                  fdle,
                                  data.frequency_points,
                                  alpha=100.0/data.ft_window)
            remaining=remaining-1

        outfile=file.replace('.nc','.png')
        outtxt=file.replace('.nc','.txt')
        self.printISF(infile=file,outfile=outfile,outtxt=outtxt)
        print outfile+'\n'+outtxt
        outfile=fdle.replace('.nc','.png')
        outtxt=fdle.replace('.nc','.txt')
        self.printDSF(infile=fdle,outfile=outfile,outtxt=outtxt)
        print outfile+'\n'+outtxt
        return False
    

    def nMoldynISF(self,T='300'):
        """calculate intermediate and final structure factors"""
        from copy import deepcopy
        from tempfile import gettempdir, mktemp
        from MMTK import Collection
        from nMoldyn.misc import inputFileRead
        from my_nMoldyn.core import IncoherentScatteringFunction
        from nMoldyn.core import ScatteringFunctionFFT
        isf=IncoherentScatteringFunction
        template="""from MMTK import * 
title = 'Incoherent Scattering Function' 
trajectory = ['*TRAJ*',] 
log_file = 'isf_*BTRAJ*' 
time_info = (0, 49999, 1) 
atoms = {'Protein.0': ['Hydrogen']} 
deuter = *DEUTER*
q_vector_set = ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0], 0.100000, 50, None) 
weights = 'incoherent' 
ft_window = 50 
units_q = 1/Units.nm 
time_steps = 2000 
frequency_points = 2000 
output_files = {'fft': 'ISF_SPECT_*ISFOUT*', 'isf': 'ISF_*ISFOUT*'}"""
        pM=self.solution['T_'+T]['nMoldyn']
        isfdir=joink(pM['simdir'],'isf')
        os.system('/bin/mkdir -p %s'%(isfdir))
        isf_spectrdir=joink(pM['simdir'],'isf_spectr')
        os.system('/bin/mkdir -p %s'%(isf_spectrdir))
        query=template.replace( '*TRAJ*',pM['nc'] )
        query=query.replace( '*BTRAJ*',basename(pM['nc']) )
        query=query.replace( '*ISFOUT*',basename(pM['nc']) )
        if self.isDeutd:
            query=query.replace("*DEUTER*","{'Protein.0': ['Hydrogen']}")
        else:
            query=query.replace("*DEUTER*","None")
        inputf=mktemp() ; open(inputf,'w').write(query)
        data = inputFileRead(inputf)
        os.system('/bin/rm '+inputf)
        #calculate scattering for each atom, store in appropriate file
        remaining=len(data.atoms)
        for atom in data.atoms.atomList():
            print 'proccessing %s-%s (%d remaining)'%(self.id,atom.name,
                                                      remaining)
            file='%s/%04d.nc'%(isfdir,atom.index)
            atoms=Collection()
            atoms.addObject(atom)
            ISF=isf(data.trajectory,
                    qVect=deepcopy(data.q_vector_set),
                    ncDataFN=file,
                    atoms=atoms,
                    bincoh=1.0,
                    timeInfo=data.time_info,
                     nsteps=data.time_steps)
            ScatteringFunctionFFT(file,
                                  '%s/%04d.nc'%(isf_spectrdir,atom.index),
                                  data.frequency_points,
                                  alpha=100.0/data.ft_window)
            remaining=remaining-1
        return False

    def getSIIScF(self,nt=200,T='300',scale='log10',simLength='',task=''):
        """calculate Self-Intermediate Incoherent Scattering Function
        nt:number of time points for which to calculate the function"""
        from scatter.scatter import SIIScF,storeSIIScF,FourierSIIScf,store_sQW
        from amber.amber10 import crd
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        if task=='0-2meV':
            """target frequency range is [0,2]meV.
            We will calculate up to 8meV, in order to avoid cut-off effects
            hv_max=8meV ==>  t_min=4.136meV*psec/2meV~0.5psec
            Since Prod.times['twx']=0.02psec --> skipframes=24
            """
            skipframes=24
            print 'loading 1 frame every ',str(1+skipframes),' frames'
            #dt: time in between loaded frames, in picoseconds
            dt=(1+skipframes)*Prod.times['twx']
            tl=(crdobj.nframe-10)*Prod.times['twx']
            nt=int(tl/dt)
            tt=[]; ww=[]; dw=4.136/(nt*dt) #dw in meV when dt in psec
            for i in range(nt): tt.append(i*dt); ww.append(i*dw)
            qq=[0.23, 0.43, 0.63, 0.83, 1.03, 1.23, 1.43, 1.63, 1.83]
            crdobj.loadframes(skip=skipframes)
            scattering=SIIScF(qq,crdobj.frames,dt,tl,nt)
            avscattering=scattering.sum(axis=0)/len(scattering)
            outf=jk(pH['simdir'],'SIIScF_2meV.dat')
            storeSIIScF(outf,avscattering,qq,tt)
            def convDCS(t): return convT(t,type='corrected')
            sQW=FourierSIIScf(avscattering,tt=tt,resolFun=convDCS)
            outf=jk(pH['simdir'],'ScF_2meV.dat')
            store_sQW(outf,sQW,qq,ww)
            return False
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
        pH['SIIScF'+scale]=crdf.replace('.crd.gz','.SIIScF.'+scale+'.dat')
        open(pH['SIIScF'+scale],'w').write(buf)
        return True

    def get_intSIIScF_sts(self,T='300',simLength='',tl=40,recover=False,
                          task='init',q_i=None,q_f=None,outf=None):

        """calculations related to the integral of the
        Self-Intermediate Incoherent Scattering Function"""
        iq_i=-1; iq_f=-1
        qq=[0.5, 1.0, 1.5, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
        nq=len(qq)
        if q_i and q_f:
            """now do the integration in between these two values"""
            q_i=float(q_i); q_f=float(q_f); iq=0; iq_i=0; iq_f=0
            for q in qq:
                if q==q_i: iq_i=iq
                if q==q_f: iq_f=iq
                iq+=1
        if task=='init':
            """Initialize integrals of the Self-Intermediate
            Incoherent Scattering Function and the structure factor
            (convoluted with the resolution function) for every
            hydrogen. Calculate also the average over all hydrogens
            """
            from scatter.scatter import intSIIScF,storeSIIScF,store_intSIIScF,fourier
            from amber.amber10 import crd
            tl=float(tl)
            Prod=self.solution['T_'+T]['Prod'+simLength]
            pH=self.solution['T_'+T]['onlyH'+simLength]
            skipframes=0
            print '#loading 1 frame every ',str(1+skipframes),' frames'
            #dt: time in between loaded frames, in picoseconds
            dt=(1+skipframes)*Prod.times['twx']
            nt=int(tl/dt); print '#number of time points=',nt
            tt=[]; ww=[]; dw=4.136/(nt*dt) #dw in meV when dt in psec
            for i in range(nt): tt.append(i*dt); ww.append(i*dw)
            crdf=pH['crd.fit2prev']
            crdobj=crd(crdf,self.Nhydr)
            crdobj.loadframes(skip=skipframes)
            scattering=intSIIScF(qq,crdobj.frames,dt,tl,nt)
            #store all intermediate factors
            pH['int_isf']=joink(pH['simdir'],'int_isf')
            os.system('mkdir -p '+pH['int_isf'])
            for i in range(self.Nhydr):
                outfile=joink(pH['int_isf'],'%04d.dat'%(i+1,))
                storeSIIScF(outfile,scattering[i],qq,tt)
            print '#store the average intermediate structure factor'
            pH['intSIIScF_sts']=crdf.replace('.crd.gz','.intSIIScF_sts.dat')
            avscattering=scattering.sum(axis=0)/len(scattering) 
            storeSIIScF(pH['intSIIScF_sts'],avscattering,qq,tt)
            crdobj=None #deallocate memory
            print '#perform Fourier transform'
            gaussianenvelope=[]; factor=1.0/convT(tt[0])
            for tti in tt:
                #gaussianenvelope.append(factor*convT(tti))
                gaussianenvelope.append(1.0) #do not convolve
            for iatom in range(self.Nhydr):
                sc=scattering[iatom]
                for iq in range(nq):
                    sc[iq]=fourier(sc[iq],conv=gaussianenvelope).real
            print '#store all structure factors'
            pH['int_dsf']=joink(pH['simdir'],'int_dsf')
            os.system('mkdir -p '+pH['int_dsf'])
            for i in range(self.Nhydr):
                outfile=joink(pH['int_dsf'],'%04d.dat'%(i+1,))
                storeSIIScF(outfile,scattering[i],qq,ww)            
            print '#store the average structure factor'
            pH['intSq_sts']=crdf.replace('.crd.gz','.intSq_sts.dat')
            avscatteringSq=numpy.zeros(nq*nt).reshape(nq,nt)
            for iq in range(nq): avscatteringSq[iq]=fourier(avscattering[iq]*numpy.array(gaussianenvelope)).real
            storeSIIScF(pH['intSq_sts'],avscatteringSq,qq,ww)
            if iq_i>=0 and iq_f>=0:
                fc=3.0/(q_f**3-q_i**3); rr='%05.2f-%05.2f'%(q_i,q_f)

                print '#integrate intermediate and structure factor'
                print '#in the %s q-range'%rr                

                outf=joink(pH['simdir'],'int_isf/%s-int_isf.dat'%rr)
                isc=fc*(avscattering[iq_f]-avscattering[iq_i])
                store_intSIIScF(outf,isc,tt)
                
                outf=joink(pH['simdir'],'int_dsf/%s-int_dsf.dat'%rr)
                isc=fc*(avscatteringSq[iq_f]-avscatteringSq[iq_i])
                store_intSIIScF(outf,isc,ww)

        elif task=='SIIScF_sts corrected resolution':
            """We recalculate the average structure factor with the
               NIST resolution function"""
            from scatter.scatter import store_intSIIScF,fourier
            #read the integral of the average intermediate struct factor.
            #This is integral from Q=0 to any Q
            inp='ref15-shel.00001_50000.fit2prev.H.intSIIScF_sts.dat'
            pinp = open( inp );
            pinp.readline() #this is just a comment line
            pinp.readline() #the list of Q values
            tt = [ ]
            avscattering =[ '' ] * nq
            for iq in range( nq ): avscattering[ iq ] = [ ]
            for l in pinp.readlines():
                vals = [ float(x) for x in l.split() ]
                tt.append( float( vals[0] ) )
                for iq in range( nq ):
                    avscattering[ iq ].append( vals[ 1+iq ] )
            nt = len( tt )
            dt = 0.02 #time between frames, in picoseconds
            avscattering = numpy.array(avscattering).reshape(nq,nt)
            #prepare resolution function
            factor=1.0/convT( dt, type='corrected' )
            gaussianenvelope=[]
            for tti in tt:
                gaussianenvelope.append(factor*convT(tti,type='corrected'))
            #calculate the average structure factor
            avscatteringSq=numpy.zeros(nq*nt).reshape(nq,nt)            
            for iq in range(nq):
                avscatteringSq[iq]=fourier(avscattering[iq]*numpy.array(gaussianenvelope)).real
            #now calculate the integral between the two Q values passed
            fc=3.0/(q_f**3-q_i**3)
            outf='int_dsf/corrected.dat'
            isc=fc*(avscatteringSq[iq_f]-avscatteringSq[iq_i])
            #prepare list of energies
            ww=[]; dw=4.136/(nt*dt) #dw in meV when dt in psec
            for i in range(nt): ww.append(i*dw)
            #store the integrated structure factor
            store_intSIIScF(outf,isc,ww)
            return False
        
        elif task=='fourier or ISSF':
            """We calculate the fourier transform of the ISSF.dat"""
            from scatter.scatter import store_intSIIScF,fourier
            #read the integral of the average intermediate struct factor.
            #This is integral from Q=0 to any Q
            inp='IISF.dat'
            pinp = open( inp );
            pinp.readline() #this is just a comment line
            pinp.readline() #this is just a comment line
            pinp.readline() #this is just a comment line
            tt = [ ]
            IISF = [ ]
            for l in pinp.readlines():
                vals = [ float(x) for x in l.split() ]
                tt.append( float( vals[0] ) )
                IISF.append( float( vals[1] ) )
            nt = len( tt ) #number of time points
            dt = 0.02 #time between frames, in picoseconds
            #prepare list of energies
            ww=[]; dw=4.136/(nt*dt) #dw in meV when dt in psec
            for i in range(nt): ww.append(i*dw)
            #fourier and store
            ISF = fourier( IISF ).real
            buf0  = '#Incoherent Structure Factor\n'
            buf0 += '#w(ueV) Q_1(A^(-1)) Q_2 ...\n'
            buf0 += '#Q=        3.50000\n'
            buf1 = ''
            buf2 = ''
            for i in range( nt ):
                buf1 += '%9.3f %8.5f\n'%( ww[i], ISF[i] );
                buf2 += '%9.3f %8.5f\n'%( ww[i], abs(ISF[i]) );
            open( 'ISF.dat', 'w' ).write( buf0 + buf1 )
            open( 'ISF_abs.dat', 'w' ).write( buf0 + buf2 )
            return False
        return True
        
    def getSIIScF_sts(self,T='300',simLength='',tl=40,recover=False,
                      task=''):
        """calculate Self-Intermediate Incoherent Scattering Function
        nt:number of time points for which to calculate the function
        tl: total time desired, in picoseconds"""
        from scatter.scatter import SIIScF,storeSIIScF,fourier,tt2ww
        if task=='quasiH harmonic approx':
            """Fourier{ 1-1/6*q^2*MSD^2(t) }, with q=2.0A^{-1}"""
            from utilities.readingWritingFiles import read_to_cols,\
                 write_lists
            qH_p=self.solution['T_'+T]['onlyH'+simLength]['qH']
            assert 'msd2' in qH_p.keys(),\
                   'ERROR getSIIScF_sts, quasiH harmonic '+\
                   'approx: no msd2 info\n'
            isfd=jk(qH_p['sd'],'isf'); qH_p['isf']=isfd #store I(q,t)
            dsfd=jk(qH_p['sd'],'dsf'); qH_p['dsf']=dsfd #store S(q,w)
            os.system('mkdir -p %s %s'%(isfd,dsfd))
            for imode in range(1,1+3*self.Nhydr):
                inpf=jk(qH_p['msd2'],'%04d.dat'%imode)
                cols=read_to_cols(inpf)
                y=[str(1-float(x)/3.0) for x in cols[2]]; cols[2]=y
                outf=jk(qH_p['isf'],'%04d.dat'%imode)
                cmd='#t(psec) I(q=2A^-1,t)=1-q^2*MSD^2/6'
                write_lists(outf,[cols[1],cols[2]],comment_line=cmd)
                inpf=outf
                outf=jk(qH_p['dsf'],'%04d.dat'%imode)
                self.getFourierTransform(inpf,outf)
            return True
        from amber.amber10 import crd
        tl=float(tl)
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']
        nt=int(tl/dt) ; print 'number of time points=',nt
        tt=[]
        for i in range(nt): tt.append(i*dt)
        qq=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] #Angstroms^(-1)
        nq=len(qq)
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        crdobj.loadframes(skip=skipframes)
        scattering=SIIScF(qq,crdobj.frames,dt,tl,nt)
        #store all intermediate factors
        pH['isf']=joink(pH['simdir'],'isf')
        os.system('mkdir -p '+pH['isf'])
        for i in range(self.Nhydr):
            outfile=joink(pH['isf'],'%04d.dat'%(i+1,))
            storeSIIScF(outfile,scattering[i],qq,tt)
        print 'store the average intermediate structure factor'
        pH['SIIScF_sts']=crdf.replace('.crd.gz','.SIIScF_sts.dat')
        storeSIIScF(pH['SIIScF_sts'],
                    scattering.sum(axis=0)/len(scattering),
                    qq,tt)
        crdobj=None #deallocate memory
        print 'perform Fourier transform'
        gaussianenvelope=[]; factor=1.0/convT(tt[0])
        for tti in tt: gaussianenvelope.append(factor*convT(tti))
        for iatom in range(self.Nhydr):
            sc=scattering[iatom]
            for iq in range(nq):
                sc[iq]=fourier(sc[iq],conv=gaussianenvelope).real
        print 'store all structure factors'
        pH['dsf']=joink(pH['simdir'],'dsf')
        os.system('mkdir -p '+pH['dsf'])
        for i in range(self.Nhydr):
            outfile=joink(pH['dsf'],'%04d.dat'%(i+1,))
            storeSIIScF(outfile,scattering[i],qq,tt)
        print 'store the average structure factor'
        pH['Sq_sts']=crdf.replace('.crd.gz','.Sq_sts.dat')
        storeSIIScF(pH['Sq_sts'],
                    scattering.sum(axis=0)/len(scattering),
                    qq,tt2ww(tt))
        return True

    def getAvSIIScF_sts(self,isfdir,outfile,pdbfile=None):

        """get average Self-Intermediate Incoherent Scattering
        Function (SIIScF), by providing a directory containing the
        atomic SIIScF's.

        If we provide the pdbfile, it should contain the incoherent
        scattering lengths in place of the B-factors. In addition, the
        order in which the atoms of the pdbfile show up should be the
        same as the ordering of the atomic SIIScF files contained in
        isfdir
        
        """
        import numpy
        #prepare the structure scattering factors
        bfs=numpy.ones(self.Nhydr);
        if pdbfile:
            i=0
            for line in open(pdbfile).readlines():
                bfs[i]=float(line[60:66]); i+=1
        bfs=(bfs*bfs)/(bfs*bfs).sum()
        #calculate the average
        buf=''
        nq=10 #number of different q-vectors
        nt=2000 #number of time-points
        dt=0.02 #delta-time, in picoseconds between consec time-points
        av=numpy.zeros(nt*nq).reshape(nt,nq)
        #pdb.set_trace()
        i=0
        while i<self.Nhydr:
            isf=joink(isfdir,'%04d.dat'%(i+1,)); #print isf
            it=0
            for line in open(isf).readlines():
                if '#' in line: buf=line
                iq=0
                for q in line.split()[1:]:
                    av[it,iq]=av[it,iq]+float(q)*bfs[i]
                    iq+=1
                if '#' not in line: it+=1
            i+=1
        #write to file
        t=0.0
        for it in range(nt):
            buf+='%7.3f'%(t,)
            for iq in range(nq): buf+=' %8.6f'%(av[it,iq],)
            buf+='\n'; t+=dt
        open(outfile,'w').write(buf)
        return False
        
    def getFourierTransform(self,infile,outfile,col=-1,ss=521.0):
        """Fourier transform
        required first column contains time, in picoseconds
        col: column where to get the data (def: last column, index begins
             at zero)
        ss: temporal resolution, in picoseconds (def: 'BASIS'=521ps)"""
        from copy import copy
        from scatter.scatter import fourier
        col=int(col); ss=float(ss)
        isf=[]; ts=[]
        for line in open(infile,'r').readlines():
            if line[0]=='#': continue #this is a comment line
            items=[float(x) for x in line.split()]
            ts.append(items[0]); isf.append(items[col])
        #find frequency unit (inverse of last time point. Time units
        #assumed to be in picoseconds)
        nu=1.0/float(open(infile,'r').readlines()[-1].split()[0])
        #list of frequencies, in mili electron volts
        factor=4.136 #from 1/picosec to 10E-3eV)
        omegas=[]    #list of energies, in 10E-3eV
        for i in range(len(isf)):
            t=ts[i]
            omegas.append(i*nu*factor)
            isf[i]*=convT(t,ss=ss)/convT(0,ss=ss) #"convolute"
        #the self-intermediate structure factor, averaged over all
        #q-directions, is symmetric with respect to time reversal
        #(I(-t)=I(t)). We cannot do a fourier transform from -T to T,
        #but only from 0 to T. This limitation implies we retain twice
        #the real part of the fourier transform.
        dsf=2*fourier(isf).real
        #print to file
        buf='#fourier transform of '+basename(infile)+'\n'
        for i in range(len(isf)): buf+='%f %f\n'%(omegas[i],dsf[i])
        open(outfile,'w').write(buf)
        return False

    def getVelCorr(self,tl=40,T='300',simLength=''):
        """calculate velocity autocorrelation function"""
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
        nt=int(tl/dt) ; print 'number of time points=',nt
        crdobj.loadframes(skip=skipframes)
        scattering=autoCorr(crdobj.frames,dt,tl,nt)
        #store all velocity autocorrelations
        pH['velAcorr_dir']=joink(pH['simdir'],'velAcorr')
        os.system('mkdir -p '+pH['velAcorr_dir'])
        for i in range(self.Nhydr):
            sc=scattering[i] ; buf='#time vel-autocorr\n'
            for it in range(nt): buf+='%f %f\n'%(it*dt,sc[it])
            outfile=joink(pH['velAcorr_dir'],'%04d.dat'%(i+1,))
            open(outfile,'w').write(buf)
        #average for all atoms
        sc=scattering.sum(axis=0)/self.Nhydr ; buf='#time vel-autocorr\n'
        for it in range(nt): buf+='%f %f\n'%(it*dt,sc[it])
        pH['velACorr']=crdf.replace('.vel.H.gz','.velACorr.dat')
        open(pH['velACorr'],'w').write(buf)
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
        from scatter.scatter import fourier
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
            ss=fourier(ss)
            N=2
            for x in ss[1:]: #avoid first entry, which is delta function
                buf+=str(N*factor)+' '+str(x)+'\n'
                N+=1
            line=pt.readline()
        pH['Sq_sts']=pH['SIIScF_sts'].replace('SIIScF','Sq')
        open(pH['Sq_sts'],'w').write(buf)
        return True

    def get_msd2(self,T='300',tl=40,simLength='',task=''):
        """calculate time-dependent average atomic fluctuation
        tl: total time desired, in picoseconds"""
        from scatter.scatter import msd2
        from amber.amber10 import crd
        tl=float(tl)
        Prod=self.solution['T_'+T]['Prod'+simLength]
        if task=='quasiH':
            from scatter.scatterModes import msd2 #overwrite msd2
            qH_p=self.solution['T_'+T]['onlyH'+simLength]['qH']
            crdf=qH_p['crd']; nd=3*self.Nhydr
            crdobj=crd(crdf,self.Nhydr) #assume modes grouped in triads
            skipframes=0
            #dt: time in between loaded frames, in picoseconds
            dt=(1+skipframes)*Prod.times['twx']
            nt=int(tl/dt) ; print 'number of time points=',nt
            print 'loading 1 frame every ',str(1+skipframes),' frames'
            crdobj.loadframes(skip=skipframes)
            #npart should be self.Nhydr, nxyz should be 3
            nfr,npart,nxyz=crdobj.frames.shape 
            #there's only one coordinate per mode, thus we reshape frames
            scattering=msd2(crdobj.frames.reshape(nfr,npart*nxyz),dt,tl,nt)
            #store all mean fluctuations
            qH_p['msd2']=joink(qH_p['sd'],'msd2')
            os.system('mkdir -p '+qH_p['msd2'])
            for i in range(nd):
                ss=scattering[i] ; buf='#time MSD^2\n'
                for it in range(nt): buf+='%f %f\n'%(dt*it,ss[it])
                outfile=joink(qH_p['msd2'],'%04d.dat'%(i+1,))
                open(outfile,'w').write(buf)
            #store average over all atoms
            ss=scattering.sum(axis=0)/nd
            buf='#time MSD^2\n'
            for it in range(nt): buf+='%f %f\n'%(dt*it,ss[it])
            qH_p['av_msd2']=jk(qH_p['sd'],'msd2.dat')
            open(qH_p['av_msd2'],'w').write(buf)        
            return True
        elif task=='gnmodes':
            g_p=self.solution['T_'+T]['onlyH'+simLength]['gn']
            crdf=g_p['crd']
            crdobj=crd(crdf,self.Nhydr)
            skipframes=0
            #dt: time in between loaded frames, in picoseconds
            dt=(1+skipframes)*Prod.times['twx']
            nt=int(tl/dt) ; print 'number of time points=',nt
            print 'loading 1 frame every ',str(1+skipframes),' frames'

            crdobj.loadframes(skip=skipframes)

            scattering=msd2(crdobj.frames,dt,tl,nt) #scattering[Nhydr][nt]
            #store all mean fluctuations
            g_p['msd2']=joink(g_p['sd'],'msd2')
            os.system('mkdir -p '+g_p['msd2'])
            for i in range(self.Nhydr):
                ss=scattering[i] ; buf='#time MSD^2\n'
                for it in range(nt): buf+='%f %f\n'%(dt*it,ss[it])
                outfile=joink(g_p['msd2'],'%04d.dat'%(i+1,))
                open(outfile,'w').write(buf)
            #store average over all atoms
            ss=scattering.sum(axis=0)/self.Nhydr
            buf='#time MSD^2\n'
            for it in range(nt): buf+='%f %f\n'%(dt*it,ss[it])
            g_p['av_msd2']=jk(g_p['sd'],'msd2.dat')
            open(g_p['av_msd2'],'w').write(buf)        
            return True
        
        #calculate for atoms    
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['crd.fit2prev']        
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']
        nt=int(tl/dt) ; print 'number of time points=',nt
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        crdobj.loadframes(skip=skipframes)
        scattering=msd2(crdobj.frames,dt,tl,nt)
        #store all mean fluctuations
        pH['msd2']=joink(pH['simdir'],'msd2')
        os.system('mkdir -p '+pH['msd2'])
        for i in range(self.Nhydr):
            ss=scattering[i] ; buf='#time MSD^2\n'
            for it in range(nt): buf+='%f %f\n'%(dt*it,ss[it])
            outfile=joink(pH['msd2'],'%04d.dat'%(i+1,))
            open(outfile,'w').write(buf)
        #store average over all atoms
        ss=scattering.sum(axis=0)/self.Nhydr
        buf='#time MSD^2\n'
        for it in range(nt): buf+='%f %f\n'%(dt*it,ss[it])
        pH['av_msd2']=crdf.replace('.crd.gz','.msd2.dat')
        open(pH['av_msd2'],'w').write(buf)        
        return True

    def humpness_msd2(self,T='300',simLength='',task=''):
        """find MSD2 humpness below 1ps for every MSD2(t) plot"""       
        from stats.statistics import linreg
        pH=self.solution['T_'+T]['onlyH'+simLength]
        msd2_dir=joink(pH['simdir'],'msd2')
        pH['humpness']=joink(msd2_dir,'humpness.dat')
        humpf=pH['humpness']; nat=self.Nhydr
        if task=='quasiH':
            qH_p=self.solution['T_'+T]['onlyH'+simLength]['qH']
            msd2_dir=joink(qH_p['sd'],'msd2')
            qH_p['humpness']=jk(msd2_dir,'humpness.dat')
            humpf=qH_p['humpness']; nat=3*self.Nhydr
        elif task=='gnmodes':
            gn_p=self.solution['T_'+T]['onlyH'+simLength]['gn']
            msd2_dir=joink(gn_p['sd'],'msd2')
            gn_p['humpness']=jk(msd2_dir,'humpness.dat')
            humpf=gn_p['humpness']
        buf='#i-hydr atom-# slope intercept corr-coeff  tmax msd2max hump\n'
        humps=[]
        for i in range(nat):
            infile=joink(msd2_dir,'%04d.dat'%(i+1,))
            t=[] ; t0=0.0 ; msd2=[] ; m0=0.0
            xl= 1.0 #1picosec
            xL=10.0 #10picosecs
            for line in open(infile).readlines()[1:]:
                x,y=line.split(); x=float(x) ; y=float(y)
                if x<xl:
                    if y>m0:
                        m0=y ; t0=x
                elif x<xL:
                    t.append(float(x)); msd2.append(float(y))
                else:
                    break #do not store values for t>xL
            #linear regression in the [1.0,10.0]psec time window
            regr=linreg(t,msd2); a=regr['slope']; b=regr['intercept']
            hump=m0-(b+a*t0); humps.append(hump)
            k=self.unsolv['hydr']['atomnos'][i]
            buf+='%3d %4d %f %f %f %f %f %f\n'%(i+1,k,a,b,regr['corrcoef'],
                                            t0,m0,hump)
        pdb.set_trace()
        open(humpf,'w').write(buf)
        
        if not task:
            #create PDB with humpness values as B-factors.
            #Non-H atoms are given a 0.0 humpness
            buf='' ; i=0; minH=min(humps); maxH=max(humps)
            for line in chomp(open(self.unsolv['preminPDB']).readlines()):
                if line[0:4]=='ATOM':
                    l=len(line)
                    if len(line)<66: line=line+' '*(66-l)
                    bf=50.0
                    if line[13]=='H' and i<self.Nhydr:
                        bf=100*(humps[i]-minH)/(maxH-minH) ;i+=1
                    line=line[0:60]+'%5.2f'%(bf,)+line[66:]
                buf+=line+'\n'            
            pH['pdbhump']=joink(pH['msd2'],'humpness.pdb')
            open(pH['pdbhump'],'w').write(buf)
        return True

    def fract_hump(self,T='300',simLength='',task=''):
        """find fraction of H's (D's) important for the hump in the
        MSD^2(t) plot"""
        pH=self.solution['T_'+T]['onlyH'+simLength]
        humpf=pH['humpness']; msd2_dir=pH['msd2']
        if task=='gnmodes':
            humpf=pH['gn']['humpness']; msd2_dir=pH['gn']['msd2']
        #store the humps
        humps={}
        for line in open(humpf).readlines()[1:]:
            l=line.split() ; i=int(l[0]) ; h=float(l[-1]) ; humps[i]=h
        #sort in decreasing humpness
        pairs=[ (h,i) for (i,h) in humps.items()]
        pairs.sort(); pairs.reverse()
        #create average MSD^2 plot by increasing the number of hydrogens
        q=len(open(joink(msd2_dir,'0001.dat')).readlines()[1:])
        n=len(pairs); MSD2=[0.0]*q; bufs=['']*q
        for pair in pairs:
            ii=pair[1]
            ls=open(joink(msd2_dir,'%04d.dat'%(ii,))).readlines()[1:]
            j=0
            for l in ls:
                (t,m)=l.split(); #print j
                MSD2[j]+=float(m)
                j+=1
        for j in range(q):
            MSD2[j]/=n; bufs[j]+=' %5.3f'%(MSD2[j],)
        for pair in pairs[0:-1]:
            ii=pair[1];print ii
            ls=open(joink(msd2_dir,'%04d.dat'%(ii,))).readlines()[1:]
            j=0
            for l in ls:
                (t,m)=l.split()
                MSD2[j]=(MSD2[j]*n-float(m))/(n-1)
                bufs[j]+=' %5.3f'%(MSD2[j],)
                j+=1
            n-=1
        buf='#t MSD2...\n'+'\n'.join(bufs)
        open( joink(msd2_dir,'fract_hump.dat'),'w').write(buf)
        return False

    def fract_boson_peak(self,T='300',simLength='',q=2.0,task='atoms'):
        """find fraction of particles (atoms/modes) important for the
        Boson peak in the structure factor."""
        nd=3*self.Nhydr #number of quasi-harmonic modes
        pH=self.solution['T_'+T]['onlyH'+simLength]
        Prod=self.solution['T_'+T]['Prod'+simLength]
        if task=='atoms':
            """We create structure factor plots by consecutively
            removing one atom, starting with the atom with highest
            MSD2-excess, then continue with atoms having decreasing
            MSD2-excess values."""        
            humpf=pH['humpness'];  dsf_dir=pH['dsf']
            #store the MSD2-excess values
            humps={}
            for line in open(humpf).readlines()[1:]:
                l=line.split() ; i=int(l[0]) ; h=float(l[-1]) ; humps[i]=h
            #sort in decreasing MSD2-excess values
            pairs=[ (h,i) for (i,h) in humps.items()]
            pairs.sort(); pairs.reverse()
            #find column appropriate to selected q value
            qs=[ float(x) for x in \
                 open(joink(dsf_dir,'0001.dat')).readline().split()[1:] ]
            ncol=0 
            while qs[ncol]!=q: ncol+=1
            ncol+=1#time occuppies zeroth column, hence increase index
            #find number of time points, and store time points
            datf=joink(dsf_dir,'0001.dat')
            ts=[ float(l.split()[0]) for l in open(datf).readlines()[1:] ]
            nt=len(ts)
            #start 
            sss=[]; ss=numpy.zeros(nt)
            for pair in pairs:
                datf=joink(dsf_dir,'%04d.dat'%int(pair[1]))
                ls=open(datf).readlines()[1:]
                for il in range(nt):
                    ss[il]+=float(ls[il].split()[ncol])
                sss.append(ss.copy())

            for ip in range(len(sss)): sss[ip]/=(ip+1.0)
            buf='#scattering factor by consecutive addition of particles\n'
            buf+='#with increasing MSD2-excess\n'
            buf+='#time s_1 s_2 ...s_{nparticles}\n'
            for it in range(nt):
                buf+=' %6.3f'%ts[it]
                for ip in range(len(sss)): buf+=' %7.4f'%sss[ip][it]
                buf+='\n'
            outf=joink(dsf_dir,'fract_boson_peak_atoms.dat')
            open(outf,'w').write(buf)
        elif task=='quasiH 0001_1125':
            from amber.amber10 import crd
            qH_p=pH['qH']; qH_p['filtered trajs']=[]
            #store the MSD2-excess values, will serve for ranking of modes
            assert 'humpness' in qH_p.keys(),'ERROR: no humpness file\n'
            humpf=qH_p['humpness'];
            humps={}
            for line in open(humpf).readlines()[1:]:
                l=line.split() ; i=int(l[0]) ; h=float(l[-1]) ; humps[i]=h
            #sort in drecreasing MSD2-excess values
            pairs=[ (h,i) for (i,h) in humps.items()]
            pairs.sort(); pairs.reverse()

            ftdir=jk(qH_p['sd'],'filtered_trajs')
            if not isthere(ftdir): os.system('/bin/mkdir -p '+ftdir)
            
            #add incrementally modes
            cutoff=nd;
            while cutoff<2*nd:
                if cutoff>nd: cutoff=nd
                #find modes we include
                filter=numpy.zeros(nd)
                for i in range(cutoff):
                    (hump,imode)=pairs[i]; #print pairs[i]
                    filter[imode-1]=1 #add this mode. Recall C convention
                #filter trajectory keeping only components from
                #included modes
                outf=jk(qH_p['sd'],'filtered_trajs/0001_%04d.crd'%cutoff)
                print outf
                #for purposses of debugging. Check if outf exists and is
                #kosher
                proceed=True
                if isthere(outf+'.gz'):
                    crdobj=crd(outf+'.gz',self.Nhydr)
                    try: crdobj.loadframes(skip=0)
                    except: crdobj.nframe=0
                    if crdobj.nframe==50000:
                        print 'FINISHED '+str(crdobj.nframe)
                        proceed=False
                    else:
                        print 'ERROR: nfr='+str(crdobj.nframe)
                    crdobj.close(task='remove dump file')
                if proceed: #wrong trajectory
                    self.quasiHforH(T=T,simLength=simLength,filter=filter,
                                    outf=outf,
                                    task='filterTrajectory onlyH')
                #self.quasiHforH(T=T,simLength=simLength,filter=filter,
                #                outf=outf,task='filterTrajectory onlyH')
                qH_p['filtered trajs'].append(outf+'.gz')
                cutoff*=2
            return False#True
        elif task=='quasiH':
            """create structure factors in two ways:

            (1) consecutively adding one mode at a time. 
            Modes are ranked by decreasing MSD2-excess values
            Assumed files containing structure factor for the modes
            only have two columns. First column for energy, second
            column for structure factor with q=2.0A^{-1}

            The way to do it is to add to the current trajectory
            coordinates the contribution from additional modes, then with
            this modified trajectory, recalculate the intermediate
            structure factor, then the structure factor.

            We will add increasing number of modes: 1, 2, 4, 8, 16, 32,
            64, 128, 512, all modes.
            
            (2) consecutively removing one mode at a time from the
                average structure factor.

            The way to do it is to add to the current trajectory
            coordinates the contribution from additional modes, then with
            this modified trajectory, recalculate the intermediate
            structure factor, then the structure factor.

            We will removoe increasing number of modes: 1, 2, 4, 8, 16,
            32, 128, 512, all-1 modes
            """
            from amber.amber10 import crd
            qH_p=pH['qH']; qH_p['filtered trajs']=[]
            #store the MSD2-excess values, will serve for ranking of modes
            assert 'humpness' in qH_p.keys(),'ERROR: no humpness file\n'
            humpf=qH_p['humpness'];
            humps={}
            for line in open(humpf).readlines()[1:]:
                l=line.split() ; i=int(l[0]) ; h=float(l[-1]) ; humps[i]=h
            #sort in drecreasing MSD2-excess values
            pairs=[ (h,i) for (i,h) in humps.items()]
            pairs.sort(); pairs.reverse()

            ftdir=jk(qH_p['sd'],'filtered_trajs')
            if not isthere(ftdir): os.system('/bin/mkdir -p '+ftdir)
            
            #add incrementally modes
            cutoff=1;
            while cutoff<2*nd:
                if cutoff>nd: cutoff=nd
                #find modes we include
                filter=numpy.zeros(nd)
                for i in range(cutoff):
                    (hump,imode)=pairs[i]; #print pairs[i]
                    filter[imode-1]=1 #add this mode. Recall C convention
                #filter trajectory keeping only components from
                #included modes
                outf=jk(qH_p['sd'],'filtered_trajs/0001_%04d.crd'%cutoff)
                print outf
                self.quasiHforH(T=T,simLength=simLength,filter=filter,
                                outf=outf,task='filterTrajectory onlyH')
                qH_p['filtered trajs'].append(outf+'.gz')
                cutoff*=2
                
            #reduce incrementally modes
            cutoff=1
            while cutoff<nd:
                if cutoff>=nd: cutoff=nd-1
                #find modes we include
                filter=numpy.ones(nd)
                for i in range(cutoff):
                    (hump,imode)=pairs[i]; #print pairs[i]
                    filter[imode-1]=0 #remove mode. Recall C convention
                #filter trajectory keeping only components from
                #included modes
                outf=jk(qH_p['sd'],
                        'filtered_trajs/%04d_%04d.crd'%(nd,cutoff))
                print outf                    
                self.quasiHforH(T=T,simLength=simLength,filter=filter,
                                outf=outf, task='filterTrajectory onlyH')
                qH_p['filtered trajs'].append(outf+'.gz')
                cutoff*=2

            #Calculate now structure factors
            from scatter.scatter import SIIScF,storeSIIScF,FourierSIIScf,tt2ww,store_sQW
            from utilities.readingWritingFiles import write_from_numpy
            tl=40; dt=Prod.times['twx']; nt=int(tl/dt); tt=[]
            for i in range(nt): tt.append(i*dt)
            qq=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] #Angstroms^(-1)
            for crdf in qH_p['filtered trajs']:
                crdobj=crd(crdf,self.Nhydr)
                crdobj.loadframes(skip=0)
                #scattering[nxyz][nq][nt]
                scattering=SIIScF(qq,crdobj.frames,dt,tl,nt)
                #average over all hydrogen atoms
                scattering=scattering.sum(axis=0)/len(scattering)
                #output to file
                outf=crdf.replace('.crd.gz','.isf')
                storeSIIScF(outf,scattering,qq,tt)
                #find fourier transform of intermediate structure factor
                sQW=FourierSIIScf(scattering,tt=tt,resolFun=convT_BASIS)
                #ouput to file
                outf=outf.replace('.isf','.sqw')
                store_sQW(outf,sQW,qq,tt2ww(tt))
                crdobj.close(task='remove dump file')
            return True
        elif task=='quasiH harmonic approx':            
            """create structure factors in the quasi-harmonic
            approximation in two ways:
            (1) consecutively removing one mode at a time from the
                average structure factor.
            (2) consecutively adding one mode at a time. 
            Modes are ranked by decreasing MSD2-excess values
            Assumed files containing structure factor for the modes
            only have two columns. First column for energy, second
            column for structure factor with q=2.0A^{-1}
            """
            qH_p=pH['qH']
            #store the MSD2-excess values, will serve for ranking of modes
            assert 'humpness' in qH_p.keys(),'ERROR: no humpness file\n'
            humpf=qH_p['humpness']; dsf_dir=qH_p['dsf']
            humps={}
            for line in open(humpf).readlines()[1:]:
                l=line.split() ; i=int(l[0]) ; h=float(l[-1]) ; humps[i]=h
            #sort in drecreasing MSD2-excess values
            pairs=[ (h,i) for (i,h) in humps.items()]
            pairs.sort(); pairs.reverse()
            #find column appropriate to selected q value
            ncol=1#energy occupies zeroth column, hence increase index
            #find number of energy points, and store
            datf=joink(dsf_dir,'0001.dat')
            ts=[ float(l.split()[0]) for l in open(datf).readlines()[1:] ]
            nt=len(ts)
            #start adding modes, ranked by MSD2-excess
            sss=[]; ss=numpy.zeros(nt)
            for pair in pairs:
                datf=joink(dsf_dir,'%04d.dat'%int(pair[1]))
                ls=open(datf).readlines()[1:]
                for il in range(nt):
                    ss[il]+=float(ls[il].split()[ncol])
                sss.append(ss.copy())
            for ip in range(len(sss)): sss[ip]/=(ip+1.0)
            #write to file
            buf='#scattering factor by consecutive addition of modes\n'
            buf+='#with decreasing MSD2-excess\n'
            buf+='#E(meV) s_1 s_2 ...s_{%s}\n'%len(pairs)
            for it in range(nt):
                buf+=' %6.3f'%ts[it]
                for ip in range(len(sss)): buf+=' %7.4f'%sss[ip][it]
                buf+='\n'
            outf=joink(dsf_dir,'fract_boson_peak_towards_mean.dat')
            open(outf,'w').write(buf)
            #start removing modes from the cummulative
            sss=[ss.copy(),]
            for pair in pairs[:-1]:
                datf=joink(dsf_dir,'%04d.dat'%int(pair[1]))
                ls=open(datf).readlines()[1:]
                ls=open(datf).readlines()[1:]
                for il in range(nt):
                    ss[il]-=float(ls[il].split()[ncol])
                sss.append(ss.copy())
            L=len(sss)
            for ip in range(L): sss[ip]/=(L-ip)
            #write to file
            buf='#scattering factor by consecutive removal of modes\n'
            buf+='#with decreasing MSD2-excess\n'
            buf+='#E(meV) s_%s s_%s ...s_1\n'%(L,L-1)
            for it in range(nt):
                buf+=' %6.3f'%ts[it]
                for ip in range(len(sss)): buf+=' %7.4f'%sss[ip][it]
                buf+='\n'
            outf=joink(dsf_dir,'fract_boson_peak_away_from_mean.dat')
            open(outf,'w').write(buf)
                                          
        return False
    
    def getEISF(self,T='300',simLength='',task=''):
        """calculate Self-Intermediate Incoherent Scattering Function
        nt:number of time points for which to calculate the function"""
        if not task: task=='EISF'
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
        if task=='EISF':
            scattering=EISF(qq,crdobj.frames)
        elif task=='EISF2':
            scattering=EISF2(qq,crdobj.frames)
        elif task=='msd2 from EISF2':
            from scatter.scatter import msd2_from_EISF2
            pH['msd2 from EISF2']=msd2_from_EISF2(qq,crdobj.frame)
            return True
        buf=''
        for (q,s) in scattering: buf+=str(q)+' '+str(s)+'\n'
        key=task
        pH[key]=crdf.replace('.crd.gz','.'+key+'.dat')
        open(pH[key],'w').write(buf)
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
        Prod=self.solution['T_'+T]['Prod'+simLength]
        Prod.rmsCentr()
        return True

    def gnmodes(self,T='300',simLength='',task='init'):
        """calculate quasiharmonic modes in the isotropic gaussian
        network approximation"""
        from networks.GN3D.gn3D import gn3D
        from amber.Prod import Prod
        from amber.amber10 import crd
        if 'onlyH' in task: #work for only hydrogens trajectory
            pH=self.solution['T_'+T]['onlyH'+simLength]
            nat=pH['Nhydr']; gn=gn3D(nat); gn_p=None
            sd=joink(pH['simdir'],'gn'); 
            if task=='onlyH init':
                os.system('mkdir -p '+sd)
                pH['gn']={'sd':sd,}; gn_p=pH['gn']
            else:
                gn_p=pH['gn']
                gn.pickleLoad( open( jk(sd,gn_p['obj']) ) )
                if task=='onlyH projectTrajectory':
                    crdf=pH['crd.fit2prev']; crdobj=crd(crdf,nat)
                    gn.loadM(numpy.ones(nat))  #assume hydrogen-mass=one
                    gn.findVarCovar(crdobj.next)#var-cov matrix
                    crdobj.rewind()            #rewind to first frame
                    pout=open(jk(gn_p['sd'],'crd'),'w')
                    pout.write('#trajectory of gaussian network modes\n')
                    gn.projectTrajectory(crdobj.next,pout)
                    os.system('gzip -f '+jk(sd,'crd'))
                    gn_p['crd']=jk(sd,'crd.gz')#compress trajectory
                elif task=='onlyH chisquare':
                    crdf=jk(sd,gn_p['crd']) #trajectory of modes
                    #note nat==Nhydr, even though this is a trajectory
                    #of nodes and there are 3*Nhydr. It is no problem
                    #if we later do frame.ravel(), since
                    #frame.shape=(Nhydr,3)                    
                    crdobj=crd(crdf,nat)    #trajectory object
                    gn.chisquare(crdobj.next,crdobj.rewind)
                    gn_p['chisq']='\n'.join([str(x) for x in gn.chisq])
            #save updated gn object to file 
            gn_p['obj']=jk(sd,'gn.dump');
            gn.pickleDump( open(gn_p['obj'],'w') )
            return True
             
        Prod=self.solution['T_'+T]['Prod'+simLength];
        strip_p=Prod.strip; nat=strip_p['nat']; gn=gn3D(nat); gn_p=None
        sd=joink(Prod.simdir,'strip/gn');
        if task=='init': #init storage of gnmodes
            os.system('mkdir -p '+sd)
            strip_p['gn']={'sd':sd}; gn_p=strip_p['gn']
        else:
            gn_p=strip_p['gn'];
            gn.pickleLoad( open( jk(sd,gn_p['obj']) ) )
            if task=='projectTrajectory':
                crdf=joink(Prod.simdir,strip_p['fit2centr'])
                crdobj=crd(crdf,nat)#coordinate trajectory object
                gn.loadM(self.masses[0:nat]) #load masses
                gn.findVarCovar(crdobj.next) #mass weighted var-cov matrix
                crdobj.rewind()              #rewind to first frame
                #calculate normal modes and project trajectory onto modes
                pout=open(joink(gn_p['sd'],'crd'),'w')
                pout.write('#trajectory of gaussian network modes\n')
                gn.projectTrajectory(crdobj.next,pout)
                os.system('gzip -f '+jk(sd,'crd'))
                gn_p['crd']='crd.gz'           #compress trajectory
            elif task=='chisquare':
                crdf=jk(sd,gn_p['crd']) #trajectory of modes
                crdobj=crd(crdf,nat)    #trajectory object
                gn.chisquare(crdobj.next,crdobj.rewind)
        #save updated gn object to file
        gn_p['obj']='gn.dump'         #pickle dump file of gn object
        gn.pickleDump( open( jk(sd,gn_p['obj']),'w' ) )
        return True
        

    def electrostatics(self,T='300',simLength='',task='force on FE3'):
        pr=self.solution['T_'+T]['Prod'+simLength]
        if task=='force on FE3':
            script="""#!/bin/bash
            
prodd=_PRODD_
dir=_RSTD_
in=_IN_

/bin/cp $in $prodd/$dir
cd $prodd
echo '#electrostatic force for atom 789 (FE3)' > $dir/force.dat
echo '#counter snapshot fx fy fz' >> $dir/force.dat
n=1
for i in `seq 0 _NSNAP_`;do
  echo \"snapshot $i\"
  sander -O -i $dir/electrostatics.in -p curr/top -c $dir/_ID_.rst.$n -ref curr/ref
  v=`grep \"NUMERICAL, ANALYTICAL FORCES\" mdout -A 3|tail -3|tr -s ' '|cut -d ' ' -f 3|tr '\\n' ' '`
  line=`printf '%04d %05d %13.8f %13.8f %13.8f' $i $n $v`
  echo $line >> $dir/force.dat
  let "n += _JUMP_"
/bin/rm mdout  #clean-up
done"""
            pdir=pr.simdir            
            dict={'_PRODD_':pdir,
                  '_NSNAP_':'4999',  #'999'
                  '_JUMP_' :'10',    #'50'
                  '_RSTD_' :'prev/00001_50000_10_rst', #50_rst'
                  '_TOP_'  :'curr/top',
                  '_REF_'  :'curr/ref',
                  '_IN_'   :jk(pcd,'electrostatics.in'),
                  '_ID_'   :self.id
                  }
            for k,v in dict.items(): script=script.replace(k,v)
            (a,tmpf)=mkstemp(dir='/projects/tmp');
            open(tmpf,'w').write(script)
            os.system('source %s && /bin/rm %s'%(tmpf,tmpf))
            return False
        elif task=='force norm on FE3':
            ff=jk(pr.simdir,'prev/00001_50000_10_rst/force.dat')
            avf=numpy.zeros(3); n=0
            buf=''
            for l in open(ff).readlines():
                if l[0]=='#':
                    if 'counter' in l: l=l.strip()+' norm\n'
                    buf+=l
                else:
                    Z=numpy.array([float(x) for x in l.split()][-3:])
                    avf+=Z; n+=1
                    buf+=l.strip()+' %6.3f\n'%numpy.linalg.norm(Z)
            avf/=n
            buf='#average force components: '+str(avf)+'\n'+buf
            open(ff,'w').write(buf)
            return False
        elif task=='force norm deviation on FE3':
            buf=''; avf=None
            ff=jk(pr.simdir,'prev/00001_50000_50_rst/force.dat')
            for l in open(ff).readlines():
                bf=l
                if 'average force components' in l:
                    m=re.compile('\[([-\d\.\s]+)').search(l).group(1)
                    avf=numpy.array([float(x) for x in m.split()])
                if 'counter' in l:
                    bf=l.strip()+' norm-dev\n'
                if l[0]!='#':
                    Z=numpy.array([float(x) for x in l.split()][2:5])-avf
                    bf=l.strip()+' %6.3f\n'%numpy.linalg.norm(Z)
                buf+=bf
            a,tmpf=mkstemp(dir='/projects/tmp')
            open(tmpf,'w').write(buf)
            print tmpf,str(avf)
            #open(ff,'w').write(buf)
            return False
        elif task=='electrostatic work and power':
            from amber.amber10 import ptraj,evecs,top,crd
            #instantiate topology object
            topf=jk(pr.simdir,'strip/ref15-shel.top')
            topobj=top(topf)
            nat=topobj.ITITL['NATOM']
            #instantiate evecs object and initialize some stuff
            nmode=200 #do calculations only for the first 200 modes
            evecsf=jk(pr.simdir,'strip/evecs.all.dat')
            evecsobj=evecs(evecsf)
            evecsobj.loadModes( nq=nmode )  #don't load all modes
            evecsobj.massProcess( topobj.AMASS ) #pass atomic mass
            sM=evecsobj.sM
            #instantiate trajectory object
            (a,prevcrdf)=mkstemp(dir='/projects/tmp')
            crdf='strip/ref15-shel.00001_50000.fit2centr.crd.gz'
            crdf=jk(pr.simdir,crdf)
            os.system('gunzip -c %s > %s'%(crdf,prevcrdf))
            print prevcrdf
            crdobj=crd(prevcrdf,nat)
            nframe=crdobj.nframe
            #load forces on Fe iron
            ff=jk(pr.simdir,'prev/00001_50000_10_rst/force.dat')
            Fshort=[]
            for line in open(ff).readlines()[3:]:
                f=[float(x) for x in line.split()[2:5]] #force
                Fshort.append( numpy.array(f) )
            #there are only 1000 forces, but nframe=50000. We interpolate
            delta=nframe/len(Fshort)
            pdb.set_trace()
            Fshort.append(Fshort[-1])
            fprev=Fshort[0]; F=[];
            for fnext in Fshort[1:]:
                for i in range(delta):
                    F.append( fprev + ((fnext-fprev)*i)/delta )
                fprev=fnext
            #load first frame
            crdobj.loadframe()
            prevfr=crdobj.frame
            prevF=F[0]
            power=numpy.zeros((1+nmode)*nframe).reshape(1+nmode,nframe)
            work=numpy.zeros((1+nmode)*nframe).reshape(1+nmode,nframe)
            #cycle other frames
            for iframe in range(1,nframe):
                crdobj.loadframe()
                nextfr=crdobj.frame
                nextF=F[iframe]
                #displacement in between consecutive time points
                dfr=nextfr-prevfr
                dX=dfr[-1] #last atom is the Fe atom
                #average force in between consecutive time points
                avF=(prevF+nextF)/2.0
                #calculate elementary work, store in first array element
                P=(dX*avF).sum()
                power[0][iframe]=P
                work[0][iframe]=P+work[0][iframe-1]
                #calculate elementary work due to the first 200 modes
                dfr=dfr.ravel()
                for imode in range(nmode):
                    Q=(evecsobj.qs[imode]).m.ravel() #quasiharmonic mode
                    #{sM*Q} is an orthonormal set in 3N-dimensional space
                    dq=(dfr*sM*Q).sum() #projection of drf onto mode
                    dX=(sM*Q*dq).reshape(nat,3) #coordinate component
                    dX=dX[-1]         #last atom is the Fe atom
                    P=(dX*avF).sum()   #elementary work, or power
                    power[1+imode][iframe]=P
                    work[1+imode][iframe]=P+work[1+imode][iframe-1]
                #update previous frame and force
                prevfr=nextfr
                prevF=nextF 
            #output results
            outdir=jk(pr.simdir,'strip/work_power_10')
            os.system('/bin/mkdir -p '+outdir)
            #calculate average and standar deviations
            pav=power.mean(axis=1)
            wav=work.mean(axis=1)
            psd=power.std(axis=1)
            wsd=work.std(axis=1)
            buf='#NOTE: mode=0 refer to work of total displacement, not a mode\n'
            buf+='#imode      <P>        std(P)       <W>         std(W)\n'
            for imode in range(1+nmode):
                buf+='%3d %12f %12f %12f %12f\n'%(imode,pav[imode],psd[imode],wav[imode],wsd[imode])
            outf=jk(outdir,'work_power.dat')
            open(outf,'w').write(buf)
            #output one file per mode
            for im in range(1+nmode):
                buf='#mode %d\n# frame    Power        Work\n'%(1+im)
                for ifr in range(nframe):
                    buf+='%5d %12f %12f\n'%(ifr,power[im][ifr],
                                            work[im][ifr])
                outf=jk(outdir,'%03d.dat'%(im,))
                open(outf,'w').write(buf)
            #clean-up!
            crdobj.close(task='remove dump file')
            os.system('/bin/rm '+prevcrdf)
            return False
        elif task=='covariance mode work to overall work':
            nmode=200 #do calculations only for the first 200 modes
            inpdir=jk(pr.simdir,'strip/work_power_10')
            #load overall work
            L=open( jk(inpdir,'000.dat') ).readlines()[2:]
            w0=numpy.array( [ float(x.split()[-1]) for x in L ] )
            #now cycle over the modes
            cov = [ ]
            for imode in range(1,nmode+1):
                #load mode work
                L=open( jk(inpdir,'%03d.dat'%imode) ).readlines()[2:]
                w=numpy.array( [ float(x.split()[-1]) for x in L ] )
                cov.append( numpy.cov(w0,w)[0][1] )
            #write covariances
            buf ='#covariance between mode work and overall work\n'
            buf+='#mode covariance\n'
            for imode in range(nmode):
                buf+='%03d %9.3f\n'%(1+imode,cov[imode] )
            open(jk(inpdir,'cov_W0_Wi.dat'),'w').write(buf)
            return False
        elif task=='estimate EP change from average force':
            from amber.amber10 import evecs
            #load average force
            ff=jk(pr.simdir,'prev/00001_50000_50_rst/force.dat')
            buf=open(ff).read()
            m=re.compile('\[([-\d\.\s]+)').search(buf).group(1)
            avf=numpy.array([float(x) for x in m.split()])
            #load mode standard deviations, which are the eigenvalues
            evecsf=joink(pr.simdir,'strip/evecs.all.dat')
            ev=evecs(evecsf)
            #now calculate the EP change for every node
            buf="""#estimate electrost pot change from average force
# for every quasi-harmonic mode
# dq*numpy.cdot(avf,e)
# avf: average electrostatic force vector in the trajectory
# dq : displacement of the quasi-harmonic mode
# e  : vector indicating the direction of the displacement of the Fe atom
#      in the quasi-harmonic mode.
# mode-number EP-change mode-displ\n"""
            for q in ev:     #loop in the quasi-harmonic modes
                dq=1/ev.q.e  #displacement
                e=ev.q.m[-1] #FE vector displacement
                buf+='%04d %11.8f %8.6f\n'\
                      %(q.i,abs(dq*numpy.dot(avf,e)),dq)
            open(jk(pr.simdir,'strip/EP_estimate.dat'),'w').write(buf)
            return False
        return False
                    
    def volume(self,T='100',simLength='',):
        Prod=self.solution['T_'+T]['Prod'+simLength];
        
        
        
            
        
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
    pdb.set_trace()
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

def tview(id):
    
    """viewer of object properties in terminal"""

    cl.loadCAN(id).tview()

def view(id):

    """viewver of object properties"""

    cl.loadCAN(id).view()
    

def foo(*kargs,**kwargs):
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
    
def workDSSP(*kargs,**kwargs):
    if 'combined' in kwargs.keys():
        duty=kwargs['combined']
        if duty=='HB diff map':
            from utilities.readingWritingFiles import write_from_numpy
            from plot.twoD.gnuplot import pm3dgrey
            from plot.twoD.xmgrace import binPlot
            for (m,t) in (('1fhm','ref15-shel'),
                          ('1fhm_D','ref15-shel_D')):
                meso=cl.loadCAN(m).dssp['dobj'].hbm #for mesophilic
                print 'meso='+str(meso.sum())
                ther=cl.loadCAN(t).dssp['dobj'].hbm #for thermophilic
                print 'ther='+str(ther.sum())
                buf=''; diff=numpy.zeros(52*52).reshape(52,52)
                for i in range(51): 
                    for j in range(i+1,52): # j > i
                        #structural alignment #by TMalign results in
                        #meso[i] aligning to #ther[i+1]
                        if ther[i][j]<0 or meso[i+1][j+1]<0:
                            buf+='%d %d\n%d %d\n'%(i,j,j,i)
                            d=ther[i][j]-meso[i+1][j+1]
                            #Upper triangle (j>i) shows HBbonds more
                            #stable in the thermophilic. Lower triange
                            #(j<i) shows HBbonds more stable in the
                            #mesophilic. The point of this plot is to
                            #have a white background, corresponding to
                            #no HBbond
                            if d<0:   diff[i][j]=-1*d
                            elif d>0: diff[j][i]=d
                            #print '%2d %2d %5.2f %5.2f %5.2f'%(i+1,j+1,d,ther[i][j],meso[i+1][j+1])
            
                return False
                outf='%s/analysis/%s_%s_HBnative.dat'%(pjd,m,t)
                open(outf,'w').write(buf)
                outf='%s/analysis/%s_%s_diffHBmap.dat'%(pjd,m,t);

                #binPlot(diff,outf,absmin=0.001)
                pdb.set_trace()
                #print '\n'.join([ str(x) for x in diff.sum(axis=1) ])
                #print 'diff='+str(diff.sum())
    else:
        cl.cycleOverList('workDSSP',*kargs,**kwargs)

def volume(*kargs,**kwargs):
     cl.cycleOverList('volume',*kargs,**kwargs)

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

def relax(*kargs,**kwargs):
    """minimize and relax the premin pdb structure"""
    cl.cycleOverList('relax',*kargs,**kwargs)

def are_orthogonal(*kargs,**kwargs):
    """check if quasiharmonic modes are orthogonal"""
    cl.cycleOverList('are_orthogonal',*kargs,**kwargs)

def quasiharmonic(*kargs,**kwargs):
    """perform quasiharmonic analysis"""
    cl.cycleOverList('quasiharmonic',*kargs,**kwargs)

def quasiHforH(*kargs,**kwargs):
    """obtain quasiharmonic modes with my own library"""
    cl.cycleOverList('quasiHforH',*kargs,**kwargs)

def getNvecs(*kargs,**kwargs):
    cl.cycleOverList('getNvecs',*kargs,**kwargs)

def quasiharmonic_movies(*kargs,**kwargs):
    cl.cycleOverList('quasiharmonic_movies',*kargs,**kwargs)
    
def generate_movies(*kargs,**kwargs):
    cl.cycleOverList('generate_movies',*kargs,**kwargs)
    
def projectCRDontoEVECS(*kargs,**kwargs):
    """For every snapshot, find the projection of each trajectory
    snapshot onto the quasi-harmonic modes"""
    cl.cycleOverList('projectCRDontoEVECS',*kargs,**kwargs)

def deuteratePDB(*kargs,**kwargs):
    """modify preminPDB by hydrogen-deuterium exchange and rewrite
    to parsedPDb"""
    cl.cycleOverList('deuteratePDB',*kargs,**kwargs)

def corr_evecs_rmsd(*kargs,**kwargs):
    cl.cycleOverList('corr_evecs_rmsd',*kargs,**kwargs)

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

def get_masses(*kargs,**kwargs):
    """load atomic masses"""
    cl.cycleOverList('get_masses',*kargs,**kwargs)

def dihedralHist(*kargs,**kwargs):
    """store the dihedral distribution for each atom"""
    cl.cycleOverList('dihedralHist',*kargs,**kwargs)

def getZMatrix(*kargs,**kwargs):
    """calcualate trajectory of Z-Matrix"""
    cl.cycleOverList('getZMatrix',*kargs,**kwargs)
    
def get_mass(*kargs,**kwargs):
    """retrieve mass of the protein alone"""
    cl.cycleOverList('get_mass',*kargs,**kwargs)

def needed_waters(*kargs,**kwargs):
    """print number of needed waters molecules for
    a given mass(H20)-to-mass(protein) ratio"""
    cl.cycleOverList('needed_waters',*kargs,**kwargs)

def prepare_raw_pdb(*kargs,**kwargs):
    """Extract as many waters as needed from the solvated premin
    in order to set the hydration ratio to preselected value"""
    cl.cycleOverList('prepare_raw_pdb',*kargs,**kwargs)
    
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

def get_intSIIScF_sts(*kargs,**kwargs):
    """calculate different scattering-related functions"""
    cl.cycleOverList('get_intSIIScF_sts',*kargs,**kwargs)

def getAvSIIScF_sts(*kargs,**kwargs):
    
    """get average Self-Intermediate Incoherent Scattering
    Function (SIIScF), by providing a directory containing the
    atomic SIIScF's.
    
    If we provide the pdbfile, it should contain the incoherent
    scattering lengths in place of the B-factors. In addition, the
    order in which the atoms of the pdbfile show up should be the
    same as the ordering of the atomic SIIScF files contained in
    isfdir
    """
    cl.cycleOverList('getAvSIIScF_sts',*kargs,**kwargs)
    
def getVelCorr(*kargs,**kwargs):
    """calculate velocity autocorrelation function"""
    cl.cycleOverList('getVelCorr',*kargs,**kwargs)

def getDensStates(*kargs,**kwargs):
    """fourier transform of velACorr_sts"""
    cl.cycleOverList('getDensStates',*kargs,**kwargs)

def getSq_sts(*kargs,**kwargs):
    """fourier transform SIIScF_sts"""
    cl.cycleOverList('getSq_sts',*kargs,**kwargs)

def importProd(*kargs,**kwargs):
    """import Prod structure directory"""
    cl.cycleOverList('importProd',*kargs,**kwargs)

def extractHfromPDB(*kargs,**kwargs):
    """extract hydrogen lines from preminPDB"""
    cl.cycleOverList('extractHfromPDB',*kargs,**kwargs)

def getHinRes(*kargs,**kwargs):
    """make a list of hydrogen numbers in every residue"""
    cl.cycleOverList('getHinRes',*kargs,**kwargs)

def hydrWork(*kargs,**kwargs):
    """create trajectory of kinetic energy"""
    cl.cycleOverList('hydrWork',*kargs,**kwargs)

def Ktraj(*kargs,**kwargs):
    """create trajectory of kinetic energy"""
    cl.cycleOverList('Ktraj',*kargs,**kwargs)

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

def diff_MSD_res(T='360',simLength='10ns'):

    """calculate difference in MSD by residue between mesophilic
    and thermophilic rubredoxin, after structural superposition"""
    from math import sqrt,pi
    for m,t in ( ('1fhm','ref15-shel'), ):
        cm=cl.loadCAN(m); ct=cl.loadCAN(t)
        msd=cm.solution['T_'+T]['Prod'+simLength].simdir
        tsd=ct.solution['T_'+T]['Prod'+simLength].simdir
        mfile=joink(msd,'strip/'+m+'.00001_50000.fit2centr.bfres.dat')
        tfile=joink(tsd,'strip/'+t+'.00001_50000.fit2centr.bfres.dat')
        ml=[]; tl=[]; factor=sqrt(3./8.)/pi
        for line in open(mfile,'r').readlines():
            ml.append(sqrt(float(line.split()[-1]))*factor)
        for line in open(tfile,'r').readlines():
            tl.append(sqrt(float(line.split()[-1]))*factor)
        #recall residue 'i' of meso aligned to 'i+1' of thermo
        for i in range(1,len(ml)-1): print i,ml[i]-tl[i+1]
        

def getResolution(*kargs,**kwargs):
    """create gaussian for resolution
    The original data was taken at 32 10-6, the BASIS data was
    collected at 3 10-6, these are fwhm resolutions"""
    fwhm={'original':32E-6,'BASIS':3E-6}
    z=2.354820045
    sigma={'original':z*32E-6,'BASIS':z*3E-6}

def get_msd2(*kargs,**kwargs):
    """calculate Self-Intermediate Incoherent Scattering Function
    nt:number of time points for which to calculate the function"""
    cl.cycleOverList('get_msd2',*kargs,**kwargs)

def humpness_msd2(*kargs,**kwargs):
    """find MSD2 humpness below 1ps for every MSD2(t) plot"""
    cl.cycleOverList('humpness_msd2',*kargs,**kwargs)

def fract_hump(*kargs,**kwargs):
    """find fraction of H's (D's) important for the hump in the
    MSD^2(t) plot"""
    cl.cycleOverList('fract_hump',*kargs,**kwargs)

def fract_boson_peak(*kargs,**kwargs):
    """find fraction of particles (atoms/modes) important for the
    Boson peak in the structure factor"""
    cl.cycleOverList('fract_boson_peak',*kargs,**kwargs)

def av_msd2(*kargs,**kwargs):
    """average mean displacements"""
    cl.cycleOverList('av_msd2',*kargs,**kwargs)
    
def av_H_msd(*kargs,**kwargs):
    """average mean square displacement of H or D from trajectory
    fitted to centroid"""
    cl.cycleOverList('av_H_msd',*kargs,**kwargs)

def rmsCentr(*kargs,**kwargs):
    cl.cycleOverList('rmsCentr',*kargs,**kwargs)

def gnmodes(*kargs,**kwargs):
    cl.cycleOverList('gnmodes',*kargs,**kwargs)

def electrostatics(*kargs,**kwargs):
    cl.cycleOverList('electrostatics',*kargs,**kwargs)

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

def checknMoldynISF(*kargs,**kwargs):
    """calculate intermediate and final structure factors"""
    cl.cycleOverList('checknMoldynISF',*kargs,**kwargs)

def nMoldynISF(*kargs,**kwargs):
    """calculate intermediate and final structure factors"""
    cl.cycleOverList('nMoldynISF',*kargs,**kwargs)

def findWeights(*kargs,**kwargs):
    """find incoherent scattering weights"""
    cl.cycleOverList('findWeights',*kargs,**kwargs)

def printISF(*kargs,**kwargs):
    cl.cycleOverList('printISF',*kargs,**kwargs)

def printDSF(*kargs,**kwargs):
    cl.cycleOverList('printDSF',*kargs,**kwargs)

def averageISF(*kargs,**kwargs):
    cl.cycleOverList('averageISF',*kargs,**kwargs)

def getFourierTransform(*kargs,**kwargs):
    cl.cycleOverList('getFourierTransform',*kargs,**kwargs)

def serveKenSq_sts(inf='',outf=''):
    """reformat *fit2prev.H.Sq_sts.dat to give it to Ken"""
    #read input file
    #pdb.set_trace()
    pin=open(inf); pin.readline();
    qq=pin.readline().split()[1:]
    ww=[]; ss=[];
    for i in range(len(qq)): ss.append([])
    for l in pin.readlines():
        print l.strip()
        items=l.split(); ww.append(items[0])
        for i in range(len(qq)): ss[i].append(items[1+i])
    #write output file
    buf=''
    for iq in range(len(qq)):
        buf+='# Q=%sA^-1\n#\n# E(meV) S(q,w)\n'%qq[iq]
        for iw in range(len(ww))[0:968]: buf+='%s %s\n'%(ww[iw],ss[iq][iw])
        buf+='&\n'
    open(outf,'w').write(buf)

def diffusion(task=None,inp=None,inp2=None,out=None):
    if task=='estimate diffusion of hydrogens':
        """estimate the diffusion coefficient by calculating the
        average of the MSD in the last 100ps and then dividing by the
        simulation lenght"""
        from utilities.readingWritingFiles import read_to_cols
        a = read_to_cols(inp, outFmt='numpy.array', xtype='float')
        a = a.transpose()
        #start at iat=2, since the first two columns are for time and
        #the MSD resulting for averaging the MSD of all
        #hydrogens. Based on this average, select chunk 15000:20000
        #pdb.set_trace()
        T = inp2 #temperature
        factors = {'260':{'f':1./35, 'i':15000, 'e':20000},
                   '300':{'f':1./15, 'i':2, 'e':7520},
                   '360':{'f':1./17, 'i':2, 'e':8372}
                   }
        f = factors[T]['f']; i=factors[T]['i']; e=factors[T]['e']
        buf = '# hydrogen D_t'
        for iat in range(2,len( a )):
            msd = a[iat]
            dc = f * numpy.average(msd[i:e])#units=10^(-5)cm^2/s
            buf += '%3d %6.3f\n'%(iat-1,dc)
        open('diffusionfit2prevH.dat', 'w').write( buf )
    elif task=='find diffussion cutoff':
        from utilities.readingWritingFiles import read_to_cols
        T = inp2
        factors = { '260':{'f':3.3611, 'co':0.21},
                    '300':{'f':1.9694, 'co':0.04},
                    '360':{'f':3.4150, 'co':0.02}
                    }
        f = factors[T]['f']; co = factors[T]['co']
        a = read_to_cols(inp, outFmt='numpy.array', xtype='float')
        a = a.transpose()

        D = numpy.array( a[1] ) #diffussion coefficients
        #sort from big to small
        perm = numpy.argsort(D)[::-1]
        D = f*D[perm]
        #pdb.set_trace()
        #start adding up
        dd = 0.; n = 1
        pdb.set_trace()
        for d in D:
            if d < co:
                print dd/n
                return False
            dd += d
            n += 1
    return False

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
