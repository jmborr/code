#!/usr/bin/python

import os,sys,re,pdb,numpy,pickle
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir
from utilities.small_utilities import chomp,Bye,junkName
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
joink=os.path.join; jk=joink
isthere=os.path.exists
basename=os.path.basename

my_project='eflow'                   #substitute 'my_project'
projectsdir='/projects'
pcd=joink(codedir,'projects',my_project)  #code directory
pjd=joink(projectsdir,my_project)         #projects directory
pjdout=joink(projectsdir,my_project,'out')
collf=joink(pjd,'varsets.can')
canlistf=joink(pjd,'out','list')
repod=joink(pjd,'repod')

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
        Prod=self.__dict__['T_270']['10ns']['Prod']
        del Prod.strip['crdfile.res']
        return True

    def gensimdir(self):
        """simulation directory"""
        self.simdir=joink(pjdout,self.id)
        if not isthere(self.simdir):
            os.system('/bin/mkdir -p '+self.simdir)
        return True

    def storeL(self):
        """store protein length"""
        self.L=53
        return True

    def importSomeFiles(self):
        """import PDB, topology, and reference files"""
        #solvated configuration
        self.parmfile=joink(self.simdir,self.id+'.top')
        self.crdfile=joink(self.simdir,self.id+'.crd')
        self.preminPDB=joink(self.simdir,self.id+'.premin.pdb')
        #dry configuration
        self.unsolv={}
        self.unsolv['parmfile']=joink(self.simdir,self.id+'.uns.top')
        self.unsolv['crdfile']=joink(self.simdir,self.id+'.uns.crd')
        self.unsolv['preminPDB']=joink(self.simdir,
                                       self.id+'.uns.premin.pdb')
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

    def genPDBid(self):
        self.pdbid=self.id[0:4]
        return True

    def getAtomSymbol(self):
        """Store atomic symbol"""
        self.atomSymbol=[]
        for x in open(self.preminPDB).readlines():
            if x[0:5]=='ATOM ': self.atomSymbol.append( x[13] )
        return True
    
    def neighborlist(self):
        """find neighbor atom contact list"""
        from pdbfiles.smallThings import neighborList
        #find native contact map
        neig,dist,xyz=neighborList(open(self.unsolv['preminPDB']),co=6.0)
        self.contmap={'all':{'co':6.0,'map':neig,'distances':dist}}
        return True

    def qH(self,T='270',subdir='',task='all atoms'):
        """use my quasi harmonic library"""
        from amber.Prod import Prod
        from quasiH.quasi import modes
        from amber.amber10 import crd
        Prod=self.__dict__['T_'+T]['Prod']
        if subdir: Prod=self.__dict__['T_'+T][subdir]['Prod']
        nat=Prod.strip['nat']
        if 'all atoms' in task: #work for only hydrogens trajectory
            nd=nat*3; qH=modes(nd); qH_p=None
            sd=jk(Prod.simdir,'strip','qH');
            if 'init' in task:
                os.system('mkdir -p '+sd)
                Prod.strip['qH']={'sd':sd}; qH_p=Prod.strip['qH']
            else:
                qH_p=Prod.strip['qH']
                qH.pickleLoad( open(jk(sd,qH_p['obj'])) )
                if 'findVarCovar' in task:
                    M=[]
                    for m in self.masses[0:nat]: M+=[m]*3
                    qH.loadM(M)
                    crdf=jk(Prod.simdir,Prod.strip['fit2centr'])
                    crdobj=crd(crdf,nat)
                    qH.findVarCovar(crdobj.next)
        #save updated object to file 
        qH_p['obj']=jk(sd,'qH.dump');
        qH.pickleDump( open(qH_p['obj'],'w') )
        return True
        
    def quasiharmonic(self,T='270',subdir='',nvecs=1000,usevel=False):
        """obtain quasiharmonic modes for CA, heavy-atoms,
        hydrogen-atoms, and all-atoms
        usevel: modes for velocities, instead of coordinates"""
        from amber.Prod import Prod
        from amber.amber10 import ptraj
        pr=self.__dict__['T_'+T]['Prod']
        if subdir: pr=self.__dict__['T_'+T][subdir]['Prod']
        filetype='fit2centr'
        if usevel: filetype='velfile' #vels, not coords!
        script=junkName() ; traj=junkName()
        os.system('gunzip -c %s > %s'%\
                  (joink(pr.simdir,pr.strip[filetype]),traj) )
        masks={'CA':'@CA','heavyat':'!@H*','H':'@H*','all':''}
        for (maskk,maskv) in masks.items():
            print 'analyzing quasiharmonic modes for '+maskk+' atoms'
            junkf=junkName()
            open(script,'w').write("""trajin %s 1 50000 5
matrix mwcovar name mwcovar %s out mwcovar.dat byatom
analyze matrix mwcovar out evecs.dat thermo vecs %d
analyze matrix mwcovar out %s thermo vecs 0"""%(traj,maskv,nvecs,junkf))
            cmd='ptraj10 %s < %s > eigen.dat '%\
                 (self.unsolv['parmfile'],script)
            os.system(cmd+' && /bin/rm  '+junkf)
            mw='strip/mwcovar.'+filetype+'.'+maskk+'.dat'
            eg='strip/evecs.'  +filetype+'.'+maskk+'.dat'
            ev='strip/eigen.'  +filetype+'.'+maskk+'.dat'
            os.system( 'mv mwcovar.dat '+joink(pr.simdir,mw) )
            os.system( 'mv evecs.dat '  +joink(pr.simdir,eg)  )
            os.system( 'mv eigen.dat '  +joink(pr.simdir,ev)  )
            pr.strip['mwcovar.'+filetype+'.'+maskk]=mw
            pr.strip['evecs'   +filetype+'.'+maskk]=eg
            pr.strip['eigen'   +filetype+'.'+maskk]=ev
        pr.dump() #store these dictionary items
        os.system('/bin/rm %s %s '%(traj,script))
        return True

    def quasiharmonicCovarAll(self,T='270',subdir=''):
        """obtain all quasiharmonic modes
        usevel: modes for velocities, instead of coordinates"""
        from amber.amber10 import ptraj
        Prod=self.__dict__['T_'+T]['Prod']
        if subdir: Prod=self.__dict__['T_'+T][subdir]['Prod']
        filetype='fit2centr'
        script=junkName() ; traj=junkName()
        os.system('gunzip -c %s > %s'%\
                  (joink(Prod.simdir,Prod.strip[filetype]),traj) )
        masks={'all':'',}; nvecs=3*Prod.strip['nat'] #all eigenvalues!!!
        for (maskk,maskv) in masks.items():
            print 'analyzing quasiharmonic modes for '+maskk+' atoms'
            junkf=junkName()
            open(script,'w').write("""trajin %s 1 50000 5
matrix covar name covar %s out covar.dat byatom
analyze matrix covar out evecs.dat thermo vecs %d
analyze matrix covar out %s thermo vecs 0"""%(traj,maskv,nvecs,junkf))
            cmd='ptraj10 %s < %s > eigen.dat '%\
                 (self.unsolv['parmfile'],script)
            os.system(cmd+' && /bin/rm  '+junkf)
            mw='strip/covar.'+filetype+'.'+maskk+'.dat'
            eg='strip/evecs.'  +filetype+'.'+maskk+'.dat'
            ev='strip/eigen.'  +filetype+'.'+maskk+'.dat'
            os.system( 'mv covar.dat '+joink(Prod.simdir,mw) )
            os.system( 'mv evecs.dat '  +joink(Prod.simdir,eg)  )
            os.system( 'mv eigen.dat '  +joink(Prod.simdir,ev)  )
            Prod.strip['covar.'+filetype+'.'+maskk]=mw
            Prod.strip['evecs'   +filetype+'.'+maskk]=eg
            Prod.strip['eigen'   +filetype+'.'+maskk]=ev
        Prod.dump() #store these dictionary items
        os.system('/bin/rm %s %s '%(traj,script))
        return True

    def genIsDeutd(self):
        """flag as deuterated"""
        self.isDeutd=False
        if self.id in ('dhp__D',): self.isDeutd=True
        return True

    def importProd(self,T='300',subdir='',recover=False):
        """import Prod structure directory"""
        from amber.Prod import Prod
        if recover: recover=True
        TT='T_'+T
        if TT not in self.__dict__.keys(): self.__dict__[TT]={}
        prodir=joink(self.simdir,TT,subdir,'Prod')
        pr=Prod(self.id,prodir,module_name='amber.Prod')
        if not recover: 
            pr.importStripTopCrd(self.unsolv['parmfile'],
                                 self.unsolv['crdfile'])
            pr.count_frames(mode='unsolvated') #check number of frames
            #several RMS calculations to check consistency of simulation
            #pr.getAvg()       #should already be done in jaguar
            pr.rmsJump(outcrd=True)#consec rms and output fitted trajectory
            pr.rmsRef()        #rms with respect to reference
            pr.centr_and_avg() #find centroid
            pr.rmsCentr()      #rms with respect to centroid
            pr.rmsFirst()      #rms with respect to first frame
            pr.fit2centr()     #fit trajectory to centroid
        if subdir:
            self.__dict__[TT][subdir]={'Prod':pr}
        else:
            self.__dict__[TT]['Prod']=pr #store
        return True

    def findResIndexes(self):
        """for each residue, find the list of atom indexes
        NOTE: first atom index is 0, NOT 1
        NOTE: we create and extra residue with the 24 catalytic site
              atoms right before the Na+ ions. We remove the 24 catalytic
              site atoms from the DHP and NDP residues"""
        nres=162 #protein residues + DHP + NPD + catalytic-site
        ptin=open(self.preminPDB,'r'); l=ptin.readline()
        resindexes=[]; currires=None; catalytic=False
        while l:
            if l[0:5]=='ATOM ':
                iat=int(l[6:11])-1; ires=l[22:26]; resname=l[17:20]
                if not currires:
                    currires=ires #initialize residue index
                    group=[] #initialize list of atom residues
                if ires!=currires:
                    resindexes.append(group)
                    currires=ires
                    group=[]
                    if resname=='Na+' and not catalytic:
                        """introduce the catalytic site residue"""
                        for index in self.hot_spot['indexes']:
                            group.append(index)
                        resindexes.append(group); catalytic=True; group=[];
                if iat not in self.hot_spot['indexes']: group.append(iat)
            l=ptin.readline()        
        resindexes.append(group) #enter last group
        self.resinfo={'indexes':resindexes,'nres':nres}
        if 'resinidexes' in self.__dict__.keys():
            del self.__dict__['resinidexes']
        return True

    def keepHeavy(self,T='300',subdir='1ns'):
        """generate trajectory containing only the heavy atoms"""
        print 'generate trajectory containing only the heavy atoms'
        from amber.amber10 import top,crd,ptraj
        t_p=self.__dict__['T_'+T]
        if subdir: t_p=self.__dict__['T_'+T][subdir]
        Prod=t_p['Prod']
        #create appropriate subdirectory and dictionary
        simdir=jk(self.simdir,'T_'+T,subdir,'ha');
        os.system('mkdir -p '+simdir)
        t_p['ha']={'sd':simdir}; ha_p=t_p['ha']
        #find masses
        M=[]
        for m in self.masses[0:Prod.strip['nat']]:
            if m>2.3: M.append(m)
        ha_p['nat']=len(M); ha_p['masses']=numpy.array(M)
        #retrieve heavy atom coordinate trajectory
        crdf=joink(Prod.simdir,Prod.strip['fit2centr'])
        outf=jk(ha_p['sd'],'crd')
        ptr=ptraj(self.unsolv['parmfile'],crdf)
        ptr.strip(entities=['@H*',]) #strip all hydrogens
        ptr.go(outcrd=outf,outfile='./junk.log')
        os.system('gzip -f '+outf); ha_p['crd']=outf+'.gz'
        #retrieve heavy atom velocity trajectory
        crdf=joink(Prod.simdir,Prod.strip['velfile'])
        outf=jk(ha_p['sd'],'vel')
        ptr=ptraj(self.unsolv['parmfile'],crdf)
        ptr.strip(entities=['@H*',]) #strip all hydrogens
        ptr.go(outcrd=outf,outfile='./junk.log')
        os.system('gzip -f '+outf);  ha_p['vel']=outf+'.gz'
        return True

    def keepAll(self,T='300',subdir='1ns'):
        """mimic keepHeavy construct"""
        t_p=self.__dict__['T_'+T]
        if subdir: t_p=self.__dict__['T_'+T][subdir]
        pdb.set_trace()
        Prod=t_p['Prod']
        simdir=jk(Prod.simdir,'strip')
        t_p['aa']={'sd':simdir,
                   'nat':Prod.strip['nat'],
                   'masses':self.masses[0:Prod.strip['nat']],
                   'crd':jk(Prod.simdir,Prod.strip['fit2centr']),
                   'vel':jk(Prod.simdir,Prod.strip['velfile']),
                   }
        return True
        
    def findHAIndexes(self):
        """list of heavy atom indexes"""
        self.haidx=[]
        for l in open(self.unsolv['preminPDB']).readlines():
            if l[:5]=='ATOM ' and l[13]!='H':
                self.haidx.append(int(l[7:11]))
        return True

    def gnmodes(self,T='300',subdir='1ns',task='heavy'):
        """calculate quasiharmonic modes in the isotropic gaussian
        network approximation"""
        from networks.GN3D.gn3D import gn3D
        from amber.Prod import Prod
        from amber.amber10 import crd
        #resolve subdir option
        t_p=self.__dict__['T_'+T]
        if subdir: t_p=self.__dict__['T_'+T][subdir]
        Prod=t_p['Prod']
        if 'all atoms' in task: #modes for heavy atoms
            aa_p=Prod.strip; nat=aa_p['nat']; gn=gn3D(nat); gn_p=None
            sd=jk(Prod.simdir,'strip','gn');
            if 'init' in task:    #initialize modes directory
                os.system('mkdir -p '+sd);
                aa_p['gn']={'sd':sd}; gn_p=aa_p['gn']
                aa_p['crd']=jk(Prod.simdir,aa_p['fit2centr']);
                aa_p['masses']=self.masses[0:nat]
            else:                 #load initialized modes info
                gn_p=aa_p['gn']
                gn.pickleLoad( open( jk(sd,gn_p['obj']) ) )
                if 'findVarCovar' in task:
                    crdf=aa_p['crd']
                    crdobj=crd(crdf,nat)#coordinate trajectory object
                    gn.loadM(aa_p['masses']) #load masses
                    gn.findVarCovar(crdobj.next)#mass weighted varcov matr
                if 'findForceConstants' in task:
                    crdf=aa_p['crd']
                    crdobj=crd(crdf,nat)#coordinate trajectory object
                    gn.loadM(aa_p['masses']) #load masses
                    gn.findVarCovar(crdobj.next)#mass weighted varcov matr
                    gn.findModes()
                    gn.findForceConstants()
                elif 'projectTrajectory' in task:
                    crdf=aa_p['crd']
                    crdobj=crd(crdf,nat)#coordinate trajectory object
                    gn.loadM(aa_p['masses']) #load masses
                    gn.findVarCovar(crdobj.next)#mass weighted varcov matr
                    crdobj.rewind()             #rewind to first frame
                    #calculate modes and project trajectory onto modes
                    pout=open(joink(gn_p['sd'],'crd'),'w')
                    pout.write('#trajectory of gaussian network modes\n')
                    gn.projectTrajectory(crdobj.next,pout)
                    os.system('gzip -f '+jk(sd,'crd'))
                    gn_p['crd']='crd.gz'           #compress trajectory
                elif 'chisquare' in task:
                    crdf=jk(sd,gn_p['crd']) #trajectory of modes
                    crdobj=crd(crdf,nat)    #trajectory object
                    gn.chisquare(crdobj.next,crdobj.rewind)
        elif 'heavy atoms' in task: #modes for heavy atoms
            ha_p=t_p['ha']; nat=ha_p['nat']; gn=gn3D(nat); gn_p=None
            sd=jk(ha_p['sd'],'gn');
            if 'init' in task:    #initialize modes directory
                os.system('mkdir -p '+sd);
                ha_p['gn']={'sd':sd}; gn_p=ha_p['gn']
            else:                 #load initialized modes info
                gn_p=ha_p['gn']
                gn.pickleLoad( open( jk(sd,gn_p['obj']) ) )
                if 'findForceConstants' in task:
                    crdf=joink(ha_p['crd'])
                    crdobj=crd(crdf,nat)#coordinate trajectory object
                    gn.loadM(ha_p['masses']) #load masses
                    gn.findVarCovar(crdobj.next)#mass weighted varcov matr
                    gn.findModes()
                    gn.findForceConstants()
                elif 'projectTrajectory' in task:
                    crdf=joink(ha_p['crd'])
                    crdobj=crd(crdf,nat)#coordinate trajectory object
                    gn.loadM(ha_p['masses']) #load masses
                    gn.findVarCovar(crdobj.next)#mass weighted varcov matr
                    crdobj.rewind()             #rewind to first frame
                    #calculate modes and project trajectory onto modes
                    pout=open(joink(gn_p['sd'],'crd'),'w')
                    pout.write('#trajectory of gaussian network modes\n')
                    gn.projectTrajectory(crdobj.next,pout)
                    os.system('gzip -f '+jk(sd,'crd'))
                    gn_p['crd']='crd.gz'           #compress trajectory
                elif 'chisquare' in task:
                    crdf=jk(sd,gn_p['crd']) #trajectory of modes
                    crdobj=crd(crdf,nat)    #trajectory object
                    gn.chisquare(crdobj.next,crdobj.rewind)
        #save updated gn object to file
        gn_p['obj']='gn.dump'         #pickle dump file of gn object
        gn.pickleDump( open( jk(sd,gn_p['obj']),'w' ) )
        return True


    def resCM(self,T='270',subdir='10ns'):
        """calculate trajectory for coordinates and velocities of the
        center of mass of each residue"""
        from amber.amber10 import crd
        Prod=self.__dict__['T_'+T]['Prod']
        if subdir: Prod=self.__dict__['T_'+T][subdir]['Prod']
        nat=Prod.strip['nat']; nres=self.resinfo['nres']
        resindexes=self.resinfo['indexes'][0:nres]
        masses=self.masses[0:nat]
        #output coordinates trajectory file
        Prod.strip['fit2centr.res']=joink('strip','residues.crd')
        outf=joink(Prod.simdir,Prod.strip['fit2centr.res'])
        
        #input coordinates trajectory file.
        crdf=joink(Prod.simdir,Prod.strip['fit2centr'])
        #produce the output trajectory file
        crdobj=crd(crdf,nat)       
        ptout=open(outf,'w'); ptout.write('nres %d\n'%(nres,))
        results=crdobj.coalesce(ptout,resindexes,masses)
        os.system('gzip -f '+outf); Prod.strip['fit2centr.res']+='.gz'
        

        #input velocities trajectory file
        velf=joink(Prod.simdir,Prod.strip['velfile'])
        #output velocities trajectory file
        Prod.strip['velfile.res']=joink('strip','residues.vel')
        outf=joink(Prod.simdir,Prod.strip['velfile.res'])
        #produce the output trajectory file
        velobj=crd(velf,nat)
        ptout=open(outf,'w'); ptout.write('nres %d\n'%(nres,))
        results=velobj.coalesce(ptout,resindexes,masses)
        os.system('gzip -f '+outf); Prod.strip['velfile.res']+='.gz'
        
        #store residue masses
        self.resinfo['masses']=results['M']
        Prod.dump()
        return True

    def clusterize(self,T='270',subdir='10ns'):
        """cluster residue trajectory snapshots"""
        from amber.amber10 import crd
        Prod=self.__dict__['T_'+T]['Prod']
        if subdir: Prod=self.__dict__['T_'+T][subdir]['Prod']
        nres=self.resinfo['nres']
        crdf=joink(Prod.simdir,Prod.strip['fit2centr.res'])
        crdobj=crd(crdf,nres)
        crdobj.findCentroid()

        
    def importHotSpotInfo(self,massf=None,indexf=None):
        """import hot spot properties"""
        import numpy
        if 'hot_spot' not in self.__dict__.keys(): self.hot_spot={}
        #import PDB indexes for atoms
        if indexf:
            self.hot_spot['indexes']=[]
            index=0 #indexes start with 0, instead of 1
            for line in open(indexf).readlines():
                if '1' in line: self.hot_spot['indexes'].append(index)
                index+=1
        #import masses
        if massf:
            M=[]
            allM=open(massf).read().replace('\n',' ').split()
            for i in range(len(allM)): allM[i]=float(allM[i])
            for index in self.hot_spot['indexes']: M.append(allM[index])
            self.hot_spot['masses']=numpy.array(M)
            self.masses=numpy.array(allM)
        return True

    def printK(self,T='300',subdir=None):

        """obtain kinetic energy plot of hot spot, random set of
        atoms, and all atoms"""
        
        from amber.amber10 import crd
        from random import random
        from stats.statistics import average,standard_dev
        import numpy
        #K2T=9.654E-4 #from amber kinetic energy units to degrees
        K2T=199.1
        TT='T_'+T ;
        if subdir:
            Prod=self.__dict__[TT][subdir]['Prod']
        else:
            Prod=self.__dict__[TT]['Prod']
            
        Prod.strip['Kf']=joink('strip','K.dat') #file to store K data
        
        Natoms=len(self.hot_spot['indexes'])
        nat=Prod.strip['nat']
        randomIdx=[] ; rmass=[]
        for i in range(Natoms):
            index=int(random()*nat)
            randomIdx.append(index)
            rmass.append(self.masses[index])
        rmass=numpy.array(rmass)
        dt=Prod.times['twx'] ; t=0 ; N=0
        buf='#time(psec) T(K) randomT allT\n'
        velf=joink(Prod.simdir,Prod.strip['velfile'])
        velobj=crd(velf,nat)
        Ts={'hot-spot':[],'random':[],'all':[]}
        N=0
        while velobj.loadframe():
            N+=1 ; t+=dt ; v2=[] ; rv2=[] #;print 'N=',N
            for index in self.hot_spot['indexes']:
                v2.append((velobj.frame[index]*velobj.frame[index]).sum())
            K=numpy.dot(self.hot_spot['masses'],numpy.array(v2))
            T=K*K2T/Natoms ; Ts['hot-spot'].append(T)
            for index in randomIdx:
                rv2.append((velobj.frame[index]*velobj.frame[index]).sum())
            rK=numpy.dot(rmass,numpy.array(rv2))
            rT=rK*K2T/Natoms ; Ts['random'].append(rT)
            aK=(self.masses[0:nat]*\
                (velobj.frame*velobj.frame).sum(axis=1)).sum()
            aT=aK*K2T/nat ; Ts['all'].append(aT)
            buf+='%lf %lf %lf %lf\n'%(t,T,rT,aT)
        open(joink(Prod.simdir,Prod.strip['Kf']),'w').write(buf)
        print '%d frames read'%N
        print 'Kinetic energies in '+Prod.strip['Kf']
        print 'averages: T=%lf rT=%lf aT=%lf'%(\
            average(Ts['hot-spot']), average(Ts['random']),\
            average(Ts['all']) )
        print 'standar devs: T=%lf rT=%lf aT=%lf'%\
              (standard_dev(Ts['hot-spot']),standard_dev(Ts['random']),\
               standard_dev(Ts['all']) )
        Prod.dump()
        return True

    def Ktraj(self,T='270',subdir=None):
        """create trajectory of kinetic energy"""
        from amber.amber10 import crd
        from pdbfiles.smallThings import insertBfact
        from utilities.readingWritingFiles import write_from_numpy
        Prod=self.__dict__['T_'+T]['Prod']
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

    def Kcov(self,T='270',subdir=None,mask='all'):
        """covariance matrix of kinetic energy"""
        from amber.amber10 import crd
        from utilities.readingWritingFiles import write_from_numpy
        Prod=self.__dict__['T_'+T]['Prod']
        if subdir: Prod=self.__dict__['T_'+T][subdir]['Prod']
        results=None; rootkey=None; npart=0
        if mask=='all':
            rootkey='Kcov'
            velf=joink(Prod.simdir,Prod.strip['velfile']);
            npart=Prod.strip['nat']; velobj=crd(velf,npart)
            results=velobj.Kcov(self.masses[0:npart])
        elif mask=='res':
            rootkey='Kcov.res'; npart=self.resinfo['nres'];
            velf=joink(Prod.simdir,Prod.strip['velfile.res']);
            results=crd(velf,npart).Kcov(self.resinfo['masses'][0:npart])
        #store results in files
        pdb.set_trace()
        for key in results.keys():
            kkey=rootkey+'.'+key
            Prod.strip[kkey]=joink('strip',rootkey+'_'+key+'.dat')
            outf=joink(Prod.simdir,Prod.strip[kkey]);
            pt=open(outf,'w'); pt.write('#%s nparticles=%d\n'%(key,npart))
            write_from_numpy(pt,results[key],format=' %9.5f',
                             ncols=8,comment='')
        Prod.dump()
        return True
        
    def KderivTraj(self,T='270',subdir=None):
        """create trajectory of kinetic energy time-derivate"""
        from amber.amber10 import crd
        #select appropriate Prod object
        TT='T_'+T ;
        if subdir:
            Prod=self.__dict__[TT][subdir]['Prod']
        else:
            Prod=self.__dict__[TT]['Prod']
        #open file to store values
        Kdcrd=joink(Prod.simdir,'strip','Kd.crd');
        if isthere(Kdcrd): os.system('/bin/rm '+Kdcrd)
        Kdpt=open(Kdcrd,'a')
        #calculate derivative of kinetic energy
        velf=joink(Prod.simdir,Prod.strip['velfile'])
        nat=Prod.strip['nat']
        crd(velf,nat).dK_dt(self.masses[0:nat],outptf=Kdpt)
        #gzip trajectory and cleanup
        Prod.strip['Kdcrd']=joink('strip','Kd.crd')
        os.system('gzip -f '+joink(Prod.simdir,Prod.strip['Kdcrd']))
        Prod.strip['Kdcrd']+='.gz'
        Prod.dump()
        os.system('/bin/rm '+junkf)
        return True

    def KderivSnapshot(self,T='270',subdir=None,outf=None):
        """print a snapshot of Kderiv"""
        from utilities.readingWritingFiles import read_to_numpy
        from pdbfiles.smallThings import insertBfact
        if not outf: outf=junkName()+'.pdb'
        Prod=self.__dict__['T_'+T]['Prod']#select appropriate Prod object
        if subdir: Prod=self.__dict__['T_'+T][subdir]['Prod']
        Kpt=open('/projects/eflow/out/dhp/T_270/constE/Prod/strip/Kd.crd')
        buf=insertBfact(self.unsolv['preminPDB'],
                        read_to_numpy(Kpt,Prod.strip['nat']))
        open(outf,'w').write(buf)
        print 'Kderiv PDB snapshot printed to '+outf
        return False

    def KdcorrMatrix(self,T='300',subdir=None):
        """calculate the correlation matrix of Kderivative
        c_{ij)=<Kd_i(t)*Kd_j(t)>_t/(sigma_i * sigma_j)
        sigma_i=sqrt[ <(Kd_i(t)^2>_t - <Kd_i(t)>^2) ], same for sigma_j
        """
        from pdbfiles.smallThings import insertBfact
        from utilities.readingWritingFiles import read_to_numpy,\
             write_from_numpy
        #select appropriate Prod object
        TT='T_'+T ;
        if subdir:
            Prod=self.__dict__[TT][subdir]['Prod']
        else:
            Prod=self.__dict__[TT]['Prod']
        #allocate needed arrays
        nat=Prod.strip['nat']
        kd=numpy.zeros(nat); kd2=numpy.zeros(nat);
        kc=numpy.zeros(nat*nat).reshape(nat,nat)
        #open kd file and fill arrays
        junkf=junkName();
        cmd='gunzip -c '+joink(Prod.simdir,Prod.strip['Kdcrd']+' > '+junkf)
        os.system(cmd); Kpt=open(junkf)
        curr=read_to_numpy(Kpt,nat); nframe=1
        while curr.any():
            kd=kd+curr; kd2+=curr*curr; print 'nframe=',nframe
            kc+=numpy.dot(curr.reshape(nat,1),curr.reshape(1,nat))
            curr=read_to_numpy(Kpt,nat); nframe+=1
        os.system('/bin/rm '+junkf) #cleaning
        #save the average derivative
        kd=kd/nframe
        Prod.strip['Kdav']=joink('strip','Kd_average.dat')
        pt=open( joink(Prod.simdir,Prod.strip['Kdav']),'w' )
        write_from_numpy(pt,kd,format=' %7.4f',ncols=10,comment='#<dK/dt>')
        kc=kc/nframe-numpy.dot(kd.reshape(nat,1),kd.reshape(1,nat))
        #save the deviations
        kd=numpy.sqrt( kd2/nframe - kd*kd )
        Prod.strip['Kddev']=joink('strip','Kd_deviation.dat')
        pt=open( joink(Prod.simdir,Prod.strip['Kddev']),'w' )
        write_from_numpy(pt,kd,format=' %7.5f',ncols=10,
                         comment='#sqrt( <(dK/dt)^2> - (<dK/dt>)^2 )')
        Prod.strip['KddevPDB']=joink('strip','Kd_deviation.pdb')
        buf=insertBfact(self.unsolv['preminPDB'],kd)
        open(joink(Prod.simdir,Prod.strip['KddevPDB']),'w').write(buf)
        #save the correlations
        kd=numpy.dot(kd.reshape(nat,1),kd.reshape(1,nat)) #sigma_i*sigma_j
        kc=kc/kd #final correlation
        Prod.strip['Kdcorr']=joink('strip','Kd_correlation.dat')
        pt=open( joink(Prod.simdir,Prod.strip['Kdcorr']),'w' )
        write_from_numpy(pt,kc,ncols=nat,format=' %8.5f',
                         comment='#corr matrix of dK/dt')
        #update the dump file
        Prod.dump()
        return True

    def Kdcorrdecrease(self,T='300',subdir=None):
        """output one row of KdcorrMatrix into pdb as bfactor"""
        from utilities.readingWritingFiles import read_to_numpy
        from pdbfiles.smallThings import insertBfact
        #select appropriate Prod object
        TT='T_'+T ;
        if subdir:
            Prod=self.__dict__[TT][subdir]['Prod']
        else:
            Prod=self.__dict__[TT]['Prod']
        #allocate needed arrays
        nat=Prod.strip['nat']; I=1547
        pin=open(joink(Prod.simdir,Prod.strip['Kdcorr']));
        for i in range(I): pin.readline() #avoid first I lines
        kdcorr=read_to_numpy(pin,nat); kdcorr*=-1.0
        kdcorr[I-1]=0.0
        buf=insertBfact(self.unsolv['preminPDB'],kdcorr)
        open('./junk.pdb','w').write(buf)
        return False

    def diagonalizeKdcorr(self,T='300',subdir=None):
        from numpy.linalg import eigh
        from utilities.readingWritingFiles import read_to_numpy
        from pdbfiles.smallThings import insertBfact
        #select appropriate Prod object
        TT='T_'+T ;
        if subdir:
            Prod=self.__dict__[TT][subdir]['Prod']
        else:
            Prod=self.__dict__[TT]['Prod']
        #open kdcrr file and fill arrays
        nat=Prod.strip['nat']
        pt=open( joink(Prod.simdir,Prod.strip['Kdcorr']) )
        kc=read_to_numpy(pt,nat*nat,shape=(nat,nat))
        [evals,evecs]=eigh(kc); perm=numpy.argsort(evals)
        evals=evals[perm]; evecs=evecs[:,perm]
        Prod.strip['Kdcorr.evec']=joink('strip','Kd_correlation.evec.dat')
        pout=open(joink(Prod.simdir,Prod.strip['Kdcorr.evec']),'w')
        n=0
        for eval in evals:
            pout.write('# %5d %f\n'%(n+1,evals[n]))
            pout.write(' '.join([str(x) for x in evecs[n]]) )
            pout.write('\n')
            n+=1
        pout.close()
        #save as bfactors
        Kdcorrd=joink(Prod.simdir,'strip','Kdcorrd')
        os.system('/bin/mkdir -p '+Kdcorrd)
        for i in range(1,201):
            pdbout=joink(Kdcorrd,'%05d.pdb'%(i,))
            buf='REMARK  eval=%f\n'%(evals[i*-1],)
            buf+=insertBfact(self.unsolv['preminPDB'],evecs[i*-1])
            open(pdbout,'w').write(buf)
        return True
    
    def KcrossCorr(self,T='300',subdir=None):
        """calculate the time cross-correlation function"""
        #select appropriate Prod object
        TT='T_'+T ;
        if subdir:
            Prod=self.__dict__[TT][subdir]['Prod']
        else:
            Prod=self.__dict__[TT]['Prod']
        
        
    def Kflow_conductance(self,T='300'):

        """calculate the integrated cross-correlation of the
        time-derivative of the kinetic energy between two atoms"""

        #first create the trajectory of kinetic energies
        
        #read prev, curr, next frames

        #calculate the kinetic energy time-derivative as (K_next-K_prev)/2

        #calculate the 
        pass

    def Kflux(self,T='300',subdir=None):
        """calculate flux trajectory and obtain PCA components"""
        from amber.amber10 import crd
        from utilities.readingWritingFiles import write_from_numpy
        from pdbfiles.smallThings import neighborList
        #select appropriate Prod object
        TT='T_'+T ;
        if subdir:
            Prod=self.__dict__[TT][subdir]['Prod']
        else:
            Prod=self.__dict__[TT]['Prod']
        #find native contact map
        neig,dist,xyz=neighborList(open(self.unsolv['preminPDB']),co=4.0)
        #open velocity and coordinate trajectories
        nat=Prod.strip['nat']
        velf=joink(Prod.simdir,Prod.strip['velfile']);
        velobj=crd(velf,nat)
        #open flux trajectory
        Prod.strip['crdKflux']=joink('strip','kflux.crd')
        fluxf=joink(Prod.simdir,Prod.strip['crdKflux']); pt=open(fluxf,'w')
        #calculate trajectory flux of kinetic energy
        av,dev,n=velobj.kflux(self.masses[0:nat],self.VdW[0:nat],
                              neig,dist,xyz,p)
        pt.close(); os.system('gzip -f '+fluxf)
        #store net flux
        Prod.strip['Kflux.av']=joink('strip','kflux.av.crd')
        pt=open(joink(Prod.simdir,Prod.strip['Kflux.av']),'w')
        write_from_numpy(pt,av.ravel(),format=' %7.3f',ncols=10)
        #update the pickle file
        Prod.dump()
        return True

    def gneflux(self,T='300',subdir=None,mask='heavy atoms'):
        """calculate net flux in the trajectory, in the isotropic
        gaussian network approximation """
        from amber.amber10 import crd
        from networks.GN3D.gn3D import gn3D
        from utilities.readingWritingFiles import write_from_numpy
        from pdbfiles.smallThings import insertBfact
        #select appropriate Prod object
        t_p=self.__dict__['T_'+T]
        if subdir: t_p=self.__dict__['T_'+T][subdir]
        Prod=t_p['Prod']
        if mask=='all atoms':
            aa_p=t_p['aa']; nat=aa_p['nat'] #number of heavy atoms
            #load gaussian network
            gn_p=aa_p['gn']; gn=gn3D(nat)
            gn.pickleLoad( open( jk(gn_p['sd'],gn_p['obj']) ) )
            #calculate net flux
            crdobj=crd(aa_p['crd'],nat); velobj=crd(aa_p['vel'],nat)
            flux=gn.eflux(crdobj.next,velobj.next)
            #save flux to file
            aa_p['gnflux']=jk(aa_p['sd'],'gnflux.dat')
            pt=open(aa_p['gnflux'],'w')
            cmt='#flux in Gaussian network approx. for heavy atoms'
            write_from_numpy(pt,flux,format=' %9.5f',ncols=nat,comment=cmt)
            #save absorvances to file ind dat format
            absorv=flux.sum(axis=1)
            print 'total absorvance=',str(absorv.sum())
            aa_p['gnabs.dat']=jk(aa_p['sd'],'gnabs.dat')
            pt=open(aa_p['gnabs.dat'],'w')
            cmt='#absorvances in Gaussian network approx. for heavy atoms'
            write_from_numpy(pt,absorv,format=' %9.5f',ncols=1,comment=cmt)
            #save absorvances to file ind pdb format
            buf=insertBfact(self.unsolv['preminPDB'],absorv)
            ha_p['gnabs.pdb']=jk(ha_p['sd'],'gnabs.pdb')
            open(ha_p['gnabs.pdb'],'w').write(buf)
        elif mask=='heavy atoms':
            ha_p=t_p['ha']; nat=ha_p['nat'] #number of heavy atoms
            #load gaussian network
            gn_p=ha_p['gn']; gn=gn3D(nat)
            gn.pickleLoad( open( jk(gn_p['sd'],gn_p['obj']) ) )
            #calculate net flux
            crdobj=crd(ha_p['crd'],nat); velobj=crd(ha_p['vel'],nat)
            flux=gn.eflux(crdobj.next,velobj.next)
            #save flux to file
            ha_p['gnflux']=jk(ha_p['sd'],'gnflux.dat')
            pt=open(ha_p['gnflux'],'w')
            cmt='#flux in Gaussian network approx. for heavy atoms'
            write_from_numpy(pt,flux,format=' %9.5f',ncols=nat,comment=cmt)
            #save absorvances to file ind dat format
            absorv=flux.sum(axis=1)
            print 'total absorvance=',str(absorv.sum())
            ha_p['gnabs.dat']=jk(ha_p['sd'],'gnabs.dat')
            pt=open(ha_p['gnabs.dat'],'w')
            cmt='#absorvances in Gaussian network approx. for heavy atoms'
            write_from_numpy(pt,absorv,format=' %9.5f',ncols=1,comment=cmt)
            #save absorvances to file ind pdb format
            abx=numpy.zeros(Prod.strip['nat'])#init all absorvances to zero
            #change absorvances of heavy atoms to appropriate values
            for i in range(nat): abx[ self.haidx[i]-1 ]=absorv[i]
            #now insert as b-factors in the pdb
            buf=insertBfact(self.unsolv['preminPDB'],abx)
            ha_p['gnabs.pdb']=jk(ha_p['sd'],'gnabs.pdb')
            open(ha_p['gnabs.pdb'],'w').write(buf)
        return True
    
    def eflux(self,T='300',subdir=None,mask='residues'):
        """calculate net flux in the trajectory, in the quasi-harmonic
        approximation
        mask: all atoms ('atoms') or just residues ('residues')?"""
        from amber.amber10 import crd
        from utilities.readingWritingFiles import write_from_numpy
        from pdbfiles.smallThings import insertBfact
        Prod=self.__dict__['T_'+T]['Prod']#select appropriate Prod object
        if subdir: Prod=self.__dict__['T_'+T][subdir]['Prod']
        elif mask=='residues':
            nres=self.resinfo['nres'];
            crdobj=crd(joink(Prod.simdir,Prod.strip['fit2centr.res']),nres)
            velobj=crd(joink(Prod.simdir,Prod.strip['velfile.res']),nres)
            results=crdobj.eflux(velobj,masses=numpy.array(None))
            #store next flux
            key='fit2centr.res.eflux'
            Prod.strip[key]=joink('strip',key+'.dat')
            pt=open(joink(Prod.simdir,Prod.strip[key]),'w')
            comment='#net flux antisymmetrix matrix (162x162)'
            write_from_numpy(pt,results['fnet'],format=' %7.3f',
                             ncols=10,comment=comment)
            #store absorvance as .dat file
            key='fit2centr.res.absorvance'
            Prod.strip[key]=joink('strip',key+'.dat')
            pt=open(joink(Prod.simdir,Prod.strip[key]),'w')
            comment='#net absorvance'
            write_from_numpy(pt,results['absorvance'],format=' %7.3f',
                             ncols=1,comment=comment)
            #store absorvance as .pdb file            
            bfactors=numpy.zeros(Prod.strip['nat'])
            for ires in range(nres):
                bf=results['absorvance'][ires]
                for index in self.resinfo['indexes'][ires]:
                    bfactors[index]=bf
            buf=insertBfact(self.unsolv['preminPDB'],bfactors)
            Prod.strip[key]=joink('strip',key+'.pdb')
            open(joink(Prod.simdir,Prod.strip[key]),'w').write(buf)
            #save dump file
            Prod.dump()
        elif mask=='atoms': #work in progress!
            nat=Prod.strip['nat']
            crdobj=crd(joink(Prod.simdir,Prod.strip['fit2centr']),nat)
            velobj=crd(joink(Prod.simdir,Prod.strip['velfile']),nat)
            results=crdobj.eflux(velobj,masses=numpy.array(None))
            #store next flux
            key='fit2centr.eflux'
            Prod.strip[key]=joink('strip',key+'.dat')
            pt=open(joink(Prod.simdir,Prod.strip[key]),'w')
            comment='#net flux antisymmetrix matrix (%dx%d)'%(nat,nat)
            write_from_numpy(pt,results['fnet'],format=' %7.3f',
                             ncols=10,comment=comment)
            #store absorvance as .dat file
            key='fit2centr.absorvance'
            Prod.strip[key]=joink('strip',key+'.dat')
            pt=open(joink(Prod.simdir,Prod.strip[key]),'w')
            comment='#net absorvance'
            write_from_numpy(pt,results['absorvance'],format=' %7.3f',
                             ncols=1,comment=comment)
            #store absorvance as .pdb file            
            bfactors=numpy.zeros(Prod.strip['nat'])
            for iat in range(nat):
                bfactors=results['absorvance'][iat]
            buf=insertBfact(self.unsolv['preminPDB'],bfactors)
            Prod.strip[key]=joink('strip',key+'.pdb')
            open(joink(Prod.simdir,Prod.strip[key]),'w').write(buf)
            #save dump file
            Prod.dump()
            return True
    
    def meanFieldFlux(self,T='300',subdir=None,task=''):
        """calculate mean flux according to
                J_{i,j}=-w_{i,j}*k_B*(T_i-T_j)
                w_{i,j}^2=k_B*(T_i+T_j)/(2*m_i*m_j*<r_i*r_j>)
                r_i: vector displacement of atom 'i'
        """
        from utilities.readingWritingFiles import write_from_numpy
        from pdbfiles.smallThings import insertBfact
        from math import log10
        pdb.set_trace()
        T_p=self.__dict__['T_'+T]
        if subdir: T_p=self.__dict__['T_'+T][subdir]
        Prod=T_p['Prod']; nat=Prod.strip['nat']
        if 'gnmodes all atoms' in task:
            #load variance covariance matrix, find displacements
            gn_p=Prod.strip['gn']
            gn=pickle.load( open(jk(gn_p['sd'],gn_p['obj'])) )
            vc=gn.varcov
            msd=numpy.sqrt(vc.diagonal()/gn.M) #remember is mass weighted
            #load native contact map
            map=self.contmap['all']['map']
            ds=self.contmap['all']['distances']
            #filter map by cutoff: VdW_i+msd[i]+VdW_j+msd[j]
            cnts=[]
            for i in range(nat):
                ds_i=ds[i]; map_i=map[i]; cnt=[]; a=self.VdW[i]+0.5*msd[i]
                for j in range(len(ds_i)):
                    if ds_i[j]<a+self.VdW[j]+0.5*msd[j]:
                        cnt.append( map_i[j] )
                cnts.append(cnt)
            #load average kinetic energies, transform to Temperatures
            pt=open( jk(Prod.simdir,Prod.strip['Kavdat']) )
            Ts=numpy.array([float(x) for x in pt.readlines()])
            #calculate interatomic flux
            f=numpy.zeros(nat*nat).reshape(nat,nat)
            for iat in range(nat):
                cnt_i=cnts[iat]; T_i=Ts[iat]
                if self.atomSymbol[iat]=='H': T_i*=1.409#SHAKE :(
                for j in range(len(cnt_i)):
                    jat=cnt_i[j]; v=vc[iat][jat]
                    if v<=0: continue #negative force constant
                    T_j=Ts[jat]
                    if self.atomSymbol[jat]=='H': T_j*=1.409#cause of SHAKE
                    avT=(T_i+T_j)/2.0; w=numpy.sqrt(avT/v)
                    #f[i][j]: flow arriving to 'i' from 'j'
                    f[iat][jat]=w*(T_j-T_i); f[jat][iat]=-f[iat][jat]
            #output to file
            outf=jk(gn_p['sd'],'meanFieldFlux.dat'); pt=open(outf,'w')
            gn_p['meanFieldFlux']=outf;
            comment='#Mean Field Flux'
            write_from_numpy(pt,f,format=' %9.5f',ncols=8,comment=comment)
            pt.close()
            #calculate interatomic absorvances
            avs=f.sum(axis=1)
            outf=jk(gn_p['sd'],'meanFieldAbsov.dat'); pt=open(outf,'w')
            gn_p['meanFieldAbsov']=outf;
            comment='#Mean Field Absorvances'
            write_from_numpy(pt,avs,format=' %f',ncols=1,comment=comment)
            pt.close()
            outf=jk(gn_p['sd'],'meanFieldAbsov.pdb'); pt=open(outf,'w')
            gn_p['meanFieldAbsov.pdb']=outf;
            buf=insertBfact(self.unsolv['preminPDB'],(-1)*avs)
            open(outf,'w').write(buf)
            #calculate interatomic Cond
            trf=numpy.abs(f).sum(axis=1)
            outf=jk(gn_p['sd'],'meanFieldCond.dat'); pt=open(outf,'w')
            gn_p['meanFieldCond']=outf;
            comment='#Mean Field Cond'
            write_from_numpy(pt,trf,format=' %f',ncols=1,comment=comment)
            pt.close()
            outf=jk(gn_p['sd'],'meanFieldCond.pdb'); pt=open(outf,'w')
            gn_p['meanFieldCond.pdb']=outf;
            buf=insertBfact(self.unsolv['preminPDB'],avs)
            open(outf,'w').write(buf)
            #calculate interresidue flux
            nres=self.resinfo['nres']
            rf=numpy.zeros(nres*nres).reshape(nres,nres)
            for ires in range(nres-1):
                idx=self.resinfo['indexes'][ires] #atoms in ires
                for jres in range(ires+1,nres):
                    jdx=self.resinfo['indexes'][jres] #atoms in jres
                    a=0; #add atomic fluxes
                    for iat in idx:
                        for jat in jdx: a+=f[iat][jat]
                    rf[ires][jres]=a; rf[jres][ires]=-a
            outf=jk(gn_p['sd'],'meanFieldFluxRes.dat'); pt=open(outf,'w')
            gn_p['meanFieldFluxRes']=outf;
            comment='#Mean Field Flux'
            write_from_numpy(pt,rf,format=' %9.5f',ncols=8,comment=comment)
            pt.close()
            #calculate interresidue absorvances
            ravs=rf.sum(axis=1)
            outf=jk(gn_p['sd'],'meanFieldAbsovRes.dat'); pt=open(outf,'w')
            gn_p['meanFieldAbsovRes']=outf;
            comment='#Mean Field Absorvances'
            write_from_numpy(pt,ravs,format=' %f',ncols=1,comment=comment)
            pt.close()
            avs=numpy.zeros(nat)
            onezero=numpy.zeros(nat)
            for ires in range(nres):
                for iat in self.resinfo['indexes'][ires]:
                    avs[iat]=log10(abs(ravs[ires]))
            outf=jk(gn_p['sd'],'meanFieldAbsovRes.pdb'); pt=open(outf,'w')
            gn_p['meanFieldAbsovRes.pdb']=outf;
            #set DHP and NDP absorvance values to average (only
            #interested in the protein)
            #av=( avs.max()+avs.min() )/2.0
            #for ires in range(nres-3,0):
            #    for iat in self.resinfo['indexes'][ires]:
            #        avs[iat]=av
            buf=insertBfact(self.unsolv['preminPDB'],(-1)*avs,)
            open(outf,'w').write(buf)
            #calculate interresidue Cond (rf:interresidue flux)
            rtrf=numpy.abs(rf).sum(axis=1)
            outf=jk(gn_p['sd'],'meanFieldCondRes.dat');
            pt=open(outf,'w')
            gn_p['meanFieldCondRes']=outf;
            comment='#Mean Field Cond'
            write_from_numpy(pt,rtrf,format=' %f',ncols=1,comment=comment)
            pt.close()
            avs=numpy.zeros(nat)
            for ires in range(nres):
                for iat in self.resinfo['indexes'][ires]:
                    avs[iat]=log10(rtrf[ires])
            outf=jk(gn_p['sd'],'meanFieldCondRes.pdb');
            pt=open(outf,'w')
            gn_p['meanFieldCondRes.pdb']=outf;
            #set DHP and NDP Cond values to average (only
            #interested in the protein)
            #av=( avs.max()+avs.min() )/2.0
            #for ires in range(nres-3,0):
            #    for iat in self.resinfo['indexes'][ires]:
            #        avs[iat]=av
            buf=insertBfact(self.unsolv['preminPDB'],avs,)
            open(outf,'w').write(buf)
            
        return True
        
    def get_VdW(self):
        """retrieve VdW radii"""
        radii={'H':1.20, 'C':1.70, 'N':1.55, 'O':1.52, 'F':1.47, 'P':1.80,
              'S':1.80, 'Cl':1.75, 'Cu':1.4, 'Na':2.27}
        VdW=[]
        for l in open(self.preminPDB).readlines():
            if l[0:5]=='ATOM ':
                name=l[12:16]
                if name=='Na+ ':
                    VdW.append( radii['Na'] )
                else:
                    VdW.append( radii[ l[13] ] )
        self.VdW=numpy.array(VdW)
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
        
    def neededWaters(self,r=0.18):
        """print number of needed waters molecules for
        a given mass(H20)-to-mass(protein) ratio"""
        r=float(r)
        watermass=18.0 ; Nwat=int(self.mass/watermass)
        Nrwat=int(r*self.mass/watermass)
        if not 'powder' in self.__dict__.keys(): self.powder={}
        self.powder['nwat']=Nrwat
        return True

    def Bfactors(self,T='300',mask=':H*',simLength=''):
        """calculate Bfactors using trajectory fit to centroid.
        Centroid is calculated using only CA atoms"""
        pr=self.solution['T_'+T]['Prod'+simLength]
        pr.Bfactors()
        pr.BfactorsByRes()
        pr.BfactorsByMask(mask='@H*')
        return True

    def msdRatio(self,T='300',simLength=''):
        """ratio between computed and calculated m.s.d"""
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        ratio=[1/2.2]*self.Nhydr #default ratio
        Prod.strip['bfactmask']
        #pH[]
        return False
        

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

    def rmsCentr(self,T='300',simLength=''):
        from amber.Prod import Prod
        TT='T_'+T
        Prod=self.__dict__['T_'+T]['Prod'+simLength]
        Prod.rmsCentr()
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
    
def foo(*kargs,**kwargs):
    cl.cycleOverList('foo',*kargs,**kwargs)
    
def gensimdir(*kargs,**kwargs):
    """simulation directory"""
    cl.cycleOverList('gensimdir',*kargs,**kwargs)

def clusterize(*kargs,**kwargs):
    """import PDB, topology, and reference files"""
    cl.cycleOverList('clusterize',*kargs,**kwargs)
    
def importSomeFiles(*kargs,**kwargs):
    """import PDB, topology, and reference files"""
    cl.cycleOverList('importSomeFiles',*kargs,**kwargs)

def getinitPDB(*kargs,**kwargs):
    """retrieve PDB file"""
    cl.cycleOverList('getinitPDB',*kargs,**kwargs)

def storeL(*kargs,**kwargs):
    """store protein length"""
    cl.cycleOverList('storeL',*kargs,**kwargs)

def get_bounding_box(*kargs,**kwargs):
    """extract the bounding box from the solvated system"""
    cl.cycleOverList('get_bounding_box',*kargs,**kwargs)

def getAtomSymbol(*kargs,**kwargs):
    """Store atomic symbol"""
    cl.cycleOverList('getAtomSymbol',*kargs,**kwargs)
        
def neighborlist(*kargs,**kwargs):
    """find neighbor atom contact list"""
    cl.cycleOverList('neighborlist',*kargs,**kwargs)

def quasiharmonic(*kargs,**kwargs):
    """perform quasiharmonic analysis"""
    cl.cycleOverList('quasiharmonic',*kargs,**kwargs)

def qH(*kargs,**kwargs):
    """use my quasi harmonic library"""
    cl.cycleOverList('qH',*kargs,**kwargs)

def quasiharmonicAll(*kargs,**kwargs):
    """perform quasiharmonic analysis"""
    cl.cycleOverList('quasiharmonicAll',*kargs,**kwargs)

def Kflux(*kargs,**kwargs):
    """retrieve mass of the protein alone"""
    cl.cycleOverList('Kflux',*kargs,**kwargs)

def gneflux(*kargs,**kwargs):
    """calculate net flux in the trajectory, in the isotropic gaussian
    network approximation """
    cl.cycleOverList('gneflux',*kargs,**kwargs)

def eflux(*kargs,**kwargs):
    cl.cycleOverList('eflux',*kargs,**kwargs)

def meanFieldFlux(*kargs,**kwargs):
    """calculate mean flux according"""
    cl.cycleOverList('meanFieldFlux',*kargs,**kwargs)

def get_VdW(*kargs,**kwargs):
    """retrieve mass of the protein alone"""
    cl.cycleOverList('get_VdW',*kargs,**kwargs)

def get_mass(*kargs,**kwargs):
    """retrieve mass of the protein alone"""
    cl.cycleOverList('get_mass',*kargs,**kwargs)

def neededWaters(*kargs,**kwargs):
    """print number of needed waters molecules for
    a given mass(H20)-to-mass(protein) ratio"""
    cl.cycleOverList('neededWaters',*kargs,**kwargs)

def getVelCorr_sts(*kargs,**kwargs):
    """calculate velocity autocorrelation function"""
    cl.cycleOverList('getVelCorr_sts',*kargs,**kwargs)

def getDensStates(*kargs,**kwargs):
    """fourier transform of velACorr_sts"""
    cl.cycleOverList('getDensStates',*kargs,**kwargs)

def count_frames(*kargs,**kwargs):
    """import Prod structure directory"""
    cl.cycleOverList('count_frames',*kargs,**kwargs)

def importProd(*kargs,**kwargs):
    """import Prod structure directory"""
    cl.cycleOverList('importProd',*kargs,**kwargs)

def quasiharmonic(*kargs,**kwargs):
    """perform quasiharmonic analysis"""
    cl.cycleOverList('quasiharmonic',*kargs,**kwargs)

def findResIndexes(*kargs,**kwargs):
    cl.cycleOverList('findResIndexes',*kargs,**kwargs)

def keepAll(*kargs,**kwargs):
    """list of heavy atom indexes"""
    cl.cycleOverList('keepAll',*kargs,**kwargs)

def findHAIndexes(*kargs,**kwargs):
    """list of heavy atom indexes"""
    cl.cycleOverList('findHAIndexes',*kargs,**kwargs)

def keepHeavy(*kargs,**kwargs):
    cl.cycleOverList('keepHeavy',*kargs,**kwargs)
    """generate trajectory containing only the heavy atoms"""

def gnmodes(*kargs,**kwargs):
    """calculate quasiharmonic modes in the isotropic gaussian
    network approximation"""
    cl.cycleOverList('gnmodes',*kargs,**kwargs)
        
def resCM(*kargs,**kwargs):
    cl.cycleOverList('resCM',*kargs,**kwargs)

def importHotSpotInfo(*kargs,**kwargs):
    """import hot spot properties"""
    cl.cycleOverList('importHotSpotInfo',*kargs,**kwargs)

def printK(*kargs,**kwargs):
    """obtain kinetic energy plot of hot spot"""
    cl.cycleOverList('printK',*kargs,**kwargs)
    
def Ktraj(*kargs,**kwargs):
    """create trajectory of kinetic energy"""
    cl.cycleOverList('Ktraj',*kargs,**kwargs)
    
def Kcov(*kargs,**kwargs):
    cl.cycleOverList('Kcov',*kargs,**kwargs)
    
def KderivTraj(*kargs,**kwargs):
    """create trajectory of kinetic energy time-derivate"""
    cl.cycleOverList('KderivTraj',*kargs,**kwargs)
        
def KderivSnapshot(*kargs,**kwargs):
    """create trajectory of kinetic energy time-derivate"""
    cl.cycleOverList('KderivSnapshot',*kargs,**kwargs)

def KdcorrMatrix(*kargs,**kwargs):
    cl.cycleOverList('KdcorrMatrix',*kargs,**kwargs)
    
def Kdcorrdecrease(*kargs,**kwargs):
    cl.cycleOverList('Kdcorrdecrease',*kargs,**kwargs)

def diagonalizeKdcorr(*kargs,**kwargs):
    cl.cycleOverList('diagonalizeKdcorr',*kargs,**kwargs)
    
def Bfactors(*kargs,**kwargs):
    """calculate Bfactors using trajectory fit to centroid.
    Centroid is calculated using only CA atoms"""
    cl.cycleOverList('Bfactors',*kargs,**kwargs)

def get_msd(*kargs,**kwargs):
    """calculate Self-Intermediate Incoherent Scattering Function
    nt:number of time points for which to calculate the function"""
    cl.cycleOverList('get_msd',*kargs,**kwargs)

def rmsCentr(*kargs,**kwargs):
    cl.cycleOverList('rmsCentr',*kargs,**kwargs)

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
