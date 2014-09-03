#!/usr/bin/python

import os,sys,re,numpy,pickle,pdb
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir,projectsdir
from utilities.small_utilities import chomp,Bye,junkName
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
joink=os.path.join; jk=joink
isthere=os.path.exists
basename=os.path.basename

my_project='dhfr_evb'                     #substitute 'my_project'
pcd=joink(codedir,'projects',my_project)  #code directory
pjd=joink(projectsdir,my_project)         #projects directory
pjdout=joink(projectsdir,my_project,'out')
scd=joink(scratchdir,my_project)          #scratach directory
scd=joink(scratchdir,my_project,'out')
collf=joink(pjd,'varsets.can')            #you may substitute scd for pjd
canlistf=joink(pjd,'out','list')          #you may substitute scd for pjd
repod=joink(pjd,'repod')                  #you may substitute scd for pjd
    
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

    def importTopRef(self):
        """import teopology file"""
        if 'reactant' not in self.__dict__.keys():self.reactant={}
        if 'product' not in self.__dict__.keys():self.product={}

        self.reactant['dir']=jk(self.simdir,'reactant')
        self.product['dir']=jk(self.simdir,'product')

        self.reactant['top']=jk(self.reactant['dir'],'top')
        self.product['top']=jk(self.product['dir'],'top')

        self.reactant['unsolvated']={}
        self.product['unsolvated']={}

        dir=jk(self.reactant['dir'],'unsolvated')
        self.reactant['unsolvated']['top']=jk(dir,'top')
        self.reactant['unsolvated']['ref']=jk(dir,'ref')
        self.reactant['unsolvated']['pdb']=jk(dir,'pdb')
        return True
        
    def importProd(self,T='100',L='0.00'):
        """import Prod structure directory"""
        from amber.Prod import Prod
        #prepare object structure
        if 'temp' not in self.__dict__: self.temp={}
        if T not in self.temp.keys(): self.temp[T]={}
        if 'lambda' not in self.temp[T].keys(): self.temp[T]['lambda']={}
        plambda=self.temp[T]['lambda']
        if L not in plambda.keys(): plambda[L]={}
        pL=plambda[L]
        prodd=jk(self.simdir,'T_%s/L-%s/Prod'%(T,L))
        pr=Prod(self.id,prodd,module_name='amber.Prod')
        pr.strip['maskl']=[':WAT',':Na+']
        #import topology and reference file to Prod object
        pr.importStripTopCrd(self.reactant['unsolvated']['top'],
                             self.reactant['unsolvated']['ref'])
        pr.count_frames(mode='unsolvated') #check number of frames
        #several RMS calculations to check consistency of simulation
        #pr.getAvg()          #averages
        pr.rmsJump(outcrd=True)#consec rms and output fitted trajectory
        pr.rmsRef()        #rms with respect to reference
        pr.centr_and_avg() #find centroid
        pr.rmsCentr()      #rms with respect to centroid
        pr.rmsFirst()      #rms with respect to first frame
        pr.fit2centr()     #fit trajectory to centroid
        pL['Prod']=pr
        return True

    def importMasses(self):
        """import atomic masses of system"""       
        import numpy
        m=open(jk(self.simdir,'mass.dat'),'r').read()
        m=[ float(x) for x in m.split() ]
        self.masses=numpy.array(m)
        return True
    
    def get_VdW(self):
        """retrieve VdW radii"""
        radii={'H':1.20, 'C':1.70, 'N':1.55, 'O':1.52, 'F':1.47, 'P':1.80,
              'S':1.80, 'Cl':1.75, 'Cu':1.4, 'Na':2.27}
        VdW=[]
        for l in open(self.reactant['unsolvated']['pdb']).readlines():
            if l[0:5]=='ATOM ':
                name=l[12:16]
                if name=='Na+ ': VdW.append( radii['Na'] )
                elif name=='Cl- ': VdW.append( radii['Cl'] )
                else: VdW.append( radii[ l[13] ] )
        self.VdW=numpy.array(VdW)
        return True

    def importHotSpotInfo(self,massf=None,indexf=None):
        """import hot spot properties"""
        import numpy
        if 'hot_spot' not in self.__dict__.keys(): self.hot_spot={}
        print '#import PDB indexes for atoms'
        if indexf:
            self.hot_spot['indexes']=[]
            index=0 #indexes start with 0, instead of 1
            for line in open(indexf).readlines():
                if '1' in line: self.hot_spot['indexes'].append(index)
                index+=1
        print '#import masses for hot_spot'
        M=[]
        for index in self.hot_spot['indexes']:
            M.append(self.masses[index])
        self.hot_spot['masses']=numpy.array(M)
        return True


    def findResIndexes(self):

        """for each residue, find the list of atom indexes

        NOTE: first atom index is 0, NOT 1

        NOTE: we create and extra residue with the 24 catalytic site
        atoms right before the Na+ ions. We remove the 24 catalytic
        site atoms from the DHP and NDP residues. This is residue with
        index 161.

        NOTE: after residue 161, we continue until all ions and water
        residues are accounted for
        """
        self.reactant['solvated']={}
        self.reactant['solvated']['pdb']=self.simdir+'/reactant/solvated/'+self.id+'.premin.pdb'
        nres=162 #protein residues + DHP + NPD + energy-sink
        pdbf=self.reactant['solvated']['pdb']
        ptin=open(pdbf,'r'); l=ptin.readline()
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
        pdb.set_trace()
        return True

    def gnmodes(self,T='100',L='0.00',task='init all atoms'):
        """calculate quasiharmonic modes in the isotropic gaussian
        network approximation"""
        from networks.GN3D.gn3D import gn3D
        from amber.Prod import Prod
        from amber.amber10 import crd
        pL=self.temp[T]['lambda'][L];
        Prod=pL['Prod']
        if 'all atoms' in task: #modes for heavy atoms
            aa_p=Prod.strip; nat=aa_p['nat']; gn=gn3D(nat); gn_p=None
            sd=jk(Prod.simdir,'strip','gn');
            if 'init' in task:    #initialize modes directory
                os.system('mkdir -p '+sd);
                aa_p['gn']={'sd':sd}; gn_p=aa_p['gn']
                aa_p['crd']=jk(Prod.simdir,aa_p['fit2centr'])
                aa_p['masses']=self.masses[0:nat]
            else:                 #load initialized modes info
                gn_p=aa_p['gn']
                gn.pickleLoad( open( jk(sd,gn_p['obj']) ) )
                if 'findVarCovar' in task:
                    crdf=aa_p['crd']
                    crdobj=crd(crdf,nat)#coordinate trajectory object
                    #pdb.set_trace()                    
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
        #pdb.set_trace()        
        gn.pickleDump( open( jk(sd,gn_p['obj']),'w' ) )
        return True

    def neighborlist(self):
        """find neighbor atom contact list"""
        from pdbfiles.smallThings import neighborList
        #find native contact map
        pdbf=self.reactant['unsolvated']['pdb']
        neig,dist,xyz=neighborList(open(pdbf),co=6.0)
        self.contmap={'all':{'co':6.0,'map':neig,'distances':dist}}
        return True

    def Ktraj(self,T='100',L='0.00'):
        """create trajectory of kinetic energy"""
        from amber.amber10 import crd
        from pdbfiles.smallThings import insertBfact
        from utilities.readingWritingFiles import write_from_numpy
        pL=self.temp[T]['lambda'][L];
        Prod=pL['Prod']
        Kcrd=joink(Prod.simdir,'strip','K.crd');
        if isthere(Kcrd): os.system('/bin/rm '+Kcrd)
        Kpt=open(Kcrd,'a')
        velf=joink(Prod.simdir,Prod.strip['velfile'])
        nat=Prod.strip['nat']
        crdobj=crd(velf,nat)
        Kav,Kdev=crdobj.kinetic_energy(self.masses[0:nat],outptf=Kpt)
        #gzip trajectory
        os.system('gzip -f '+Kcrd)
        Prod.strip['Kcrd']=joink('strip','K.crd.gz')
        #save average Kinetic energy as .pdb
        buf=insertBfact(self.reactant['unsolvated']['pdb'],Kav)
        Prod.strip['Kavpdb']=joink('strip','Kav.pdb')
        open(jk(Prod.simdir,Prod.strip['Kavpdb']),'w').write(buf)
        #save average Kinetic energy as .dat
        buf='\n'.join([str(x) for x in Kav])
        Prod.strip['Kavdat']=joink('strip','Kav.dat')
        open(jk(Prod.simdir,Prod.strip['Kavdat']),'w').write(buf)
        Prod.dump()
        return True

    def meanFieldFlux(self,T='100',L='0.00',task=''):
        """calculate mean flux according to
                J_{i,j}=-w_{i,j}*k_B*(T_i-T_j)
                w_{i,j}^2=k_B*(T_i+T_j)/(2*m_i*m_j*<r_i*r_j>)
                r_i: vector displacement of atom 'i'
        """
        from utilities.readingWritingFiles import write_from_numpy
        from pdbfiles.smallThings import insertBfact
        from math import log10
        pL=self.temp[T]['lambda'][L];Prod=pL['Prod'];nat=Prod.strip['nat']
        if 'gnmodes all atoms' in task:
            print '#load variance covariance matrix, find displacements'
            gn_p=Prod.strip['gn']
            gn=pickle.load( open(jk(gn_p['sd'],gn_p['obj'])) )
            vc=gn.varcov
            msd=numpy.sqrt(vc.diagonal()/gn.M) #remember is mass weighted
            print '#load native contact map'
            map=self.contmap['all']['map']
            ds=self.contmap['all']['distances']
            print '#filter map by cutoff: VdW_i+msd[i]+VdW_j+msd[j]'
            cnts=[]
            for i in range(nat):
                ds_i=ds[i]; map_i=map[i]; cnt=[]; a=self.VdW[i]+0.5*msd[i]
                for j in range(len(ds_i)):
                    if ds_i[j]<a+self.VdW[j]+0.5*msd[j]:
                        cnt.append( map_i[j] )
                cnts.append(cnt)
            print '#load av. kinetic energies, transform to Temperatures'
            pt=open( jk(Prod.simdir,Prod.strip['Kavdat']) )
            Ts=numpy.array([float(x) for x in pt.readlines()])
            print '#calculate interatomic flux'
            f=numpy.zeros(nat*nat).reshape(nat,nat)
            for iat in range(nat):
                cnt_i=cnts[iat]; T_i=Ts[iat]
                for j in range(len(cnt_i)):
                    jat=cnt_i[j]; v=vc[iat][jat]
                    if v<=0: continue #negative force constant
                    T_j=Ts[jat]
                    avT=(T_i+T_j)/2.0; w=numpy.sqrt(avT/v)
                    #f[i][j]: flow arriving to 'i' from 'j'
                    f[iat][jat]=w*(T_j-T_i); f[jat][iat]=-f[iat][jat]
            #output to file
            outf=jk(gn_p['sd'],'meanFieldFlux.dat'); pt=open(outf,'w')
            gn_p['meanFieldFlux']=outf;
            comment='#Mean Field Flux'
            write_from_numpy(pt,f,format=' %9.5f',ncols=8,comment=comment)
            pt.close()
            print '#calculate interatomic absorvances'
            avs=f.sum(axis=1)
            outf=jk(gn_p['sd'],'meanFieldAbsov.dat'); pt=open(outf,'w')
            gn_p['meanFieldAbsov']=outf;
            comment='#Mean Field Absorvances'
            write_from_numpy(pt,avs,format=' %f',ncols=1,comment=comment)
            pt.close()
            outf=jk(gn_p['sd'],'meanFieldAbsov.pdb'); pt=open(outf,'w')
            gn_p['meanFieldAbsov.pdb']=outf;
            pdbf=self.reactant['unsolvated']['pdb']
            buf=insertBfact(pdbf,(-1)*avs)
            open(outf,'w').write(buf)
            print '#calculate interatomic Conductance'
            trf=numpy.abs(f).sum(axis=1)
            outf=jk(gn_p['sd'],'meanFieldCond.dat'); pt=open(outf,'w')
            gn_p['meanFieldCond']=outf;
            comment='#Mean Field Cond'
            write_from_numpy(pt,trf,format=' %f',ncols=1,comment=comment)
            pt.close()
            outf=jk(gn_p['sd'],'meanFieldCond.pdb'); pt=open(outf,'w')
            gn_p['meanFieldCond.pdb']=outf;
            buf=insertBfact(pdbf,avs)
            open(outf,'w').write(buf)
            print '#calculate interresidue flux'
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
            print '#calculate interresidue absorvances'
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
            buf=insertBfact(pdbf,(-1)*avs,)
            open(outf,'w').write(buf)
            print '#calculate interresidue Cond (rf:interresidue flux)'
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
            buf=insertBfact(pdbf,avs,)
            open(outf,'w').write(buf)            
        return False#True

    def prop_vs_L(self,T='100',task=''):
        """calculate evolution of some properties versus lambda value"""

        LL=['0.00','0.10','0.20','0.30','0.40',
           '0.50','0.60','0.70','0.80','0.90','1.00']
        nL=len(LL); nres=self.resinfo['nres']
        pL=self.temp[T]['lambda']
        if task=='meanFieldCondRes':            
            conds=numpy.zeros(nres*nL).reshape(nL,nres); iL=0;
            for L in LL:
                condf=pL[L]['Prod'].simdir+'/strip/gn/meanFieldCondRes.dat'
                conds[iL]=[ float(x) for x in open(condf).readlines()[1:] ]
                iL+=1
            conds=conds.transpose()
            outdir=self.simdir+'/T_'+T+'/analysis/gn/meanFieldCondRes'
            os.system('/bin/mkdir -p '+outdir)
            from utilities.readingWritingFiles import write_from_numpy
            for ires in range(nres):
                ptout=open('%s/%03d.dat'%(outdir,ires+1),'w')
                cm='#meanFieldCondRes vs lambda for resnumber %03d'%(ires+1)
                write_from_numpy(ptout,conds[ires],ncols=1,comment=cm)
        return False

        
            

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

def conductance(inp=None,inp2=None,inp3=None,out=None,out2=None,
                task=None):
    if task=='deviation to error':
        """calculate the ratio of the change in conductance from EVB0
        to EVB4 to the estimation in the error of the conductance for
        EVB0"""
        from utilities.readingWritingFiles import read_column,write_from_numpy

        #read the conductances for the chunk trajectories of EVB0
        chunks=('1_14000', '7001_21000', '14001_28000',
                '21001_35000', '28001_42000' )
        N = None
        for chunk in chunks:
            file = 'EVB0/Prod/conductance/meanFieldCondRes_%s.dat'%chunk
            conds = numpy.array( read_column(file,1,isFloat=1) )
            if not N:
                vals = conds
                vals = vals.reshape( len(vals), 1) #make it a 2D array
                N = len(conds)
            else:
                conds = conds.reshape(N,1)
                vals = numpy.append( vals,conds, axis=1)
        #read the conductances for EVB0
        file = 'EVB0/Prod/conductance/meanFieldCondRes.dat'
        conds0 = numpy.array( read_column(file,1,isFloat=1) )
        #read the conductances for EVB4
        file = 'EVB4/Eflow/unsolv/conductance/meanFieldCondRes.dat'
        conds4 = numpy.array( read_column(file,1,isFloat=1) )
        # estimate error in EVB0 conductaces as max(vals) - min(vals)
        # print ratio of change in conductance to half the estimated
        x = conds0-conds4
        y = numpy.max(vals,axis=1)-numpy.min(vals,axis=1)
        out = open('%s/analysis/deviation_to_error.dat'%pjd,'w')
        write_from_numpy(out,x/y,format='%6.3f',ncols=1,
                         comment='#ratio (EVB4-EVB9) to error')
        return False
    elif task=='dif meanFieldFluxRes':
        from utilities.readingWritingFiles import read_to_numpy,write_from_numpy
        file='%s/out/EVB0/Prod/conductance/meanFieldFluxRes.dat'%pjd
        flux0=read_to_numpy(open(file),1000000)

        file='%s/out/EVB4/Eflow/unsolv/conductance/meanFieldFluxRes.dat'%pjd
        flux4=read_to_numpy(open(file),1000000)
        flux0 = numpy.where(flux0>0,flux0,0.0) #discard negative flux
        flux4 = numpy.where(flux4>0,flux4,0.0) #discard negative flux
        dflux = flux4-flux0
        file='%s/analysis/dif_meanFieldFluxRes.dat'%pjd
        write_from_numpy(open(file,'w'),dflux,format=' %9.5f',ncols=8,
                         comment='# flux(EVB4) - flux(EVB0)')
        return False
        
def foo2(*kargs,**kwargs):
    """create seq.dat files
    Example: ... -c id=a0512,switch=hello"""
    cl.cycleOverList('foo',*kargs,**kwargs)
    
def gensimdir(*kargs,**kwargs):
    """simulation directory"""
    cl.cycleOverList('gensimdir',*kargs,**kwargs)

def importTopRef(*kargs,**kwargs):
    """import topology structure directory"""
    cl.cycleOverList('importTopRef',*kargs,**kwargs)

def importProd(*kargs,**kwargs):
    """import Prod structure directory"""
    cl.cycleOverList('importProd',*kargs,**kwargs)

def importMasses(*kargs,**kwargs):
    cl.cycleOverList('importMasses',*kargs,**kwargs)

def get_VdW(*kargs,**kwargs):
    cl.cycleOverList('get_VdW',*kargs,**kwargs)

def importHotSpotInfo(*kargs,**kwargs):
    cl.cycleOverList('importHotSpotInfo',*kargs,**kwargs)

def findResIndexes(*kargs,**kwargs):
    cl.cycleOverList('findResIndexes',*kargs,**kwargs)

def Ktraj(*kargs,**kwargs):
    cl.cycleOverList('Ktraj',*kargs,**kwargs)

def meanFieldFlux(*kargs,**kwargs):
    cl.cycleOverList('meanFieldFlux',*kargs,**kwargs)

def prop_vs_L(*kargs,**kwargs):
    cl.cycleOverList('prop_vs_L',*kargs,**kwargs)

def gnmodes(*kargs,**kwargs):
    cl.cycleOverList('gnmodes',*kargs,**kwargs)

def neighborlist(*kargs,**kwargs):
    cl.cycleOverList('neighborlist',*kargs,**kwargs)

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
