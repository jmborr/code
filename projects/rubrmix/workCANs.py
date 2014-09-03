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

my_project='rubrmix'                      #substitute 'my_project'
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
    #fwhm={'original':32,'BASIS':3}
    z=1.0/2.354820045
    sigma={'original':z*32,'BASIS':z*3}
    ss=sigma[type] ; ss=ss*ss
    return exp(-E*E/(2*ss))/sqrt(2*pi*ss)

def convT(t,type='BASIS',ss=None):
    """fourier transform of distrE
    Time units are in picoseconds"""
    from math import exp,pi,sqrt
    sigma={'original':49.0,'BASIS':521.0} #time resolution windows
    if not ss: ss=sigma[type]
    ss=2*ss*ss ; f=1/sqrt(2*pi*ss)
    return f*exp(-t*t/ss)
            
    
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

    def insertProp(self,keys,prop):
        """Example: -b "dir dir2 key prop" """
        if type(keys)==type(''): keys=[keys,]
        pt=self.__dict__
        try:
            if len(keys)>1:
                for key in keys[0:-1]: pt=pt[key]
            pt[keys[-1]]=prop
        except:
            sys.stderr.write('Could not insert property\n')
            return False
        return True
        
    def conditions2dict(self,conditions):
        """NOTE: order is relevant if directory must be created.
        Create appropriate directory for a given set of simulation
        conditions.
        Update object, every condition spawn a dictionary
        Return dictionary handle"""
        cds=conditions.split()
        finald=self.simdir+'/'+'/'.join(cds)
        os.system('/bin/mkdir -p '+finald)
        pt=self.__dict__; i=0
        try:
            for cd in cds:
                pt=pt[cd]; i+=1
        except: #extend the dictionary
            for j in range(i,len(cds)):
                pt[cds[j]]={}; pt=pt[cds[j]]
        pt['simdir']=finald
        return pt
    
    def gensimdir(self):
        """simulation directory"""
        self.simdir=joink(pjdout,self.id)
        if not isthere(self.simdir):
            os.system('/bin/mkdir -p '+self.simdir)
        return True

    def insertTopRef(self,conditions):
        """insert topology and reference files. Some naming convengtion"""
        pt=self.conditions2dict(conditions)
        pt['top']=jk(pt['simdir'],self.id+'.top')
        pt['ref']=jk(pt['simdir'],self.id+'.crd')
        return True
    
    def importProd(self,conditions='solvated T=50'):
        """import Prod structure directory"""
        from amber.Prod import Prod
        pt=self.conditions2dict(conditions)
        prodd=jk(pt['simdir'],'Prod')
        pr=Prod(self.id,prodd,module_name='amber.Prod')
        #import topology and reference file to Prod object
        pr.importStripTopCrd(self.dry['top'],self.dry['ref'])
        pr.count_frames(mode='unsolvated') #check number of frames
        #several RMS calculations to check consistency of simulation
        #pr.getAvg()       #averages, already done in Jaguar
        pr.rmsJump(outcrd=True)#consec rms and output fitted trajectory
        pr.rmsRef()        #rms with respect to reference
        pr.centr_and_avg() #find centroid
        pr.rmsCentr()      #rms with respect to centroid
        pr.rmsFirst()      #rms with respect to first frame
        pr.fit2centr()     #fit trajectory to centroid
        pr.dump()
        pt['Prod']=pr
        return True

    def keepH(self,conditions='solvated T=50'):
        """generate trajectory containing only the hydrogens"""
        print 'generate trajectory containing only the hydrogens'
        from amber.amber10 import top,crd,ptraj
        pt=self.conditions2dict(conditions)
        pt['onlyH']={'simdir':joink(pt['simdir'],'onlyH')}
        pH=pt['onlyH']
        if not isthere(pH['simdir']): os.system('mkdir -p '+pH['simdir'])
        #check number of hydrogens read
        Nhydr=top(self.dry['top']).numberAtoms(mask='@H*')
        pH['Nhydr']=Nhydr #store also here
        if 'Nhydr' in self.__dict__ and self.Nhydr!=Nhydr:
            buf='different number of hydrogens than stored!\n%d vs %d'%\
                 (self.Nhydr,Nhydr)
            Bye(buf)
        else:
            self.Nhydr=Nhydr
        Prod=pt['Prod'] #Production object
        for keyword in ('crdfile','fit2prev','fit2centr'):
            print 'processing '+keyword
            crdf=joink(Prod.simdir,Prod.strip[keyword])
            outf=basename(crdf.replace('.crd.gz','.H.crd'))
            outf=joink(pH['simdir'],outf)
            ptr=ptraj(self.dry['top'],crdf)
            ptr.strip(entities=['!@H*',]) #keep only hydrogens
            ptr.go(outcrd=outf,outfile='./junk.log')
            print 'gzip H-only trajectory'
            os.system('gzip -f '+outf)
            pH['crd.'+keyword]=outf+'.gz'
        return True

    def getEISF(self,conditions='solvated T_50',task='eisf2'):
        """calculate Self-Intermediate Incoherent Scattering Function
        nt:number of time points for which to calculate the function"""
        from amber.amber10 import crd
        pt=self.conditions2dict(conditions)
        qq=[]; iq=0.1
        while iq<2.6: qq.append(iq); iq+=0.1
        nq=len(qq)
        if task=='eisf2':
            from scatter.scatter import EISF2
            self.Nhydr=375
            os.chdir( jk(pt['simdir'],'onlyH') )
            crdobj=crd('fit2avg.H.crd.gz',self.Nhydr)
            skipframes=0; crdobj.loadframes(skip=skipframes)
            scattering=EISF2(qq,crdobj.frames)
            buf=''
            for (q,s) in scattering: buf+=str(q)+' '+str(s)+'\n'
            open('eisf2.dat','w').write(buf)
            crdobj.close()
            return False
        elif task=='eisf fit2avg.H.crd':
            from scatter.scatter import EISF
            self.Nhydr=375
            os.chdir( jk(pt['simdir'],'onlyH') )
            crdobj=crd('fit2avg.H.crd.gz',self.Nhydr)
            skipframes=0; crdobj.loadframes(skip=skipframes)
            scattering=EISF(qq,crdobj.frames)
            buf=''
            for (q,s) in scattering: buf+=str(q)+' '+str(s)+'\n'
            open('eisf.fit2avg.H.crd.dat','w').write(buf)
            crdobj.close()
            return True
        elif task=='shuffle eisf fit2avg.H.crd':
            TT=[50,75,100,125,150,175,200,225,250,275,300,325,350]
            nT=len(TT)
            scattering=numpy.zeros(nT*nq).reshape(nT,nq)
            iT=0
            for T in [str(x) for x in TT]:
                inf=jk(pt['T_'+T]['simdir'],'onlyH/eisf.fit2avg.H.crd.dat')
                iq=0
                for l in open(inf).readlines():
                    scattering[iT][iq]=float(l.split()[-1]); iq+=1
                iT+=1
            pt['eisf.fit2avg.H.crd']=jk(pt['simdir'],'eisf.fit2avg.H.crd.dat')
            buf='# T '+' '.join([str(x) for x in qq])+'\n'
            for iT in range(nT):
               buf+=' %3d'%TT[iT]
               for iq in range(nq): buf+=' %6.4f'%scattering[iT][iq]
               buf+='\n'
            open(pt['eisf.fit2avg.H.crd'],'w').write(buf)
            return True                
        elif task=='eisf H.crd':
            from scatter.scatter import EISF
            self.Nhydr=375
            os.chdir( jk(pt['simdir'],'onlyH') )
            crdobj=crd('H.crd.gz',self.Nhydr)
            skipframes=0; crdobj.loadframes(skip=skipframes)
            scattering=EISF(qq,crdobj.frames)
            buf=''
            for (q,s) in scattering: buf+=str(q)+' '+str(s)+'\n'
            open('eisf.H.crd.dat','w').write(buf)
            crdobj.close()
            return True
        elif task=='shuffle eisf H.crd':
            TT=[50,75,100,125,150,175,200,225,250,275,300,325,350]
            nT=len(TT)
            scattering=numpy.zeros(nT*nq).reshape(nT,nq)
            iT=0
            for T in [str(x) for x in TT]:
                inf=jk(pt['T_'+T]['simdir'],'onlyH/eisf.H.crd.dat')
                iq=0
                for l in open(inf).readlines():
                    scattering[iT][iq]=float(l.split()[-1]); iq+=1
                iT+=1
            pt['eisf.H.crd']=jk(pt['simdir'],'eisf.H.crd.dat')
            buf='# T '+' '.join([str(x) for x in qq])+'\n'
            for iT in range(nT):
               buf+=' %3d'%TT[iT]
               for iq in range(nq): buf+=' %6.4f'%scattering[iT][iq]
               buf+='\n'
            open(pt['eisf.H.crd'],'w').write(buf)
            return True                
        elif task=='eisf center.H.crd':
            from scatter.scatter import EISF
            self.Nhydr=375
            os.chdir( jk(pt['simdir'],'onlyH') )
            crdobj=crd('center.H.crd.gz',self.Nhydr)
            skipframes=0; crdobj.loadframes(skip=skipframes)
            scattering=EISF(qq,crdobj.frames)
            buf=''
            for (q,s) in scattering: buf+=str(q)+' '+str(s)+'\n'
            open('eisf.center.H.crd.dat','w').write(buf)
            crdobj.close()
            return True
        elif task=='shuffle eisf center.H.crd':
            TT=[50,75,100,125,150,175,200,225,250,275,300,325,350]
            nT=len(TT)
            scattering=numpy.zeros(nT*nq).reshape(nT,nq)
            iT=0
            for T in [str(x) for x in TT]:
                inf=jk(pt['T_'+T]['simdir'],'onlyH/eisf.center.H.crd.dat')
                iq=0
                for l in open(inf).readlines():
                    scattering[iT][iq]=float(l.split()[-1]); iq+=1
                iT+=1
            pt['eisf.center.H.crd']=jk(pt['simdir'],'eisf.center.H.crd.dat')
            buf='# T '+' '.join([str(x) for x in qq])+'\n'
            for iT in range(nT):
               buf+=' %3d'%TT[iT]
               for iq in range(nq): buf+=' %6.4f'%scattering[iT][iq]
               buf+='\n'
            open(pt['eisf.center.H.crd'],'w').write(buf)
            return True                
        elif task=='shuffle':
            TT=[50,75,100,125,150,175,200,225,250,275,300,325,350]
            nT=len(TT)
            scattering=numpy.zeros(nT*nq).reshape(nT,nq)
            iT=0
            for T in [str(x) for x in TT]:
                inf=pt['T_'+T]['onlyH']['eisf2']
                iq=0
                for l in open(inf).readlines():
                    scattering[iT][iq]=float(l.split()[-1]); iq+=1
                iT+=1
            pt['eisf2']=jk(pt['simdir'],'eisf2.dat')
            pout=open(pt['eisf2'],'w')
            pout.write('# T '+' '.join([str(x) for x in qq])+'\n')
            for iT in range(nT):
               buf=' %3d'%TT[iT]
               for iq in range(nq): buf+=' %6.4f'%scattering[iT][iq]
               pout.write(buf+'\n')
        elif task=='eisf shell1.H.crd':
            from scatter.scatter import EISF
            #number of hydrogens in the shell1.crd trajectory
            nH={'ref15-shel':\
                {'50':416,  '75':416, '100':416, '125':416,
                '150':416, '175':416, '200':416, '225':416,
                '250':416, '275':416, '300':416, '325':416, '350':416},
                'ref15-shel_H':\
                {'50':418,  '75':418, '100':418, '125':418,
                '150':418, '175':418, '200':418, '225':416,
                '250':418, '275':418, '300':416, '325':416, '350':418},
                'ref15-shel_D':\
                {'50':416,  '75':414, '100':416, '125':416,
                '150':416, '175':416, '200':414, '225':416,
                '250':416, '275':416, '300':416, '325':416, '350':412}
                }
            Temp=re.compile('T_(\d+)').search(conditions).group(1)
            n=nH[self.id][Temp]
            os.chdir( jk(pt['simdir'],'onlyH') )
            crdobj=crd('shell1.H.crd.gz',n,cformat='8.3f')
            skipframes=0; crdobj.loadframes(skip=skipframes)
            scattering=EISF(qq,crdobj.frames)
            buf=''
            for (q,s) in scattering: buf+=str(q)+' '+str(s)+'\n'
            open('eisf.shell1.H.crd.dat','w').write(buf)
            crdobj.close()
            return True
        elif task=='shuffle eisf shell1.H.crd':
            TT=[50,75,100,125,150,175,200,225,250,275,300,325,350]
            nT=len(TT)
            scattering=numpy.zeros(nT*nq).reshape(nT,nq)
            iT=0
            for T in [str(x) for x in TT]:
                inf=jk(pt['T_'+T]['simdir'],'onlyH/eisf.shell1.H.crd.dat')
                iq=0
                for l in open(inf).readlines():
                    scattering[iT][iq]=float(l.split()[-1]); iq+=1
                iT+=1
            pt['eisf.shell1.H.crd']=jk(pt['simdir'],'eisf.shell1.H.crd.dat')
            buf='# T '+' '.join([str(x) for x in qq])+'\n'
            for iT in range(nT):
               buf+=' %3d'%TT[iT]
               for iq in range(nq): buf+=' %6.4f'%scattering[iT][iq]
               buf+='\n'
            open(pt['eisf.shell1.H.crd'],'w').write(buf)
            return True                
               
        
    def getSIIScF_sts(self,conditions='solvated T_50',tl=80,nt=400,
                      task=''):
        """calculate Self-Intermediate Incoherent Scattering Function
        tl: total time desired, in picoseconds
        nt: number of time points to sample"""
        from scatter.scatter import SIIScF,storeSIIScF,fourier
        from amber.amber10 import crd
        tl=float(tl)
        pt=self.conditions2dict(conditions)
        Prod=pt['Prod']; pH=pt['onlyH']
        skipframes=int(max(0,tl/(nt*Prod.times['twx'])-1))
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']; dw=1.0/(nt*dt)
        tt=[]; ww=[]; 
        for i in range(nt): tt.append(i*dt); ww.append(i*dw)        
        qq=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] #Angstroms^(-1)
        nq=len(qq)
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,pH['Nhydr'])
        crdobj.loadframes(skip=skipframes)
        scattering=SIIScF(qq,crdobj.frames,dt,tl,nt)
        print '#store all intermediate factors'
        pH['isf']=joink(pH['simdir'],'isf')
        os.system('mkdir -p '+pH['isf'])
        for i in range(self.Nhydr):
            outfile=joink(pH['isf'],'%04d.dat'%(i+1,))
            storeSIIScF(outfile,scattering[i],qq,tt)
        print '#store the average intermediate structure factor'
        pH['SIIScF_sts']=crdf.replace('.crd.gz','.SIIScF_sts.dat')
        avscattering=scattering.sum(axis=0)/len(scattering) 
        storeSIIScF(pH['SIIScF_sts'],avscattering,qq,tt)
        crdobj=None #deallocate memory        
        print '#perform Fourier transform'
        gaussianenvelope=[]; factor=1.0/convT(tt[0])
        """We choose not to convolute with the resolution function
        because it produces negative values in the structure factor at
        low-frquencies, no doubt due to insufficient sampling in the
        long time-scale"""
        for tti in tt:
            #gaussianenvelope.append(factor*convT(tti))
            gaussianenvelope.append(1.0)
        for iatom in range(self.Nhydr):
            sc=scattering[iatom]
            for iq in range(nq):
                sc[iq]=fourier(sc[iq],conv=gaussianenvelope).real
        print '#store all structure factors'
        pH['dsf']=joink(pH['simdir'],'dsf')
        os.system('mkdir -p '+pH['dsf'])
        for i in range(self.Nhydr):
            outfile=joink(pH['dsf'],'%04d.dat'%(i+1,))
            storeSIIScF(outfile,scattering[i],qq,ww)
        print '#store the average structure factor'
        pH['Sq_sts']=crdf.replace('.crd.gz','.Sq_sts.dat')
        for iq in range(nq):
            avscattering[iq]=fourier(avscattering[iq]*numpy.array(gaussianenvelope)).real
            #avscattering[iq]=fourier(avscattering[iq],
            #                         conv=gaussianenvelope).real
        storeSIIScF(pH['Sq_sts'],avscattering,qq,ww)
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

def insertProp(*kargs,**kwargs):
    """insert some property"""
    cl.cycleOverList('insertProp',*kargs,**kwargs)
    
def gensimdir(*kargs,**kwargs):
    """simulation directory"""
    cl.cycleOverList('gensimdir',*kargs,**kwargs)
    
def insertTopRef(*kargs,**kwargs):
    cl.cycleOverList('insertTopRef',*kargs,**kwargs)
    
def importProd(*kargs,**kwargs):
    cl.cycleOverList('importProd',*kargs,**kwargs)
    
def keepH(*kargs,**kwargs):
    cl.cycleOverList('keepH',*kargs,**kwargs)
    
def getSIIScF_sts(*kargs,**kwargs):
    cl.cycleOverList('getSIIScF_sts',*kargs,**kwargs)
    
def getEISF(*kargs,**kwargs):
    cl.cycleOverList('getEISF',*kargs,**kwargs)

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
