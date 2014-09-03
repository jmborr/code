#!/usr/bin/python

import os,sys,re,numpy,pdb
from tempfile import mkstemp
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir,projectsdir
from utilities.small_utilities import chomp,Bye
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
pwrite=sys.stdout.write
joink=os.path.join; jk=joink
isthere=os.path.exists
basename=os.path.basename
 
my_project='meso_thermo_Rubr'
pcd=joink(codedir,'projects',my_project)  #code directory
pjd=joink(projectsdir,my_project)         #projects directory
pjdout=joink(projectsdir,my_project,'out')
collf=joink(pjd,'varsets.can')            
canlistf=joink(pjdout,'list')                
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


    def insertProp(self,keys,prop):
        """insert property 'prop' in 'self' at the level given by 'keys'"""
        pt=self.__dict__
        for key in keys[0:-1]:
            if key in pt.keys(): pt=pt[key]
            else:
                pt[key]={}; pt=pt[key]
        pt[keys[-1]]=prop
        return True

    def delProp(self,*keys):
        pt=self.__dict__
        for key in keys[0:-1]: pt=pt[key]
        del pt[keys[-1]]
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
        except: #extend the object dictionary
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

    def foo(self,conditions):
        pt=self.conditions2dict(conditions); #pt is a dictionary
        return False #True will pickledump

    def qh(self,conditions='T_300 10ns unsolv',task='init'):
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        if task=='init':
            from amber.amber10 import evecs
            pt['evecs']=evecs(jk(sd,'evecs.dat'))
            return True
        if task=='make histograms':
            """make histograms of the mode projections"""
            from amber.amber10 import evecs_projection
            obj=evecs_projection(jk(sd,'modes236.crd'))
            r=obj.fit2Normal() #create histograms, fit to gaussian
            outd=jk(sd,'fit2Normal'); os.system('mkdir -p '+outd)
            buf='#mode mean      std      chi\n'
            for i in range(236):
                buf+='%3d %11.4e %10.4e %10.4e\n'%\
                      (i+1,r['mean'][i],r['std'][i],r['chi'][i])
                buf2='#    q         Hist\n'
                for j in range(len(r['xvals'][i])):
                    buf2+='%11.4e %11.4e\n'%(r['xvals'][i][j],r['histogram'][i][j])
                open(jk(outd,'%03d.dat'%(i+1,)),'w').write(buf2)
            open(jk(outd,'fit2normal.dat'),'w').write(buf)
            obj.close() #allows pickledump
            pt['evecs_proj']=obj
            return True

    def electrostatics(self,task='force on FE3',
                       conditions='T_300 10ns solv'):
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        if task=='force on FE3':
            script="""#!/bin/bash
     
prodd=_PRODD_    #directory where we run the simulations
atom=_ATOM_      #atom number for the Fe ion
dir='rst'    #directory where restart amber files are located
forf='force.dat' #file to output forces

cd $prodd
echo \"#electrostatic force for atom $atom (FE3)\" > $forf
echo '#counter snapshot fx fy fz' >> $forf
n=1
i=1
while [ $n -le _NSNAP_ ];do
  echo \"snapshot $n\"
  
  #gunzip -c $dir/x.$n.gz > current.rst
  /bin/cp $dir/x.$n current.rst
  
  sander -O -i electrostatics.in -p top -c current.rst -ref ref
  v=`grep \"NUMERICAL, ANALYTICAL FORCES\" mdout -A 3|tail -3|tr -s ' '|cut -d ' ' -f 3|tr '\\n' ' '`
  line=`printf '%04d %05d %13.8f %13.8f %13.8f' $i $n $v`
  echo $line >> $forf
  let "n += _JUMP_"
  let "i += 1"
/bin/rm mdout current.rst  #clean-up
done"""
            atomn={'1fhh':'782','ref15-shel':'789'}
            dict={'_PRODD_':sd,
                  '_PCD_':pcd,
                  '_ATOM_':atomn[self.id],
                  '_NSNAP_':'10000',
                  '_JUMP_' :'2',
                  }
            for k,v in dict.items(): script=script.replace(k,v)
            (a,tmpf)=mkstemp(dir='/projects/tmp');
            open(tmpf,'w').write(script)
            os.system('source %s && /bin/rm %s'%(tmpf,tmpf))
            return False
        if task=='force norm on FE3':
            ff=jk(sd,'force.dat') #file containing forces
            ffn=ff+'.new'
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
            open(ffn,'w').write(buf)
            print 'DO: /bin/mv %s %s'%(ffn,ff)
            return False
        elif task=='estimate EP change from average force':
            from amber.amber10 import evecs
            #load average force
            buf=open( jk(sd,'force.dat') ).read()
            m=re.compile('\[([-\d\.\s]+)').search(buf).group(1)
            avf=numpy.array([float(x) for x in m.split()])
            #load mode standard deviations, which are the eigenvalues
            evecsf=jk(sd,'evecs.dat')
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
            open(jk(sd,'EP_estimate.dat'),'w').write(buf)
            return False

        return False

    def namd(self,conditions='T_300 100ns solv',task=None):
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        os.chdir(sd) #go to appropriate directory
        if task=='energy evaluation all':
            """generate X.namd and X.pdb for namd2 energy
            evaluation, then submit namd2 job"""
            #submit energy evaluation
            os.system('namd2 all.namd > all.log &') 
            return False
        elif task=='energy evaluation protein':
            """generate X.namd and X.pdb for namd2 energy
            evaluation, then submit namd2 job"""
            buf=''
            for l in open('all.pdb').readlines():
                if l[0:4]=='ATOM':
                    res=l[17:20]; at=l[12:16]
                    if res in ['Na+','WAT']: l=l[0:62]+'0.00'+l[66:]
                    buf+=l
                    open('protein.pdb','w').write(buf)
            #submit energy evaluation
            buf=open('all.namd').read()
            buf=buf.replace('all.pdb','protein.pdb')
            open('protein.namd','w').write(buf)
            os.system('namd2 protein.namd > protein.log &') 
            return False
        elif task=='energy evaluation protein_all':
            """generate X.namd and X.pdb for namd2 energy
            evaluation, then submit namd2 job"""
            buf=''
            for l in open('all.pdb').readlines():
                if l[0:4]=='ATOM':
                    res=l[17:20]; at=l[12:16]
                    if res in ['NA+','WAT']: l=l[0:62]+'2.00'+l[66:]
                    buf+=l
                    open('protein_all.pdb','w').write(buf)
            #create protein_all.namd                    
            buf=open('FE_all.namd').read()
            buf=buf.replace('FE_all','protein_all')
            open('protein_all.namd','w').write(buf)                    
            #submit energy evaluation
            os.system('namd2 protein_all.namd > protein_all.log &') 
            return False
        elif task=='energy evaluation FE_all':
            """generate X.namd and X.pdb for namd2 energy
            evaluation, then submit namd2 job"""
            buf=''
            for l in open('all.pdb').readlines():
                if l[0:4]=='ATOM':
                    res=l[17:20]; at=l[12:16]
                    l=l[0:62]+'2.00'+l[66:]
                    if res =='FE3': l=l[0:62]+'1.00'+l[66:]
                    buf+=l
                    open('FE_all.pdb','w').write(buf)
            #submit energy evaluation
            os.system('namd2 FE_all.namd > FE_all.log &') 
            return False
        elif task=='energy evaluation FE_WAT':
            """generate X.namd and X.pdb for namd2 energy
            evaluation, then submit namd2 job"""
            buf=''
            for l in open('all.pdb').readlines():
                if l[0:4]=='ATOM':
                    res=l[17:20]; at=l[12:16]
                    l=l[0:62]+'2.00'+l[66:]
                    if res != 'WAT': l=l[0:62]+'0.00'+l[66:]
                    if res =='FE3': l=l[0:62]+'1.00'+l[66:]
                    buf+=l
                    open('FE_WAT.pdb','w').write(buf)
            #create FE_WAT.namd                    
            buf=open('FE_all.namd').read()
            buf=buf.replace('FE_all','FE_WAT')
            open('FE_WAT.namd','w').write(buf)
            #submit energy evaluation
            os.system('namd2 FE_WAT.namd > FE_WAT.log &') 
            return False
        elif task=='energy evaluation FE_S':
            """generate X.namd and X.pdb for namd2 energy
            evaluation, then submit namd2 job"""
            buf=''
            for l in open('all.pdb').readlines():
                if l[0:4]=='ATOM':
                    res=l[17:20]; at=l[12:16]
                    l=l[0:62]+'0.00'+l[66:]
                    if res =='FE3': l=l[0:62]+'1.00'+l[66:]
                    if res=='CYM' and at==' SG ': l=l[0:62]+'2.00'+l[66:]
                    buf+=l
                    open('FE_S.pdb','w').write(buf)
            #create FE_WAT.namd                    
            buf=open('FE_all.namd').read()
            buf=buf.replace('FE_all','FE_S')
            open('FE_S.namd','w').write(buf)
            #submit energy evaluation
            os.system('namd2 FE_S.namd > FE_S.log &') 
            return False
        elif task=='energy evaluation FE_protein':
            """generate X.namd and X.pdb for namd2 energy
            evaluation, then submit namd2 job"""
            buf=''
            for l in open('all.pdb').readlines():
                if l[0:4]=='ATOM':
                    res=l[17:20]; at=l[12:16]
                    l=l[0:62]+'2.00'+l[66:]
                    if res=='FE3': l=l[0:62]+'1.00'+l[66:]
                    if res=='WAT' or res=='Na+': l=l[0:62]+'0.00'+l[66:]
                    buf+=l
                    open('FE_protein.pdb','w').write(buf)
            #create FE_WAT.namd                    
            buf=open('FE_all.namd').read()
            buf=buf.replace('FE_all','FE_protein')
            open('FE_protein.namd','w').write(buf)
            #submit energy evaluation
            os.system('namd2 FE_protein.namd > FE_protein.log &') 
            return False
        elif task=='energy evaluation FE_Na+':
            """generate X.namd and X.pdb for namd2 energy
            evaluation, then submit namd2 job"""
            buf=''
            for l in open('all.pdb').readlines():
                if l[0:4]=='ATOM':
                    res=l[17:20]; at=l[12:16]
                    l=l[0:62]+'0.00'+l[66:]
                    if res=='FE3': l=l[0:62]+'1.00'+l[66:]
                    if res=='Na+': l=l[0:62]+'2.00'+l[66:]
                    buf+=l
                    open('FE_ions.pdb','w').write(buf)
            #create FE_WAT.namd                    
            buf=open('FE_all.namd').read()
            buf=buf.replace('FE_all','FE_ions')
            open('FE_ions.namd','w').write(buf)
            #submit energy evaluation
            os.system('namd2 FE_ions.namd > FE_ions.log &') 
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
        c.gensimdir()
        c.pickleDump()
        print remaining,header
        remaining-=1


def updateCollection():

    """update the collection file"""
    
    cl.updColletion()

def tview(id):
    
    """viewer of object properties in terminal"""

    cl.loadCAN(id).tview()

    
def view(id):

    """viewver of object properties"""

    cl.loadCAN(id).view()
    
def insertProp(*kargs,**kwargs):
    kargs=(kargs[0:-1],kargs[-1])
    cl.cycleOverList('insertProp',*kargs,**kwargs)

def delProp(*kargs,**kwargs):
    cl.cycleOverList('delProp',*kargs,**kwargs)

def foo(*kargs,**kwargs):
    """create seq.dat files
    Example: ... -c id=a0512,switch=hello"""
    cl.cycleOverList('foo',*kargs,**kwargs)

def gensimdir(*kargs,**kwargs):
    """simulation directory"""
    cl.cycleOverList('gensimdir',*kargs,**kwargs)

def qh(*kargs,**kwargs):
    """quasi harmonic bundle"""
    cl.cycleOverList('qh',*kargs,**kwargs)

def electrostatics(*kargs,**kwargs):
    """electrostatics bundle"""
    cl.cycleOverList('electrostatics',*kargs,**kwargs)
            
def namd(*kargs,**kwargs):
    cl.cycleOverList('namd',*kargs,**kwargs)
            

def help(*kargs):

    """list available services"""
    
    services=['datatree','updateCollection','initCollection','view']
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
