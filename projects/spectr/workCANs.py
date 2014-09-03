#!/usr/bin/python

import pdb
import os,sys,re,math,numpy
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir,projectsdir
from utilities.small_utilities import chomp,Bye,junkName
from inputArgs.inputArgs import inpHand

"""project description goes here"""

#global variables
atomic_weight={'H':'1', 'D':'2', 'C':'12', 'N':'14', 'O':'16', 'S':'32'}
joink=os.path.join
isthere=os.path.exists
basename=os.path.basename

my_project='spectr'                   #substitute 'my_project'
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
        pdbf=None
        if self.id in ('GLY','2-GLY','4-GLY'):
            pdbf=joink(self.simdir,self.id+'.pdb')
        elif self.id=='a-GLY':
            """open the unit cell and propagate"""
            a=5.1020; b=11.9709; c=5.4575; beta=math.pi*111.70/180.0
            iGLY=0; iat=0; ires=1; nat=10; unitcell=[]
            for l in open( joink(self.simdir,'GLYCIN02.pdb') ).readlines():
                if 'HETATM' in l:
                    l=l.replace('HETATM','ATOM  ')
                    l=l.replace('    1    ','    '+`ires`+'    ')
                    unitcell.append(l.replace('UNK','ZGL'))
                    iat+=1
                    if iat==nat:
                        iat=0; ires+=1; unitcell.append('TER\n')
            pdbf=joink(self.simdir,'unitcell.pdb')
            open(pdbf,'w').write(''.join(unitcell)+'END\n')
            self.unitcell=pdbf
            buf=''
            iat=1; ires=1
            for ia in range(1,5):
                for ic in range(1,5):
                    t=[a*ia+c*ic*numpy.cos(beta),0,c*ic*numpy.sin(beta)]
                    t=numpy.array(t) #translation vector
#                    pdb.set_trace()
                    for l in unitcell:
                        if 'ATOM  ' in l:
                            r=[float(r_i) for r_i in l[30:54].split()]
                            r='%8.3f%8.3f%8.3f'%tuple(numpy.array(r)+t)
                            l=l[:30]+r+l[54:]
                            l=l[:7]+'%4d'%(iat,)+l[11:]; iat+=1
                            l=l[:24]+'%2d'%(ires,)+l[26:];
                        elif 'TER' in l:
                            ires+=1
                        buf+=l
            buf+='END\n'
            pdbf=joink(self.simdir,'a-GLY.pdb')
            open(pdbf,'w').write(buf)
        self.pdbfile=pdbf
        return True
    
    def get_nat(self):
        """find number of atoms"""
        nat=0
        for l in open(self.pdbfile).readlines():
            if l[0:5]=='ATOM ': nat+=1
        self.nat=nat
        return True
    
    def importPropsFile(self):
        """enter info from properties.dat file"""
        from utilities.readingWritingFiles import importPropFile
        self.propf=joink(self.simdir,'properties.dat')
        self.__dict__.update( importPropFile(self.propf) )
        return True

    def nwchemin(self):
        """create input job"""
        pdb.set_trace()
        #set up system
        buf='start _ID__freq\n\n'
        for k,v in {'_ID_':self.id,'_Q_':`self.zwit['charge']`}.items():
            buf=buf.replace(k,v)
        #read fractional coordinates fron cif file
        buf+="""geometry units agnstroms\n
system crystal
lat_a 5.1020 lat_b 11.9709 lat_c 5.4575
alpha 90.00 beta  111.70 gamma 90.00
end
"""        
        pt=open( joink(self.simdir,self.id+'.cif') ); l=''
        while '_atom_site_fract_z' not in l: l=pt.readline()
        for i in range(self.nat):
            a,b,x,y,z=pt.readline().split();
            buf+=' %s %8.5f %8.5f %8.5f\n'%(b,float(x),float(y),float(z))
        buf+='end\n'
        #define tasks
        buf+="""
basis
 C library sto-3g
 H library sto-3g
 N library sto-3g
end

scf
 uhf, doublet
 print low
end

title \"_ID_ : STO-3G UHF geometry optimization\"

task scf optimize

basis
 C library 6-31g*
 H library 6-31g*
 N library 6-31g*
end

title \"_ID_ : 6-31G* UHF geometry optimization\"

task scf optimize

scf; print none; end

title \"_ID_ : 6-31* UHF frequencies\"

task scf freq
"""
        buf=buf.replace('_ID_',self.id)
        #output task file
        self.nwchemin=joink(self.simdir,self.id+'.nw')
        open(self.nwchemin,'w').write(buf)
        return True
            
        
    def g98in(self,job='opt',conv='',T='4.0'):
        """generate a g98 input file for optimization or
        frequency calculation"""
        types=[]  ;  name=self.id+'_'+job  ;
        title=self.id+'  '+job  ;  hf='rhf'
        Q=self.zwit['charge']
        mult=self.zwit['multiplicity']
        if mult>0: hf='uhf'
        buf='%Mem=64MW\n%chk='+name+'.chk\n'
        if job=='opt':
            buf+='#'+hf+'/6-31G* opt SCF=QC\n\n'
            buf+=title+'\n\n'+\
                  `Q`+' '+`mult`+'\n'
            for line in open(self.pdbfile,'r').readlines():
                if 'ATOM'==line[0:4]:
                    type=line[13]
                    xyz =line[30:54]
                    buf+=type+'   '+xyz+'\n'
                    types.append(type)            
        elif job=='freq':
            buf+='#'+hf+'/6-31G* freq=readisotopes SCF=QC\n\n'+\
                  title+'\n\n'+\
                  `Q`+' '+`mult`+'\n'
            for line in open(self.optPDB,'r').readlines():
                if 'ATOM'==line[0:4]:
                    type=line[13]
                    xyz =line[30:54]
                    buf+=type+'   '+xyz+'\n'
                    types.append(type)            
            buf=buf+'\n'+T+' 1.0\n' #temperature
            for type in types:
                if self.isDeut and type=='H': type='D'
                buf+=atomic_weight[type]+'\n'
        buf+='\n'
        ing98file=joink(self.simdir,name+'.com')
        open(ing98file,'w').write(buf)
        self.__dict__['g98in.'+job]=ing98file
        return True

    def genoptPDB(self):
        """generate PDB from gaussian optimization run"""
        self.optPDB=None
        optfile=self.__dict__['g98out.opt']
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
            fin2=open(self.pdbfile,'r')  ;  line2=fin2.readline()
            while(line2[0:4]!='ATOM'): line2=fin2.readline()
            buf=''
            for i in range(5): line=fin.readline() #advance four lines
            pdb.set_trace()
            while '------' not in line:
                c,a,a,x,y,z=line.split()
                x=float(x) ; y=float(y) ; z=float(z)
                buf+=line2[0:30]+'%8.3f%8.3f%8.3f'%(x,y,z)+line2[54:]
                line=fin.readline()  ; line2=fin2.readline()
                while 'ATOM' not in line2 and line2:
                    buf+=line2; line2=fin2.readline()
            self.optPDB=joink(self.simdir,self.id+'_opt.pdb')
            open(self.optPDB,'w').write(buf)
        return True            
    
    def g98out(self,job='opt'):
        """Run g98"""
        name=self.id+'_'+job
        inf=self.__dict__['g98in.'+job]        #check input file
        if not inf: Bye('no input file')
        outf=joink(self.simdir,name+'.out') #check out file
        self.__dict__['g98out.'+job]=outf
        if isthere(outf):
            cmd='tail -1 '+outf+' |grep "Normal termination"'
            if os.system(cmd):
                os.system('/bin/rm '+outf)
                chkf=joink('./',name+'.chk')#  ;  print outf,chkf
                if isthere(chkf): os.system('/bin/rm '+chkf)
            else:
                Bye('completed')
        cmd= joink(os.environ['HOME'],'bin','g98.csh')+'  '+inf
        print 'PLEASE MANUALLY RUN: "'+cmd+'"'
        return True

    def gen_aCLIMAX(self):
        """generate aCLIMAX frequency file from gaussian output"""
        from gaussian.api import g98freq_to_aCLIMAX
        outf=joink(self.simdir,self.id+'.aCLIMAX')
        pdb.set_trace()
        g98freq_to_aCLIMAX(self.__dict__['g98out.'+job],outf)
        return False

            















    
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
        (1)add missing hydrogens
        (2)solvate box with 10 Angstroms at least around the protein
        (3)add counterions"""
        self.intleap=''
        buf='source leaprc.ff98\n'
        if self.isDeutd: buf='source leaprc.ff98Dsimple\n'
        if self.id in ('1vcx_Fe','1vcx_Fe_D','1fhm','1fhm_D',
                       'ref15-shel','ref15-shel_D'):
            buf+='loadOff FE3.lib\n'
            buf+='loadamberparams FE.frcmod\n'
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
        return False
    
    def quasiharmonic(self,T='300',nvecs=1000,simLength=''):
        """obtain quasiharmonic modes for CA, heavy-atoms,
        hydrogen-atoms, and all-atoms"""
        from amber.Prod import Prod
        from amber.amber10 import ptraj
        pr=self.solution['T_'+T]['Prod'+simLength]
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
analyze matrix mwcovar out %s thermo vecs 0"""%(traj,maskv,nvecs,junkf))
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
    
    def quasiharmonic_movies(self,T='300',simLength=''):
        """create PDB movies for each quasiharmonic mode"""
        from amber.amber10 import evecs
        Prod=self.solution['T_'+T]['Prod'+simLength]
        topfile=self.unsolv['parmfile']
        key='fit2centr.evecsall'; outdir='strip/'+key+'d'
        ev=evecs( joink(Prod.simdir,Prod.strip[key]) )
        outf=joink(Prod.simdir,outdir)
        ev.create_movies( topfile,outdir=outf )
        Prod.strip['fit2centr.evecsalld']=outdir
        Prod.dump()
        return True#False
    
    def projectCRDontoEVECS(self,T='300',simLength=''):
        """For every snapshot, find the projection of each trajectory
        snapshot onto the quasi-harmonic modes"""
        from amber.amber10 import crd,evecs,qhm_proj
        Prod=self.solution['T_'+T]['Prod'+simLength]
        crdf=joink(Prod.simdir,Prod.strip['fit2centr'])
        evecf=joink(Prod.simdir,Prod.strip['fit2centr.evecsall'])
        outf=evecf.replace('.all.dat','.all.project.dat')#store projections
        qhmp=qhm_proj(outf)
        qhmp.populate(crdf,Prod.strip['nat'],evecf)
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
    
    def genIsDeutd(self):
        """flag as deuterated"""
        self.isDeutd=False
        if self.id in ('1vcx_Fe_D','1fhm_D','ref15-shel_D'):
            self.isDeutd=True
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
        pdb.set_trace()
        os.system(cmd)
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
        elif pdb=='uns.premin':
            self.unsolv['onlyHpreminPDB']=self.unsolv['preminPDB'].replace('.pdb','.onlyH.pdb')
            open(self.unsolv['onlyHpreminPDB'],'w').write(buf)
        return False
        
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
        Prod=self.solution['T_'+T]['Prod'+simLength] #Production object
        if 'crdfile' in pH.keys():
            pH['crd.fit2centr']=pH['crdfile']
            del pH['crdfile']
        #for keyword in ('fit2centr','fit2prev'):
        for keyword in ('crdfile','fit2prev','fit2centr'):
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
        pr.BfactorsByMask(mask='@H*')
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

    def getSIIScF_sts(self,T='300',simLength='',
                      tl=40,recover=False):
        """calculate Self-Intermediate Incoherent Scattering Function
        nt:number of time points for which to calculate the function
        tl: total time desired, in picoseconds"""
        from scatter.scatter import SIIScF,storeSIIScF,fourier
        from amber.amber10 import crd
        tl=float(tl)
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']
        nt=int(tl/dt) ; print 'number of time points=',nt
        tt=[]
        for i in range(nt): tt.append(i*dt)
        qq=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] #Angstroms^(-1)
        nq=len(qq)
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
        gaussianenvelope=[]
        for tti in tt: gaussianenvelope.append(convT(tti))
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
                    qq,tt)
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
        col: column where to get the data (def: last column)
        ss: temporal resolution, in picoseconds (def: 'BASIS'=521ps)"""
        from scatter.scatter import fourier
        factor=1/4.14 #from picoseconds to 10E-6eV)
        isf=[]; ts=[]
        for line in open(infile,'r').readlines():
            if line[0]=='#': continue #this is a comment line
            items=[float(x) for x in line.split()]
            ts.append(items[0]); isf.append(items[col])
        omegas=[]
        for i in range(len(isf)):
            t=ts[i]; omegas.append(t*factor)
            env=convT(t,ss=ss)/convT(0,ss=ss) #"convolute" in time-domain
            #print t,isf[i],isf[i]*env
            isf[i]*=env            
        dsf=fourier(isf).real
        buf=''
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
        
    def get_msd2(self,T='300',tl=40,simLength=''):
        """calculate time-dependent average atomic fluctuation
        tl: total time desired, in picoseconds"""
        from scatter.scatter import msd2
        from amber.amber10 import crd
        tl=float(tl)
        Prod=self.solution['T_'+T]['Prod'+simLength]
        pH=self.solution['T_'+T]['onlyH'+simLength]
        crdf=pH['crd.fit2prev']
        crdobj=crd(crdf,self.Nhydr)
        skipframes=0
        print 'loading 1 frame every ',str(1+skipframes),' frames'
        crdobj.loadframes(skip=skipframes)
        #dt: time in between loaded frames, in picoseconds
        dt=(1+skipframes)*Prod.times['twx']
        nt=int(tl/dt) ; print 'number of time points=',nt
        scattering=msd2(crdobj.frames,dt,tl,nt) #scattering[Nhydr][nt]
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

    def humpness_msd2(self,T='300',simLength=''):
        """find MSD2 humpness below 1ps for every MSD2(t) plot"""
        from stats.statistics import linreg
        pH=self.solution['T_'+T]['onlyH'+simLength]
        pH['msd2']=joink(pH['simdir'],'msd2')
        buf='#i-hydr slope intercept corr-coeff  tmax msd2max hump\n'
        humps=[]
        for i in range(self.Nhydr):
            infile=joink(pH['msd2'],'%04d.dat'%(i+1,))
            t=[] ; t0=0.0 ; n=0 ; nl=''; msd2=[] ; m0=0.0 ; xl=1.#1picosec
            for line in open(infile).readlines()[1:]:
                x,y=line.split(); x=float(x) ; y=float(y)
                if x<xl and y>m0:
                    m0=y ; t0=x
                if not nl:
                    if x>xl: nl=n #when time becomes greater than xl?
                else:
                    t.append(float(x)); msd2.append(float(y))
                n+=1
            regr=linreg(t,msd2); a=regr['slope']; b=regr['intercept']
            hump=m0-(b+a*t0); humps.append(hump)
            buf+='%3d %f %f %f %f %f %f\n'%(i+1,a,b,regr['corrcoef'],
                                            t0,m0,hump)
        pH['humpness']=joink(pH['msd2'],'humpness.dat')
        open(pH['humpness'],'w').write(buf)
        #create PDB with humpness values as B-factors.
        #Non-H atoms are given a 0.0 humpness
        buf='' ; i=0
        for line in chomp(open(self.unsolv['preminPDB']).readlines()):
            if line[0:4]=='ATOM':
                l=len(line)
                if len(line)<66: line=line+' '*(66-l)
                bf=0.0
                if line[13]=='H' and i<self.Nhydr:
                    bf=humps[i] ;i+=1
                line=line[0:60]+'%5.2f'%(bf,)+line[66:]
            buf+=line+'\n'
        pH['pdbhump']=joink(pH['msd2'],'humpness.pdb')
        open(pH['pdbhump'],'w').write(buf)
        return True

    def fract_hump(self,T='300',simLength=''):
        """find fraction of H's (D's) important for the hump in the
        MSD^2(t) plot"""
        pH=self.solution['T_'+T]['onlyH'+simLength]
        #store the humps
        humps={}
        for line in open(pH['humpness']).readlines()[1:]:
            l=line.split() ; i=int(l[0]) ; h=float(l[-1]) ; humps[i]=h
        #sort in decreasing humpness
        pairs=[ (h,i) for (i,h) in humps.items()]
        pairs.sort(); pairs.reverse()
        #create average MSD^2 plot by increasing the number of hydrogens
        q=len(open(joink(pH['msd2'],'0001.dat')).readlines()[1:])
        n=len(pairs); MSD2=[0.0]*q; bufs=['']*q
        for pair in pairs:
            ii=pair[1]
            ls=open(joink(pH['msd2'],'%04d.dat'%(ii,))).readlines()[1:]
            j=0
            for l in ls:
                (t,m)=l.split(); #print j
                MSD2[j]+=float(m)
                j+=1
        for j in range(q):
            MSD2[j]/=n; bufs[j]+=' %5.3f'%(MSD2[j],)
        for pair in pairs[0:-1]:
            ii=pair[1];print ii
            ls=open(joink(pH['msd2'],'%04d.dat'%(ii,))).readlines()[1:]
            j=0
            for l in ls:
                (t,m)=l.split()
                MSD2[j]=(MSD2[j]*n-float(m))/(n-1)
                bufs[j]+=' %5.3f'%(MSD2[j],)
                j+=1
            n-=1
        buf='#t MSD2...\n'+'\n'.join(bufs)
        fract_hump
        open( joink(pH['msd2'],'fract_hump.dat'),'w').write(buf)
        return False
            
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

def nwchemin(*kargs,**kwargs):
    """enter info from properties.dat file"""
    cl.cycleOverList('nwchemin',*kargs,**kwargs)

def importPropsFile(*kargs,**kwargs):
    """enter info from properties.dat file"""
    cl.cycleOverList('importPropsFile',*kargs,**kwargs)

def g98in(*kargs,**kwargs):
    """generate a g98 input file for optimization or
    frequency calculation"""
    cl.cycleOverList('g98in',*kargs,**kwargs)

def genoptPDB(*kargs,**kwargs):
    """generate PDB from gaussian optimization run"""
    cl.cycleOverList('genoptPDB',*kargs,**kwargs)

def g98out(*kargs,**kwargs):
    """Run g98"""
    cl.cycleOverList('g98out',*kargs,**kwargs)

def gen_aCLIMAX(*kargs,**kwargs):
    """generate aCLIMAX frequency file from gaussian output"""
    cl.cycleOverList('gen_aCLIMAX',*kargs,**kwargs)

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

def relax(*kargs,**kwargs):
    """minimize and relax the premin pdb structure"""
    cl.cycleOverList('relax',*kargs,**kwargs)

def are_orthogonal(*kargs,**kwargs):
    """check if quasiharmonic modes are orthogonal"""
    cl.cycleOverList('are_orthogonal',*kargs,**kwargs)

def quasiharmonic(*kargs,**kwargs):
    """perform quasiharmonic analysis"""
    cl.cycleOverList('quasiharmonic',*kargs,**kwargs)

def getNvecs(*kargs,**kwargs):
    cl.cycleOverList('getNvecs',*kargs,**kwargs)

def quasiharmonic_movies(*kargs,**kwargs):
    cl.cycleOverList('quasiharmonic_movies',*kargs,**kwargs)
    
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

def get_nat(*kargs,**kwargs):
    """find number of atoms"""
    cl.cycleOverList('get_nat',*kargs,**kwargs)

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
    
def av_H_msd(*kargs,**kwargs):
    """average mean square displacement of H or D from trajectory
    fitted to centroid"""
    cl.cycleOverList('av_H_msd',*kargs,**kwargs)

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
