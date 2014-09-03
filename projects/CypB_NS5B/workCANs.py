#!/usr/bin/python

import pdb
import os,sys,re,numpy
from copy import deepcopy
from tempfile import mkstemp,mkdtemp
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

my_project='CypB_NS5B'                      #substitute 'my_project'
pcd=joink(codedir,'projects',my_project)  #code directory
pjd=joink(projectsdir,my_project)         #projects directory
pjdout=joink(projectsdir,my_project,'out')
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

def translate(task='start',alpha=1.8,inpdb=None,outpdb=None,at1=None,at2=None):
    """You have to select the translation vector by defining initial
(at1) and final (at2) atoms in complex.pdb, as well as the magnitude
of the translation (alpha)."""
    alpha=float(alpha)
    if task=='start':
        pin=open(inpdb); l=pin.readline()
        while at1 not in l: l=pin.readline() #;print l
        start=[float(x) for x in l[30:54].split()]; #print start
        pin.seek(0); l=pin.readline()
        while at2 not in l: l=pin.readline()
        end=[float(x) for x in l[30:54].split()]; #print end
        pin.close()        
    t=[]
    for i in range(3): t.append(alpha*(end[i]-start[i]))
    buf= 'REMARK translataion %6f %6f %6f\n'%(t[0],t[1],t[2])
    pin=open(inpdb); l=pin.readline(); flag=True
    while(l):
        if 'TER' in l: flag=False
        if(flag and l[0:5]=='ATOM '):
            xyz=[float(x) for x in l[30:54].split()]
            for i in range(3): xyz[i]+=t[i]
            l=l[0:30]+'%8.3f%8.3f%8.3f'%(xyz[0],xyz[1],xyz[2])+l[54:]
        buf+=l
        l=pin.readline()
    open(outpdb,'w').write(buf)
    return False

def adjust(task='by distance',delta=1.0,inpdb=None,outpdb=None):
    from math import sqrt
    delta=float(delta)
    if task=='by distance':
        """reduce distance between 1V9T/R48/CA and 1C2PA/P540/CA by
        delta by translation of 1V9T towards 1C2PA"""
        #find distance
        pin=open(inpdb); l=pin.readline()
        while l:
            if 'CA  ARG    48' in l:
                start=[float(x) for x in l[30:54].split()]
            if 'CA  PRO   701' in l:
                end=[float(x) for x in l[30:54].split()]
                break
            l=pin.readline()
        pin.close()
        #find translation
        t=[]; dd=0.0
        for i in range(3):
            d=end[i]-start[i]; dd+=d*d; t.append(d)
        d=sqrt(dd); alpha=delta/d
        for i in range(3): t[i]*=alpha
        #translate 1V9T
        buf=''; TERflag=False
        pin=open(inpdb); l=pin.readline()
        while l:
            if not TERflag and 'ATOM ' in l:
                xyz=[float(x) for x in l[30:54].split()]
                for i in range(3): xyz[i]+=t[i]
                l=l[0:30]+'%8.3f%8.3f%8.3f'%(xyz[0],xyz[1],xyz[2])+l[54:]
            if 'TER ' in l: TERflag=True #signal end of 1V9T
            buf+=l
            l=pin.readline()
        open(outpdb,'w').write(buf)
        return False
    elif task=='bfactor correlations':
        """calculate correlation coefficient between CA bfactor of
        1QUV.pdb, 1C2PA_MET.pdb, and NS5/relax_insertion/Prod/strip"""
        #Load bfactor values
        bff={'MD':pjdout+'/NS5/relax_insertion/Prod/strip/bfactor_CA.dat',
             '1C2PA':pjdout+'/1C2PA_MET/1C2PA_MET_bfactor.dat',
             '1QUV':pjdout+'/1QUV/prepare/1QUV_bfactor.dat'}
        buf='#                     a_s      b_s         r        tt       stderr\n'
        N=600
        bf={'MD':[-1.]*N, '1C2PA':[-1.]*N, '1QUV':[-1.]*N}
        for (key,val) in bff.items():
            for l in open(val):
                if l[0]=='#': continue
                ires,b=l.split(); ires=int(ires); b=float(b)
                bf[key][ires]=b
        #Find pair correlations
        bfx={}
        for key1 in ('MD', '1C2PA', '1QUV'):
            for key2 in ('MD', '1C2PA', '1QUV'):
                if key1==key2: continue
                bf1=bf[key1]; bf2=bf[key2]
                #Exclude common regions
                b1=[]; b2=[]
                for i in range(N):
                    if bf1[i]<0 or bf2[i]<0: continue
                    b1.append( bf1[i] )
                    b2.append( bf2[i] )
                #calculate correlation between b1 and b2
                from scipy import stats
                vv = [ '%8.3f'%x for x in stats.linregress(b1,b2) ]
                buf+='#%8s %8s %s\n'%(key1, key2, '  '.join(vv) )
                #adjust bfactors of 1QUV and the 10ns simulation to
                #those of #1C2PA_MET
                if key2=='1C2PA':
                    bfx[key1]=deepcopy( bf[key1] )
                    (a,b)=stats.linregress(b1,b2)[:2]
                    for i in range(N):
                        if bfx[key1][i] > 0:
                            bfx[key1][i] = b + a * bfx[key1][i]
        #print adjusted bfactors
        buf += '#Adjusted bfactors to those of 1C2PA_MET\n'
        buf += '#1C2PA_MET  1QUV  10ns-MD\n'
        for i in range(1,N):
            buf += '%3d %6.2f %6.2f %6.2f\n'\
                   %(i,bf['1C2PA'][i],bfx['1QUV'][i],bfx['MD'][i])
        #pdb.set_trace()
        print buf
        return False

    return False

def genRST(task=None,inpdb=None,outf=None):
    """generate NRM restraints"""
    if task=='1V9T_minus_T143_V152':
        #read CyPB protein
        buf=''; pin=open(joink(pjd,'prepare/1V9T_minus_T143_V152.pdb'))
        l=pin.readline()
        while l[0:3] != 'TER': buf+=l; l=pin.readline()
        open('junk.pdb','w').write(buf)
        from amber.amber10 import genRST
        buf=genRST('junk.pdb',task='by CB map',co=7.00)
        outf=joink(pjd,'prepare/1V9T_minus_T143_V152.RST')
        open(outf,'w').write(buf)
        return False
    elif task=='1V9T_minus_T93_A99':
        #read CyPB protein
        buf=''; pin=open(joink(pjd,'prepare/1V9T_minus_T93_A99.pdb'))
        l=pin.readline()
        while l[0:3] != 'END': buf+=l; l=pin.readline()
        open('junk.pdb','w').write(buf)
        from amber.amber10 import genRST
        buf=genRST('junk.pdb',task='by CB map',co=7.00)
        outf=joink(pjd,'prepare/1V9T_minus_T93_A99.RST')
        open(outf,'w').write(buf)
        return False
    elif task=='1CYN_1C2PA_noclash':
        #pdb.set_trace()
        pin=open(joink(pjd,'prepare/1CYN_1C2PA_noclash_premin.pdb'))
        l=pin.readline(); buf=''
        while l[0:3] != 'TER': buf+=l; l=pin.readline()
        open('junk.pdb','w').write(buf)
        from amber.amber10 import genRST
        buf=genRST('junk.pdb',task='by CB map',co=7.00)
        outf=joink(pjd,'prepare/1CYN_1C2PA_noclash_premin.RST')
        open(outf,'w').write(buf)
        return False
    #if we provide an input pdb file, and specify output file
    from amber.amber10 import genRST
    open(outf,'w').write( genRST(inpdb,task='by CB map',co=7.00) )
    return False

def geometries(inp=None,inp2=None,out=None,task=None):
    """analyse geometrical properties
    """
    if task=='Ct surface':
        """calculate exposed and buried surface of C-terminal"""
        import pymol
        import pymol.cmd as pcmd
        
        pymol.finish_launching()
        if not inp2: skip=0
        else: skip = int( inp2 )
        buf='#frame  ct-total Ct-exposed Ct-buried\n'
        for iframe in range(1,10000,skip):
            ff="./pdb/frame.%d"%iframe
            if not isthere(ff): break
            print 'iframe=%d'%iframe
            pcmd.load(ff)
            pcmd.set("dot_solvent",1)
            pcmd.set("dot_density",3)

            if inp=='holo':
                pcmd.select("ct", "resi 708-748")
                pcmd.create("ct2","resi 708-748")
                pcmd.select("1cyn","resi 1-178")
                pcmd.remove("1cyn")
            elif inp=='apo':
                pcmd.select("ct", "resi 530-570")
                pcmd.create("ct2","resi 530-570")
            act=pcmd.get_area("ct")
            act2=pcmd.get_area("ct2")

            buf+='%4d %8.1f %8.1f %8.1f\n'%(iframe, act2, act, act2-act)
            pcmd.delete("all")
        pcmd.quit()
        open("./ct_surf.dat","w").write(buf)
        return False
    elif task=='bl surface':
        """calculate exposed and buried surface of beta-loop"""
        import pymol
        import pymol.cmd as pcmd

        pymol.finish_launching()
        if not inp2: skip=0
        else: skip = int( inp2 )
        buf='#frame  bl-total Bl-exposed Bl-buried\n'
        for iframe in range(1,10000,skip):
            pdbf="./pdb/frame.%d"%iframe
            if not isthere( pdbf ): break
            print 'iframe=%d'%iframe
            pcmd.load(pdbf)
            pcmd.set("dot_solvent",1)
            pcmd.set("dot_density",3)

            if inp=='holo':            
                pcmd.select("bl", "resi 708-748")
                pcmd.create("bl2","resi 708-748")
                pcmd.select("1cyn","resi 1-178")
                pcmd.remove("1cyn")
            elif inp=='apo':
                pcmd.select("bl", "resi 438-457")
                pcmd.create("bl2","resi 438-457")               
                
            abl=pcmd.get_area("bl")
            abl2=pcmd.get_area("bl2")

            buf+='%4d %8.1f %8.1f %8.1f\n'%(iframe, abl2, abl, abl2-abl)
            pcmd.delete("all")
        pcmd.quit()
        open("./bl_surf.dat","w").write(buf)
        return False
    
    elif task=='P718 Phi Psi':
        
        import pymol
        import pymol.cmd as pcmd
        pymol.finish_launching()
        buf='#Frame P718-Phi P718-psi\n'
        for iframe in range(1,2501):
            pcmd.load("./pdb/frame.%d"%iframe)
            phi=pcmd.get_dihedral("resi 717 and name C", "resi 718 and name N", "resi 718 and name CA", "resi 718 and name C")
            psi=pcmd.get_dihedral("resi 718 and name N", "resi 718 and name CA", "resi 718 and name C","resi 719 and name N")
            buf+='%4d %7.2f %7.2f\n'%(iframe,phi,psi)
            pcmd.delete("all")
        open('phi_psi_P718.dat','w').write(buf)
        pcmd.quit()
        return False

    elif task=='corrplane angle':

        from math import acos,pi
        buf = ''
        v0 = None
        iframe = 1
        for l in open('./corrplane_P718_ring.dat').readlines():
            if l[0]=='#': continue
            vi = [ float(x) for x in l.split()[1:] ]
            if iframe == 1:
                v0 = vi
                buf += '%4d %6.1f\n'%(1.0,0.0)
            else:
                c = v0[0]*vi[0] + v0[1]*vi[1] + v0[2]*vi[2] 
                buf += '%4d %6.1f\n'%(iframe, 180.0*acos( abs(c) )/pi )
            iframe += 1
        open('./corrplane_P718_ring_angle.dat','w').write( buf )

    elif task=='create list pairs':

        cof = 1.0 #cut-off frequency
        
        bl_r = range(438,457)  #range of beta-loop
        ct_r = range(530,571)  #range of C-terminal
        vl_r = range(1,438)+range(458,530) #range of "vital" part
        if inp=='holo':
            bl_r    = range(616,635)
            ct_r    = range(708,748)
            vl_r = range(178,616)+range(636,708)

        ct_ct=[]
        ct_bl=[]
        ct_vl=[]
        bl_bl=[]
        bl_vl=[]
        for line in open('residue_map.dat').readlines():
            if line[0]=='#': continue
            x = line.split(); a=int(x[0]); b=int(x[1]); c=float(x[2])

            if a==b or (abs(a-b)>1 and c < cof): continue
            #note that residue_map.dat contains both (a,b) and (b,a)
            if   a in ct_r and b in vl_r: ct_vl.append( (a,b) )
            if   a in ct_r and b in ct_r:
                if a > b: continue
                ct_ct.append( (a,b) )
            elif a in ct_r and b in bl_r: ct_bl.append( (a,b) )
            elif a in bl_r and b in bl_r:
                if a > b: continue
                bl_bl.append( (a,b) )
            elif a in bl_r and b in vl_r: bl_vl.append( (a,b) )


        buf=''
        for (a,b) in ct_ct: buf+='%3d %3d\n'%(a,b)
        open('ct_ct_pairs.list','w').write(buf)

        buf=''
        for (a,b) in ct_bl: buf+='%3d %3d\n'%(a,b)
        open('ct_bl_pairs.list','w').write(buf)
           
        buf=''
        for (a,b) in ct_vl: buf+='%3d %3d\n'%(a,b)
        open('ct_vl_pairs.list','w').write(buf)
           
        buf=''
        for (a,b) in bl_bl: buf+='%3d %3d\n'%(a,b)
        open('bl_bl_pairs.list','w').write(buf)

        buf=''
        for (a,b) in bl_vl: buf+='%3d %3d\n'%(a,b)
        open('bl_vl_pairs.list','w').write(buf)

    elif task=='compare secstruct profiles':
        """compare how much the secondary structure profile for each
        residue of the C-terminal differs from apo_trans to apo_cis """
        from amber.amber10 import secstruct
        #load the secondary structure profiles from ptraj secstruc output
        s1 = secstruct(); s2=secstruct()
        s1.load(sumFile=inp); s2.load(sumFile=inp2)
        #select residues 540-570, which correspond to indexes 539-569
        sel = numpy.arange(570)
        #calculate the difference
        x = s1.varX(s2,sel1=sel,sel2=sel)
        y = numpy.diag(x)
        buf='#resNum diff-profile\n'
        for i in range(len(sel)): buf+='%3d %7.4f\n'%(1+sel[i],y[i])
        open( out, 'w' ).write(buf)
    return False

def namd(task=None,inp=None,out=None,group=None):
    if task=='prepare input files':

        #find both groups
        if inp=='holo':
            labels={'ct':range(718,745), 'bl':range(616,636),
                    'vital': range(179,616)+range(636,708),
                    'solv':range(749,10000)+['****',],
                    'ions':range(749,772),
                    'water':range(772,10000)+['****',],
                    '1cyn':range(1,179),'prot':range(1,749),
                    'all':range(1,100000)+['****',],
                    }
        elif inp=='apo':
            labels={'ct':range(530,571), 'bl':range(438,458),
                    'vital':range(1,438)+range(458,530),
                    'solv':range(571,10000)+['****',],
                    'ions':range(571,589),
                    'water':range(589,1000)+['****',],
                    'prot':range(1,749),
                    'all':range(1,10000)+['****',],
                    }
        intra = False
        a,b=group.split('-')
        if a == b: intra = True
        if a in labels.keys(): la = labels[a]
        else: la = [ int( a ), ] #just a single residue
        if b in labels.keys(): lb = labels[b]
        else: lb = [ int( b ), ]
        #create pdb files
        buf = ''; count=0
        for l in open('namd/first_frame.pdb').readlines():
            if l[:4]=='ATOM':
                ires = l[22:26]
                try:
                    ires = int(ires)
                except:
                    pass
                if ires in la:
                    l = l.replace('  0.00  0.00', '  0.00  1.00')
                elif ires in lb:
                    if intra:
                        l = l.replace('  0.00  0.00', '  0.00  1.00')
                    else:
                        l = l.replace('  0.00  0.00', '  0.00  2.00')
                #print str(count), str(ires); count+=1
            buf += l
        open('namd/%s.pdb'%group,'w').write( buf )
        #create namd input file
        if intra:
            buf = open('namd/template_intra.namd').read()
        else:
            buf = open('namd/template_inter.namd').read()
        buf = buf.replace('_TEMPLATE_',group)
        open('namd/%s.namd'%group,'w').write( buf )
    elif task=='energy evaluation':
        currd=os.getcwd()                  
        wd=mkdtemp(dir='/projects/tmp') #run namd here
        print 'Temporary directory: ',wd
        #move input files, then run namd
        cmd ='cd namd; cp %s.namd %s.pdb box.xsc %s && '%(group,group,wd)
        cmd+='cd %s && '%wd
        cmd+='ln -s %s/top top && ln -s %s/dcd dcd; '%(currd,currd)
        cmd+='namd2 %s.namd > %s.log && '%(group,group)
        cmd+='mv %s.log %s/namd/'%(group,currd)
        os.system(cmd)
        os.system('rm -rf '+wd) #remove working directory
    elif task=='parse log file':
        buf  = '#    N   E-bonded E-nonbonded       E      electr'
        buf += '       VdW\n'
        avEb=0.0; avEnb=0.0; avE=0.0
        avELECT =0.0; avVDW=0.0
        for l in open('namd/%s.log'%group).readlines():
            if l[0:6]=='ENERGY':
                TS,BOND,ANGLE,DIHED,IMPRP,ELECT,VDW,BOUNDARY,MISC,\
                KINETIC,TOTAL,TEMP,TOTAL2,TOTAL3,TEMPAVG,PRESSURE,\
                GPRESSURE,VOLUME,PRESSAVG,GPRESSAVG\
                = [ float(x) for x in l.split()[1:] ]
                eb = BOND + ANGLE + DIHED + IMPRP #bonded energy
                enb = ELECT + VDW
                e = eb + enb
                avEb+=eb; avEnb+=enb; avE+=e
                avELECT+=ELECT; avVDW+=VDW
                buf+=' %5d %10.3f %10.3f %10.3f %10.3f %10.3f\n'\
                      %(1+TS,eb,enb,e,ELECT,VDW)
        TS+=1; #print TS
        avEb/=TS; avEnb/=TS; avE/=TS; avELECT/=TS; avVDW/=TS;
        buf2='#     avEb      avEnb       avE      avELECT      avVDW\n'
        buf2+='#%10.3f %10.3f %10.3f %10.3f %10.3f\n'%\
               (avEb,avEnb,avE,avELECT,avVDW)
        open('namd/%s.dat'%group,'w').write(buf2+buf)
    elif task=='parse forces from log file':
        from math import sqrt
        buf  = '#    N   Ex Ey Ez E Vx Vy Vz V Fx Fy Fz F\n'
        avExyz=[0.,0.,0.]; avE=0.
        avVxyz=[0.,0.,0.]; avV=0.
        avFxyz=[0.,0.,0.]; avF=0.
        for l in open('namd/%s.log'%group).readlines():
            if 'PAIR INTERACTION: STEP:' in l:
                (z,a,b,N,c,vx,vy,vz,d,ex,ey,ez) = l.split()
                N=int(N)
                
                ex=float(ex);ey=float(ey);ez=float(ez);
                avExyz[0]+=ex;avExyz[1]+=ey;avExyz[2]+=ez
                e=sqrt(ex*ex+ey*ey+ez*ez)
                avE+=e

                
                vx=float(vx);vy=float(vy);vz=float(vz);
                v=sqrt(vx*vx+vy*vy+vz*vz)
                avVxyz[0]+=vx;avVxyz[1]+=vy;avVxyz[2]+=vz
                avV+=v
                
                fx=vx+ex; fy=vy+ey; fz=vz+ez
                f=sqrt(fx*fx+fy*fy+fz*fz)
                avFxyz[0]+=fx;avFxyz[1]+=fy;avFxyz[2]+=fz
                avF+=f

                buf+=' %5d'%(1+N,)
                for x in (ex,ey,ez,e,vx,vy,vz,fx,fy,fz,f): buf+=' %8.3f'%x
                buf+='\n'
        buf2 ='#     avExyz      avE       avVxyz      avV   avFxyz  avF\n'
        avExyz.append(avE); avVxyz.append(avV); avFxyz.append(avF)
        for x in avExyz+avVxyz+avFxyz: buf2+=' %8.3f'%(x/(1+N),)
        open('namd/%s_forces.dat'%group,'w').write(buf2+'\n'+buf)

    elif task =='analyze list pairs':

        from utilities.readingWritingFiles import read_column
        from scipy import polyfit
        for prefix in 'ct_bl ct_ct ct_vl bl_bl bl_vl'.split():
            buf = '#%s\n'%prefix
            buf+= '#group     <E>    E_std     DD  abs(DD)/E_std\n'
            for pair in open('%s_pairs.list'%prefix).readlines():
                group='-'.join( pair.split() )
                datf = 'namd/%s.dat'%group #.dat file
                X = [ float(x) for x in read_column(datf,1) ]
                Ns = numpy.array( X ) #Frame number
                X = [ float(x) for x in read_column(datf,4) ]
                Es = numpy.array( X ) #Energy values
                avE = numpy.average( Es )

                stE = numpy.std( Es )
                (a,b) = polyfit(Ns,Es,1)#linear regression
                dd = a * Ns[-1]         #Energy change by linear regression
                buf += '%7s %7.2f %7.2f %7.2f %5.2f\n'%(group, avE, stE, dd, abs(dd)/stE)
            open('%s_pairs.dat'%prefix, 'w').write( buf )
            pdb.set_trace()
    return False

def energetics(task=None,inp=None,out=None):
    if task=='output CT-rest polar interactions':
        """polar interactions between C-terminal and rest of NS5B
        determined using pymol and
        $out/NS5/drag2dock/Equil_II/md16/crd.md16_unsolv.pdb"""

        #backbone-backbone
        BB = ( (720,642), (719,642), (734,626), (733,626), (738,627),
               (738,629), (739,269), (740,629) )
        #backbone-sidechain
        BS = ( (711,706), (723,643), (724,643), (734,466), (736,350),
               (737,350), (741,631), (727,630), (712,705), (715,702) )
        #sidechain-sidechain
        SS = ( (737,346), (737,350), (744,636), (722,639) )
        
        #directories where energyMapsByResidue.dat files are stored
        dirs = ('Cis', 'Trans', '1c2pa/Cis', '1c2pa/Trans')
        
        #iterate over all pairs and find average total energy values
        buf = '#Holo-Cis Holo-Trans Apo-Cis Apo-trans\n'
        pairs = { 'BB':BB, 'BS':BS, 'SS':SS }
        for k,v in pairs.items():
            for (a,b) in v:
                if a > b: t=a; a=b; b=t
                buf += ' %3d %3d'%(a,b)
                for dir in dirs: #iterate over all files
                    e = 0.0
                    #recall the shift between Apo and Holo numbering
                    a2=a; b2=b
                    if '1c2pa' in dir: a2-=178; b2-=178
                    p = re.compile('\s%d\s+%d\s+'%(a2,b2))
                    
                    f = dir+'/energyMapsByResidue.dat'
                    if isthere( f ):
                        for l in open(f).readlines():
                            if p.search( l ):
                                e = float( l.split()[-1] )
                                break
                    buf += ' %8.3f'%e
                buf += '\n'
        pdb.set_trace()
        open( 'init_pairs_energetics.dat', 'w' ).write( buf )
        return False
    elif task=='output bl-rest polar interactions':
        """polar interactions between beta-loop and rest of NS5B
        determined using pymol and
        $out/NS5/drag2dock/Equil_II/md16/crd.md16_unsolv.pdb"""

        #backbone-backbone
        BB = ( (633,636), (634,637), (626,733), (626,734), (629,738),
               (627,738), (629,740), (611,616) )
        #backbone-sidechain
        BS = ( (584,622), (631,741) )
        #sidechain-sidechain
        SS = ( (373,626), (633,746), (633,695) )
        
        #directories where energyMapsByResidue.dat files are stored
        dirs = ('Cis', 'Trans', '1c2pa/Cis', '1c2pa/Trans')
        
        #iterate over all pairs and find average total energy values
        buf = '#Holo-Cis Holo-Trans Apo-Cis Apo-trans\n'
        pairs = { 'BB':BB, 'BS':BS, 'SS':SS }
        for k,v in pairs.items():
            for (a,b) in v:
                if a > b: t=a; a=b; b=t
                buf += ' %3d %3d'%(a,b)
                for dir in dirs: #iterate over all files
                    e = 0.0
                    #recall the shift between Apo and Holo numbering
                    a2=a; b2=b
                    if '1c2pa' in dir: a2-=178; b2-=178
                    f = dir+'/energyMapsByResidue_bl_rest.dat'
                    for l in open(f).readlines():
                        if l[0]=='#': continue
                        items = l.split()
                        a3=float(items[0]); b3=int(items[1])
                        if (a3==a2 and b3==b2) or (a3==b2 and b3==a2):
                            e = float( items[-1] )
                            break
                    buf += ' %8.3f'%e
                buf += '\n'
        open( 'bl_rest_energetics.dat', 'w' ).write( buf )
        return False
    return False
   
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
