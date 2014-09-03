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

my_project='Rubr_Hydr'                      #substitute 'my_project'
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
    
def preparing(inp=None,out=None,task=None):
    if task=='prune out WAT':
        """remove all WAT residues with coordinates outside of the
        primitive cell ( .0-34.374, 0.0-35.339, 0.0-44.111) """

        Lx=34.374; Ly=35.339; Lz=44.111
        tx=21.622; ty=11.115; tz=15.516
        def doWeRemove(line):
            x=float(line[30:38])
            y=float(line[38:46])
            z=float(line[46:54])
            if x+tx<0 or x+tx>Lx: return True
            if y+ty<0 or y+ty>Ly: return True
            if z+tz<0 or z+tz>Lz: return True
            return False

        inp='%s/preparing/primitive_cell_tmp.pdb'%pjd
        lines=open(inp).readlines()
        nline=len(lines)
        line=lines[0]
        iline=1
        buf=''
        while iline < nline:
            if line[13:20]=='O   WAT':
                buf2=line
                remove = False
                if doWeRemove(line): remove=True
                line=lines[iline]
                iline+=1
                buf2+=line
                if doWeRemove(line): remove=True
                line=lines[iline]
                iline+=1
                buf2+=line
                if doWeRemove(line): remove=True
                line=lines[iline] #Ter signal
                iline+=1
                if remove==False: buf+=buf2
            else:
                buf+=line
            line=lines[iline]
            iline+=1
        out='%s/preparing/primitive_cell.solv.pdb'%pjd
        open(out,'w').write(buf)
        return False

def scattering(inp=None,inp2=None,inp3=None,
               out=None,out2=None,task=None ):
    
    from scatter.scatter import SIIScF,storeSIIScF,fourier,fourierRealSymmetric,tt2ww
    from amber.amber10 import crd

    types='H-Rub_D20 D-Rub_H20'
    hs='h_20 h_37'
    Ts='T_25 T_50 T_75 T_100 T_125 T_150 T_175 T_200 T_225 T_250 T_275 T_300 T_325 T_350 '
        
    def convT(t,type='BASIS',FWHM=None,ss=None):
        """
        FWHM: full-width at half maximum in the energy domain, in micro-eV
        ss: standard deviation in the time-domain, in picoseconds"""
        from math import exp,pi,sqrt
        #time resolution windows, in picoseconds
        sigma={'dcs-nist':108.8,'original':49.0,'BASIS':521.0}
        if FWHM: 
            ss=FWHM/2.354820045   #in micro-eV
            ss=4136/(2*pi*ss)     #in picoseconds
        if not ss: ss=sigma[type]
        ss=2*ss*ss
        return exp(-t*t/ss)

    if task=='benchmark':

        """intermediate scattering function of the hydrogens in the
        first five residues (43 hydrogen atoms)"""
        
        qq=[1.0,]#Q-values
        dt=0.2   #time in between frames, in picoseconds
        tl=2000  #total time, in picoseconds
        nt=10000 #number of time points
        Nhydr=1500 #number of hydrogen atoms
        nb=43    #number of hydrogens for the benchmark
        tt=[]
        for i in range(nt): tt.append(i*dt)
        crdf=inp #trajectory file

        crdobj=crd(crdf,Nhydr)
        crdobj.loadframes(skip=0)

        scat=SIIScF(qq,crdobj.frames[:,0:nb],dt,tl,nt)
        pdb.set_trace()
        #store all intermediate factors
        os.system('mkdir -p isf')
        for i in range( nb ):
            storeSIIScF( 'isf/%04d.dat'%(i+1,), scat[i], qq, tt )
        storeSIIScF( 'isf.dat', scat.sum(axis=0)/len(scat), qq, tt )
        #store structure factors
        ww=[]; dw=4.136/(nt*dt)
        for i in range(nt): ww.append(i*dw)
        os.system('mkdir -p dsf')
        gaussianenvelope=[1.0]*len(tt)
        #gaussianenvelope=[]; factor=1.0/convT(tt[0])
        #for tti in tt: gaussianenvelope.append(factor*convT(tti))
        for iatom in range(nb):
            sc=scat[iatom]
            for iq in range(nq):
                sc[iq]=fourier(sc[iq],conv=gaussianenvelope).real
            storeSIIScF('dsf/%04d.dat'%(i+1,), scat[i], qq, ww)
        storeSIIScF('dsf.dat', scat.sum(axis=0)/len(scat), qq, ww)
        return False

    elif task=='isf':

        """intermediate scattering function of all the hydrogens"""
        
        qq=[1.0,]#Q-values
        dt=0.2   #time in between frames, in picoseconds
        tl=2000  #total time, in picoseconds
        nt=10000 #number of time points
        Nhydr=1500 #number of hydrogen atoms
        nb=43    #number of hydrogens for the benchmark
        tt=[]
        for i in range(nt): tt.append(i*dt)
        crdf=inp #trajectory file

        crdobj=crd(crdf,Nhydr)
        crdobj.loadframes(skip=0)

        scat=SIIScF(qq,crdobj.frames[:,0:nb],dt,tl,nt)
        pdb.set_trace()
        #store all intermediate factors
        os.system('mkdir -p isf')
        for i in range( nb ):
            storeSIIScF( 'isf/%04d.dat'%(i+1,), scat[i], qq, tt )
        storeSIIScF( 'isf.dat', scat.sum(axis=0)/len(scat), qq, tt )
        #store structure factors
        ww=[]; dw=4.136/(nt*dt)
        for i in range(nt): ww.append(i*dw)
        os.system('mkdir -p dsf')
        gaussianenvelope=[1.0]*len(tt)
        #gaussianenvelope=[]; factor=1.0/convT(tt[0])
        #for tti in tt: gaussianenvelope.append(factor*convT(tti))
        for iatom in range(nb):
            sc=scat[iatom]
            for iq in range(nq):
                sc[iq]=fourier(sc[iq],conv=gaussianenvelope).real
            storeSIIScF('dsf/%04d.dat'%(i+1,), scat[i], qq, ww)
        storeSIIScF('dsf.dat', scat.sum(axis=0)/len(scat), qq, ww)
        return False

    elif task=='average SF':

        """gather _A _B _C _D (or _1 _2 _3 _4) onto one """

        if inp in ('IESF','IESFv2'):
            I=inp
            batches = inp2
            batches = batches.split()
            N = len( batches )
            suf=inp3
            for type in types.split():
                for h in hs.split():
                    for T in Ts.split():
                        print '%s %s %s'%(type,h,T)
                        out='%s/%s/%s/%s'%(pjdout,type,h,T)
                        Id={};  buf='';  nc=0;
                        for batch in batches:
                            out2 = '%s/solv/%s%s.dat' %(out,I,batch)
                            for l in open(out2).readlines():
                                if l[0]=='#':
                                    if nc<2:
                                        buf += l; nc+=1
                                    continue
                                (Q,i) = l.split()
                                if Q not in Id.keys(): Id[Q]=float(i)
                                else: Id[Q] += float(i)
                        Qk=Id.keys()
                        Qk.sort()
                        for key in Qk:
                            buf += '%s %8.6f\n'%(key,Id[ key ]/N)
                        out2 = '%s/solv/%s%s.dat'%(out,I,suf)
                        open(out2,'w').write(buf)
        elif inp in ('IISF','IISF_FWHM'):
            from utilities.readingWritingFiles import read_to_cols as RC
            from utilities.readingWritingFiles import write_from_numpy
            I=inp
            batches = inp2
            batches = batches.split()
            N = len( batches )
            suf=inp3
            first=True
            for batch in batches:
                ifile = '%s%s.dat'%(I,batch)
                if first:
                    TQ,buf = RC(ifile,comment='#', separator=' ',
                                xtype='float', outFmt='numpy.array',
                                outComm=True)
                    first = False
                else:
                    tq = RC(ifile,comment='#', separator=' ',
                            xtype='float', outFmt='numpy.array',
                            outComm=False)
                    TQ += tq
            TQ /= len( batches )
            buf += write_from_numpy(None,TQ,format=' %7.5f',
                                    ncols=TQ.shape[1])

            open( '%s%s.dat'%(I,suf) ,'w').write(buf)

    elif task=='correct nan':

        a,tmpf=mkstemp(dir='/projects/tmp')
        os.system( 'grep -v nan %s > %s'%(inp,tmpf) )
        pt = open(tmpf)
        dt=0.2 #time in between frames
        buf = ''
        L1 = pt.readline()
        while L1[0] == '#':
            buf += L1
            L1 = pt.readline()
        buf += L1
        x1 = [ float(x) for x in L1.split() ]; tprev = x1[0]

        nQ = len(x1) -1   # number of Q values
        L2 = pt.readline()

        while L2:
            x2 = [ float(x) for x in L2.split() ]; tcurr = x2[0]
            Dt = tcurr - tprev

            if Dt > 1.5*dt:  # some lines missing (1.5 to
                             # avoid rounding errors)

                idt = dt
                #pdb.set_trace()                            
                while idt+dt/2 < Dt: # The dt/2 to avoid rounding errors
                    buf += ' %9.3f'%(tprev+idt,)
                    for iQ in range(nQ):
                        dQ = x2[1+iQ]-x1[1+iQ] #change in Q
                        sl = dQ / Dt        #slope in Q
                        Q =  x1[1+iQ] + sl * idt # linear extrapol
                        buf += ' %7.5f'%Q
                    buf += '\n'
                    idt += dt
            buf += L2
            L1 = L2
            x1 = deepcopy( x2 )
            tprev = tcurr
            L2 = pt.readline()
        #sys.stdout.write(buf)
        open( inp, 'w' ).write(buf) #overwrite file
        os.system('/bin/rm %s'%tmpf)
        
    elif task=='mix average SF':

        """mix the averages from the protein hydrogens(deuteriums) and
        the solvent deuteriums(hydrogens)"""

        if inp == 'IESFv2':
            from lib.utilities.readingWritingFiles import read_column
            I=inp
            nr=1500 #number of hydrog/deut in the four rubredoxins
            buf0= '#Incoherent Elastic Intermediate Structure Factor\n'
            buf0+='#Q_1(A^(-1)) IESF(Q)\n'
            for type in types.split():
                br=25.274 #incoh scat cross-sect of H in rubredoxin
                bw=4.04   #incoh scat cross-secti of deut in water
                if type=='D-Rub_H20':
                    br=4.04; bw=25.274
                for h in hs.split():
                    nw=528 #number of hydrog/deut in the water
                    if h=='h_37': nw=968
                    for T in Ts.split():
                        print '%s %s %s'%(type,h,T)
                        out='%s/%s/%s/%s'%(pjdout,type,h,T) #directory
                        inp='%s/solv/%s_Rub.dat' %(out,I)
                        Q=[float(x) for x in read_column(inp,1) ]
                        Ir=[float(x) for x in read_column(inp,2) ]
                        Ir=numpy.array( Ir )
                        inp='%s/solv/%s_wat.dat' %(out,I)
                        Iw=[float(x) for x in read_column(inp,2) ]
                        Iw=numpy.array( Iw )
                        Imix = (nr*br*Ir + nw*bw*Iw)/( nr*br + nw*bw)
                        buf=''
                        for i in range( len(Q) ):
                            buf+='%4.2f %8.6f\n'%(Q[i],Imix[i])
                        out2='%s/solv/%s_mix.dat'%(out,I)
                        open( out2, 'w' ).write( buf0+buf )
        elif inp in ('IISF','IISF_FWHM'):
            from utilities.readingWritingFiles import read_to_cols as RC
            from utilities.readingWritingFiles import write_from_numpy
            I=inp      #type of scattering function
            type=inp2  #H-Rub_D20 or D-Rub_H20
            h = inp3   #hydration (either h_20 or h_37)
            nr=1500    #number of hydrog/deut in the four rubredoxins
            br=25.274  #incoh scat cross-secti of hydrogen in rubredoxin
            bw=4.04    #incoh scat cross-secti of deut in water
            if type=='D-Rub_H20':
              br=4.04; bw=25.274 #swith values
            nw=528 #number of hydrog/deut in the water for h_20
            if '37' in h: nw=968
            #load the IISF for the Rubredoxins
            X,buf = RC( '%s_Rub.dat'%I, comment='#', separator=' ',
                        xtype='float', outFmt='numpy.array',outComm=True)
            X = X.transpose() # column index firs, row index second
            t = X[0]          #time column
            Qr = X[1:] #Q's columns
            #load the IISF for the waters
            X = RC( '%s_wat.dat'%I, comment='#', separator=' ',
                    xtype='float', outFmt='numpy.array', outComm=False )
            X = X.transpose()
            Qw = X[1:]
            #mix both IISF functions
            Qrw = (nr*br*Qr + nw*bw*Qw)/( nr*br + nw*bw)

            #output to file
            t = t.reshape(1,len(t))
            X= numpy.concatenate( (t,Qrw), axis=0 )
            nc = len(X)
            X = X.transpose()
            buf += write_from_numpy(None,X,format=' %7.5f', ncols=nc)

            open( '%s_mix.dat'%I, 'w' ).write( buf )    
        return False
    
    elif task=='remove resolution function':
        from utilities.readingWritingFiles import read_to_cols as RC
        from utilities.readingWritingFiles import write_from_numpy
        
        sf = inp      #file with intermediate scattering function
        resolf = inp2 #scattering at low T, works as resol funct
        outf = inp3
        X,buf = RC( resolf, comment='#', separator=' ',
                    xtype='float', outFmt='numpy.array',outComm=True)
        X = X.transpose() # column index firs, row index second
        t = X[0]          #time column
        R = X[1:]
        X = RC( sf, comment='#', separator=' ',
                xtype='float', outFmt='numpy.array',outComm=False)
        X = X.transpose()
        Q = X[1:] #each colum is a Q value
        Q = Q / R #divide by resolution function
        t = t.reshape(1,len(t))
        X= numpy.concatenate( (t,Q), axis=0 )
        nc = len(X)
        X = X.transpose()
        buf += write_from_numpy(None,X,format=' %7.5f', ncols=nc)
        open( outf, 'w' ).write( buf )
        return False
 
    elif task=='apply resolution function':
        from utilities.readingWritingFiles import read_to_cols as RC
        from utilities.readingWritingFiles import write_from_numpy
        sf = inp      #file with intermediate scattering function
        FWHM=float(inp2)
        X,buf = RC( sf, comment='#', separator=' ',
                xtype='float', outFmt='numpy.array',outComm=True)
        X = X.transpose()
        T = X[0]  #list of times
        R = numpy.array( [ convT(t,FWHM=FWHM) for t in T.flat ] ) #resolution function
        Q = X[1:] #each colum is a Q value
        Q = Q * R #apply the resolution function
        T = T.reshape(1,len(T))
        X= numpy.concatenate( (T,Q), axis=0 )
        nc = len(X)
        X = X.transpose()
        buf += write_from_numpy(None,X,format=' %7.5f', ncols=nc)
        open( out, 'w' ).write( buf )
        return False
   
    elif task=='IISF fourier transform':
        from utilities.readingWritingFiles import read_to_cols as RC
        from utilities.readingWritingFiles import write_from_numpy

        sf = inp      #file with intermediate scattering function
        X,buf = RC( sf, comment='#', separator=' ',xtype='float', outFmt='numpy.array',outComm=True)
        X = X.transpose()
        

        #modify header
        y = buf.strip().split('\n')
        buf = '\n'.join( [ y[0],'#E(micro-eV) Q_1(A^(-1)) Q_2 ...',y[-1] ] ) + '\n'


        #create list of energies centered at zero
        T = X[0]  #list of times. We assume t=0 is the first entry
        dt = T[1]-T[0]
        nt = len( T )
        dw = 4.136/(nt*dt)
        W = dw * numpy.arange(nt)
        W = numpy.concatenate([-W[::-1][:-1],W])
        nW = len(W)
        W = W.reshape(1,nW)


        #create the fourier transform
        Q = X[1:] #each colum is a Q value
        S = []
        for q in Q: 

            L = numpy.concatenate([q[::-1][:-1],q])
            Lx = fourierRealSymmetric( numpy.concatenate([q[::-1][:-1],q]))
            #mirror image for of IISF for negative values of time
            S.append( fourierRealSymmetric( numpy.concatenate([q[::-1][:-1],q]) ) )

        S = numpy.array( S )

        #output to file
        Y= numpy.concatenate( (W,S), axis=0 )
        nc = len(Y)
        Y = Y.transpose()

        buf += write_from_numpy(None,Y,format=' %10.6f', ncols=nc)
        open( out, 'w' ).write( buf )
        return False

    elif task=='SF versus T':

        """ create files of structure factor versus temperature"""

        I = inp  #I=IESF_Rub
        for type in types.split():
            for h in hs.split():
                Id = {}
                print '%s %s'%(type,h)
                for T in Ts.split():
                    Id[ T ] = {}

                    out='%s/%s/%s/%s/solv/%s.dat'%(pjdout,type,h,T,I)
                    for l in open( out ).readlines():
                        if l[0] == '#': continue
                        (Q,i) = l.split()
                        Id[ T ][ Q ] = i
                buf='#Incoherent Elastic Intermediate Structure Factor'
                kT = Ts.split()
                buf+='#   T   '
                x = Id[ kT[0] ]
                kQ=x.keys()
                kQ.sort()
                for Q in kQ: buf += ' %s'%Q
                buf += '\n'
                for T in kT:
                    buf+=' %3s'%T[2:]
                    x = Id[ T ]
                    kQ=x.keys()
                    kQ.sort()
                    for Q in kQ: buf += ' %s'%x[ Q ]
                    buf+='\n'
                out = '%s/%s/%s/analysis/%s.dat'%(pjdout,type,h,I)
                open( out, 'w' ).write( buf )
    elif task=='Vineyard approx':
        """evaluate the Vineyard approximation for water:
        """
        import MDAnalysis as mda
        from scatter.scatterMDA import Qlist,getFrames,scatL,IC,II,SSF
        from utilities.readingWritingFiles import write_from_numpy
        
        group = inp
        sfx = out #suffix to append to output file
        
        #create the Universe
        os.system('/bin/rm junk.prmtop junk.dcd 2> /dev/null')
        os.system('ln -s top junk.prmtop; ln -s dcd junk.dcd')
        U = mda.Universe('junk.prmtop','junk.dcd')

        QQ=Qlist(0.1, 2.0, 0.1)
        #print 'calculate the intermediate coherent scattering function'
        b=scatL(U,group)['coherent'] #coherent scattering lengths
        frames = getFrames(U,group,skip=1,format='fac')
        sf = IC(frames,b,QQ,nsampling=50,nt=100,dt=50)
        open( 'ic_%s_MDA.dat'%sfx, 'w').write(sf['buf'])

        print 'compute the static structure factor suitable for the'
        print 'Vineyard approximation'
        sf = SSF(frames,b,QQ,nsampling=200)
        open( 'ssf_%s_MDA.dat'%sfx, 'w').write(sf['buf'])
    
        print 'compute the incoherent structure factor suitable for the'
        print 'indistinguishable particle approximation'
        print 'Scoh ~ SSF - (sum_i bi^2 - II)'
        sf = II(frames,b*b,QQ,nsampling=200,nt=100,dt=50)
        ts = sf['ts']  #time values
        sfs = sf['sf'] #scattering function
        sfs = (b*b).sum()/(len(b)*len(b)) - sfs/len(b)
        buf = ''; nQ=len(QQ)
        for it in range( len(ts) ):
            s = write_from_numpy(None,sfs[it],ncols=nQ,format=' %10.7f')
            buf+=' %6.3f %s'%(ts[it],s)
        open( 'ii.2_%s_MDA.dat'%sfx, 'w').write(buf)

        print 'compute the incoherent structure factor suitable for the'
        print 'Vineyard approximation'
        print 'Scoh ~ SSF * II'
        b = numpy.ones( len(b) ) #no scattering lenghts go here
        sf = II(frames,b,QQ,nsampling=200,nt=100,dt=50)
        open( 'ii_%s_MDA.dat'%sfx, 'w').write(sf['buf'])       

        
        os.system('/bin/rm junk.prmtop junk.dcd 2> /dev/null')
    elif task=='output distances':
        import MDAnalysis as mda
        from scatter.scatterMDA import getFrames
        from MDAnalysis.analysis.distances import distance_array
        group = 'resname WAT'
        #create the Universe
        os.system('/bin/rm junk.prmtop junk.dcd 2> /dev/null')
        os.system('ln -s top junk.prmtop; ln -s dcd junk.dcd')
        U = mda.Universe('junk.prmtop','junk.dcd')
        QQ=[1.0,]
        frames = getFrames(U,group,skip=1,format='fac')
        (nfr,nat,nc) = frames.shape
        t = int( nfr/4 ) #pick a time t
        #iterate through all atoms
        t0 = int(nfr/2)  #pick a t0
        fi = frames[t0]; fj = frames[t0+t]
        dii=distance_array(fi,fj,None).diagonal()
        buf='\n'.join( [str(x) for x in dii] )
        open('junk.dat','w').write(buf)
        #iterate through as many as t0's as atoms
        i = 600        #pick one atom
        buf=''
        delta = (nfr-t)/nat; t0=0
        for it in range(nat):
            fij = frames[t0][i]-frames[t0+t][i] #particular atom
            buf+=str( numpy.sqrt((fij*fij).sum()) )+'\n'
            t0 += delta
        open('junk.2.dat','w').write(buf)

    elif task=='check static structure factor':
        import MDAnalysis as mda
        from scatter.scatterMDA import getFrames,Qlist,SSFcheck,scatL
        outputfile = out #write the structure factor to this file
        group = 'resname WAT'
        #create the Universe
        os.system('ln -s top junk.prmtop; ln -s dcd junk.dcd')
        U = mda.Universe('junk.prmtop','junk.dcd')
        QQ=Qlist(0.1, 2.0, 0.1)
        b=scatL(U,group)['coherent'] #coherent scattering lengths
        frames = getFrames(U,group,skip=1,format='fac')
        sf = SSFcheck(frames,b,QQ,nsampling=200)
        open(outputfile,'w').write(sf['buf'])

    elif task=='generate dif_sphere_V1':
        from utilities.readingWritingFiles import read_to_cols as RC

        #load the SSF function
        X = RC( inp, comment='#', separator=' ', xtype='float', outFmt='numpy.array',outComm=False)
        #Ken used Q values from 0.3 to 1.9, every 0.2. This corresponds to columns 3 to 11
        X = X.transpose() # column index firs, row index second

        #Data goes is in the [-20,20] micro-eV, but we will read only
        #between -4 micro-eV and +4 micro-eV because there are too
        #many points otherwise

        W = X[0]
        n = len(W); d=(n-1)/40; nh=(n-1)/2
        W = W[nh-d:nh+d+1]
        nW = len(W)

        Q = X[2:11] #structure factors for the different Q values
        nQ = len(Q) 

        buf ='# Number of energy transfer values\n%d\n'%len(W)
        buf+='# Number of Q transfer values\n%d\n'%nQ
        buf+='# energy transfer (micro-eV) Values:\n'
        buf+='\n'.join( ['%6.3f'%w for w in W.flat] )+'\n'
        buf+='# Q transfer (1/Angstroms) Values:\n'
        buf+='\n'.join( ['%4.2f'%q for q in (numpy.arange(0.3,2.1,0.2)).flat ] )+'\n'

        for iq in range(nQ):
            buf+='# Group n\n'
            for iw in range(nW):
                q = Q[iq]
                q = q[nh-d:nh+d+1] #add for negative values of the energy
                q = 0.9387 + 0.044787 * q #rescale to mimick values from the Rubrdoxin boson peak study
                buf+='%7.5f %7.5f\n'%(q[iw], 0.1*numpy.sqrt(q[iw]) ) #Ken estimated error as 0.1*sqrt(value)

        open(out,'w').write(buf)

    return False

def hydration(inp=None,inp2=None,inp3=None,out=None,out2=None,
              task=None):
    from utilities.readingWritingFiles import read_column as readcol
    from utilities.readingWritingFiles import read_to_cols as RC
    from utilities.readingWritingFiles import write_from_numpy

    if task=='watershell':
        buf = '#residue  <first-shell> std  <second-shell> std\n'
        nres = 54 #number of residues + FE in rubredoxin
        for ires in range(1,1+nres):
            x1=[] #first shell
            x2=[] #second shell
            for iprot in range(0,4):  #average for all four rubredoxins
                jres = ires + nres*iprot
                inp='watershell/%03d.dat'%jres
                y = readcol(inp,2,isInt=1)
                x1.append( y )
                y = readcol(inp,3,isInt=1)
                x2.append( y )
            x1 = numpy.array( x1 )
            x2 = numpy.array( x2 )
            buf +=' %3d %5.2f %5.2f %5.2f %5.2f\n'%(ires,x1.mean(),x1.std(),x2.mean(),x2.std())
        open('watershell.dat','w').write(buf)
        
    elif task=='interprotein contacts':

        buf = '# Residue pair contact number\n'
        nres = 54 #number of residues + FE in rubredoxin + 1
        map = numpy.zeros(1+nres)
        for line in open('residue_map.dat'):
            if line[0]=='#': continue
            items = line.split( )
            i = int( items[0] )
            j = int( items[1] )
            if abs(i-j) >= nres: #interprotein contact
                if i>nres:
                    i = i - (i/nres)*nres
                    if i==0: i=nres
                if j>nres:
                    j -= (j/nres)*nres
                    if j==0: j=nres
                map[i] += float( items[2] )
                map[j] += float( items[2] )
        for i in range(1,1+nres):
            buf+=' %3d %6.3f\n'%(i,map[i])
        open('residue_map.v2.dat','w').write(buf)
        
    elif task=='dielectric':
        """calculate the dielectric constant"""

        dipf = 'dipole_WAT.dat' #file containing dipole
        T = float(inp)   #temperature
        kT = 1.99E-03*T #Kcal/mol
        V = [float(x) for x in inp2.split()] #Vx Vy Vz
        V = V[0] * V[1] * V[2] #volume of the box
        dt = 0.2 #frame between numbers
        
        X = RC( dipf, comment='#', separator=' ',
                xtype='float', outFmt='numpy.array',outComm=False)
        X=X.transpose()
        
        t = X[0] * dt #time
        M= X[1:4] 
        M = M.transpose() #dipole moments of the solvent

        Mav = numpy.average(M,axis=0)  #dipole averaged over simulation
        M2  = numpy.average( M*M, axis=0)
        DM  = (M2 - Mav*Mav).sum() #fluctuation of the dipole moment
        e = 1+4.189*DM/(kT*V) #dielectric constant
        
    elif task=='water_diffusion':
        """calculate the diffussion constant"""
        from scipy import polyfit
        diff = 'water_diffusion_r.xmgr'
        X = RC( diff, comment='#', separator=' ',
                xtype='float', outFmt='numpy.array',outComm=False)
        X=X.transpose()
        t = X[0] # time
        MSD = X[1] #average diffussion
        (D,b) = polyfit(t,MSD,1) #linear regression
        if D<=0: D=1E-06  #arbitrary low limit
        open('water_diff_coeff.dat','w').write(`D`+'\n')
        
    elif task=='average shortest distance':
        """calculate average shortest distance between each hydrogen
        atom in the protein and any water oxygen. The average is taken
        over a number of frames."""

        import MDAnalysis as mda
        from myMDAnalysis.utils import shortest_distance

        #select 200 hydrogens by default
        N = 200; suffix = ''
        if inp:
            N = int(inp); suffix = '_'+inp

        #create the Universe
        os.system('/bin/rm junk.prmtop junk.dcd 2> /dev/null')
        os.system('ln -s top junk.prmtop; ln -s dcd junk.dcd')
        U = mda.Universe('junk.prmtop','junk.dcd')
       
        group1='resid 1-216 and name H*'
        group2='resname WAT and name O'
        
        #calculate shortest distance averaged over considered frames
        sd = shortest_distance(U,group1,group2,skip=100)
        sd = numpy.average(sd,axis=0)

        #sort from minimum to maximum distance
        perm = numpy.argsort( sd )

        X=U.selectAtoms(group1)

        #find number of hydrogens belonging to methyl groups for the
        #first N hydrogens and the last N hydrogens
        from seq.letters import isMethyl
        n=0; n1=0; n2=0;
        for i in range( len(X) ):
            at = X[perm][i]
            if isMethyl(at.resname,at.name):
                n  += 1
                if i < N: n1 += 1
            at = X[perm][-i]
            if i<N and isMethyl(at.resname,at.name): n2 += 1
        buf  = '#number of methyls in %d hydrogens = %d\n'%(len(sd),n)
        buf += '#number of methyls in first %d hydrogens = %d\n'%(N,n1)
        buf += '#number of methyls in last  %d hydrogens = %d'%(N,n2)
        open('n_methyls%s.dat'%suffix,'w').write(buf)

        #find number of methyl groups versus number of hydrogens as we
        #increase the distance to the solvent
        n=0
        buf = '#H-index number-of-methyls dist-to-solv\n'
        for i in range( len(X) ):
            at = X[perm][i]
            if isMethyl(at.resname,at.name): n += 1
            buf += '%4d %6.3f %3d\n'%(1+i,sd[perm][i],n)
        open('methyl_vs_all%s.dat'%suffix,'w').write(buf)
                
        #write to file
        buf ='#(1)H-index (2)pdb-number (3)averag-shortest-dist-to-solv\n'
        buf+='#(4)charge (5)name (6)residue\n'
        for i in range(len(sd)):
            at = X[i]
            buf += '%4d %4d %5.2lf %5.2f %4s %3s\n'%\
                   (i,1+at.number,sd[i],at.charge,at.name,at.resname)
        open('asd_H_solv%s.dat'%suffix,'w').write(buf)
        
        #write to file sorted from minimum to maximum distance
        buf ='#(1)H-index (2)pdb-number (3)averag-shortest-dist-to-solv\n'
        buf+='#(4)charge (5)name (6)residue\n'
        for i in range(len(sd)):
            at = X[perm][i]
            buf += '%4d %4d %5.2lf %5.2f %4s %3s\n'%\
                   (i,1+at.number,sd[perm][i],at.charge,at.name,at.resname)
        open('asd_H_solv_sorted%s.dat'%suffix,'w').write(buf)
    
        #calculate diffussion with msd2 function
        from amber.amber10 import extricateBC
        from scatter.scatter import msd2
        from scatter.scatterMDA import getFrames
        delta = 0.2    #time between frames
        tl = 1000.0    #compute msd2 from t=0 to t=tl
        nt = int(tl/delta)  #compute msd2 nt time in between t=0 and t=tl

        #diffussion of protein hydrogens
        rr = getFrames(U,group1,skip=1,format='afc')
        #for the first N atoms
        xyz = rr[perm][0:N].transpose( (1,0,2) )  #format fac
        extricateBC(xyz,U.coord.dimensions)
        diff=numpy.average( msd2(xyz,delta,tl,nt,norig=100), axis=0 )
        buf='#time MSD2\n'
        for it in range(nt):
            buf += '%6.2f %6.3f\n'%(delta*it,diff[it])
        open('diff_sda_closest%s.dat'%suffix,'w').write(buf)
        
        #for the last N atoms
        xyz = rr[perm][-N:].transpose( (1,0,2) )  #format fac
        extricateBC(xyz,U.coord.dimensions)
        diff=numpy.average( msd2(xyz,delta,tl,nt,norig=100), axis=0 )
        buf='#time MSD2\n'
        for it in range(nt):
            buf += '%6.2f %6.3f\n'%(delta*it,diff[it])
        open('diff_sda_furthest%s.dat'%suffix,'w').write(buf)

        #diffusion of water oxygens
        rr = getFrames(U,group2,skip=1,format='fac')
        xyz = rr
        extricateBC(xyz,U.coord.dimensions)
        diff=numpy.average( msd2(xyz,delta,tl,nt,norig=100), axis=0 )
        buf='#time MSD2\n'
        for it in range(nt):
            buf += '%6.2f %6.3f\n'%(delta*it,diff[it])
        open('diff_WAT_O.dat','w').write(buf)
        os.system('/bin/rm junk.prmtop junk.dcd')

    elif task=='residence times':
        """list all hydrogen--oxygen residence times according to
        function residence_times"""
        from myMDAnalysis.utils import residence_times
        import MDAnalysis as mda
        
        #create the Universe
        os.system('/bin/rm junk.prmtop junk.dcd 2> /dev/null')
        os.system('ln -s top junk.prmtop; ln -s dcd junk.dcd')
        U = mda.Universe('junk.prmtop','junk.dcd')
       
        group1='resid 1-216 and name H*'
        group2='resname WAT and name O'

        ico=3.0  #boundary for contact volume
        oco=4.0  #boundary for no-man land volume
        dt = 0.2 #time in picoseconds between consecutive frames
        rt=residence_times(U,group1,group2,ico,oco,skip=1,dt=dt)
        buf='# time RES name numb charge rest name num charge\n'
        for info in rt:
            info = tuple( info )
            #t = x[0]
            #Hres = x[1]; Hname = x[2]; Hnum = x[3]; Hq = x[4]
            #Ores = x[5]; Oname = x[6]; Onum = x[7]; Oq = x[8]
            buf +='%6.1f %3s %4s %5d %5.2f %3s %4s %5d %5.2f\n'%info

        open('residence_times.dat','w').write(buf)
    return False

def namd(task=None,inp=None,out=None,group=None):

    if task=='prepare input files':
        labels={'prot':range(1,217),
                'ions':range(217,253),
                'water':range(253,10000)+['****',],
                'solv':range(217,10000)+['****',],
                'all':range(1,100000)+['****',],
                }
        intra = False
        a,b=group.split('-')
        if a == b: intra = True
        if a in labels.keys(): la = labels[a]
        else: la = [ int( a ), ] #just a single residue number
        if b in labels.keys(): lb = labels[b]
        else: lb = [ int( b ), ]
        #create pdb files
        buf = ''; count=0
        for l in open('namd/pdb').readlines():
            if l[:4]=='ATOM':
                ires = l[22:26]
                try:
                    ires = int(ires)
                except:
                    pass
                if ires in la:
                    l = l[:62]+'1.00\n'
                if ires in lb:
                    if intra:
                        l = l[:62]+'1.00\n'
                    elif ires not in la:
                        l = l[:62]+'2.00\n'
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
    elif task=='parse energies from log file':
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
        open('namd/%s_energies.dat'%group,'w').write(buf2+buf)
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
