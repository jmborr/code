#!/usr/bin/python

import os,sys,re,numpy,pdb,argparse
from tempfile import mkstemp,mkdtemp
from utilities.can import CAN,CANl
from utilities.codedir import codedir,scratchdir,projectsdir
from utilities.small_utilities import chomp,Bye

"""project description goes here"""

#global variables
pwrite=sys.stdout.write
joink=os.path.join; jk=joink
isthere=os.path.exists
basename=os.path.basename
 
my_project='RdCp_RdPf'
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

    def qh(self,conditions='T_300 10ns unsolv',task='init',
           inpf=None,inpf2=None,inpf3=None,inpf4=None,
           outf=None,outf2=None,outf3=None):
        from amber.amber10 import mwcovar,crd,evecs,evecs_projection
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        #find out temperature
        T=float(re.compile('T_(\d+)').search(conditions).group(1))
        if task=='init':
            from amber.amber10 import evecs
            pt['evecs']=evecs(jk(sd,'evecs.dat'))
            return True
        elif task=='make histograms':
            """make histograms of the mode projections"""
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
        elif task=='mode std':
            """read a mode projection file and output the average and
            standard deviation and the cummulative standard deviation"""
            av=[]; std2=[]; il=0; N=0; Nstd2=0.0
            for l in open(inpf).readlines()[2:]:
                Q2=0.0
                lst=[float(x) for x in l.split()]
                N=len(lst)
                if len(av)==0:
                    av=[0.0]*N; std2=[0.0]*N
                for i in range(1,N):
                    q2=lst[i]*lst[i]
                    Q2+=q2
                    av[i]+=lst[i]; std2[i]+=q2                    
                il+=1; #print Q2
            #output to file
            buf='#mode (std-dev)^2 average cumm-(std-dev)^2\n'
            for i in range(1,N):
                av[i] = av[i]/il
                std2[i]= std2[i]/il - av[i]*av[i]
                Nstd2+=std2[i]
                buf+='%5d %6.3f %8.5f %7.3f\n'%(i,std2[i],av[i],Nstd2)
            open(outf,'w').write(buf)
            print "total std^2=%f\n"%Nstd2
            return False
        elif task=='check CA evecs':
            from math import sqrt
            obj=evecs(inpf); obj.loadModes()
            obj2=obj
            if inpf2:
                obj2=evecs(inpf2); obj2.loadModes()
            #
            ##rotate modes according to rotation R
            R=numpy.array( 0.26560286,  0.76209860,  0.59047510,
                          -0.92655922,  0.37099580, -0.06204935,
                          -0.26635151, -0.53062966,  0.80466703 )
            R=R.reshape(3,3)
            obj.rotate_modes(R)
            #
            ##dot product between the modes
            #m=12.01 #mass of CA
            #for I in range(obj.nq):
            #    for J in range(obj2.nq):
            #        s=numpy.dot(obj.qs[I].m.ravel(),obj2.qs[J].m.ravel())
            #        print "I=%3d J=%3d s=%f"%(I,J,m*m*s*s)
            return False
        elif task=='check CA fluctuations':
            obj=crd(inpf,52)
            avf=obj.avframe()
            obj.loadframes()
            for f in obj.frames:
                f = f - avf
                print (f.ravel()*f.ravel()).sum()
            return False
        elif task=='CA modes similarity':
            #load evecs files
            eA=evecs(inpf);  eA.loadModes()
            eB=evecs(inpf2); eB.loadModes()
            #list of masses, all are carbons
            M = [ 12.01 ] * 3*eA.natoms
            #calculate similarity
            D,H,buf = eA.HessSimilarity(eB,M,M,RMS=True)
            open(outf,'w').write(buf)
            #output the absolute value of the dot product
            buf='#absolute value of dot product\n'
            for i in range(len(M)):
                for j in range(len(M)):
                    buf += '%4d %4d %6.4f\n'%(1+i,1+j,D[i][j])
            open(outf2,'w').write(buf)
            return False
        elif task=='mode align':
            #load evecs files
            eA=evecs(inpf);  eA.loadModes()
            eB=evecs(inpf2); eB.loadModes()
            #align modes of the first to the second
            eA.align_modes_to(eB)
            #overwrite average structure
            eA.avgX=eB.avgX
            #print to file
            open(outf,'w').write( eA.write2string() )
            return False
        elif task=='output global fluctuation':
            """output total fluctuations from mwcovar """
            m1 = mwcovar( inpf )  
            print m1.total_fluct()
            return False
        elif task=='output eigenvalues':
            """output eigenvalues of mwcovar """
            m1 = mwcovar( inpf )  
            m1.findModes()
            buf= '\n'.join([str(x) for x in m1.evals])+'\n&'; sys.stdout.write( buf )
            return False
        elif task=='rotate mwcovar':
            """given matrix C, transform by R^t * C * R where R is a
            rotation matrix from superposing average conformation
            stored in inpf to average conformation stored in inpf2"""
            from tertiary.rmsd import rmsd
            from scipy.linalg import block_diag
            evAx = evecs(inpf).avgX   #average structure of inpf
            evBx = evecs(inpf2).avgX
            N = evAx.shape[0] #number of particles
            #R (or U) takes coords in the evA frame of reference to the evB frame of reference
            r, R, t = rmsd( evAx, evBx ) #evAx --> R * (evAx-CMA) + CMB
            L = [R,]*N           #create a list of R repeated N times
            U = block_diag( *L ) #generate block diagonal matrix. Treat L as a list of 3x3 matrixes
            """assuming covar "C" for RdCp, U^t * C * U is covar matrix for RdCp in the frame of reference of RdPf
            x = U * y, where "y" are coords in the RdPf frame of reference and "x" are coords in the RdCp frame
            x^t * C^-1 * x = y^t * (U^t * C^-1 * U ) * y
            U^t * C^-1 * U is covar matrix for RdCp but in the RdPf frame of reference"""
            m = mwcovar( inpf3 ) 
            m.mwc = numpy.dot(U.transpose(), numpy.dot( m.mwc, U ) )
            m.sav2fil( outf, format='full precision' )
            return False
        elif task=='renormalize mwcovar':
            """renormalize and save mwcovar"""
            m1 = mwcovar( inpf )
            m1.mwc *= float( inpf2 )
            m1.sav2fil( outf, format='full precision' )
            return False
        elif task=='substract correlation matrixes':
            """note that we may have to renormalize the particle
            fluctuations of m2 by those of m1"""
            m1 = mwcovar( inpf )
            m2 = mwcovar( inpf2 )
            m1.mwc = m1.mwc - m2.mwc * (m1.total_fluct()/m2.total_fluct())
            m1.sav2fil( outf, format='full precision' )
            print 'saving inpf - inpf2 to inpf3'
            m1.sav2fil( outf2, format='as 3 columns, isotropic average' )
            print 'saving inpf - inpf2 as three columns in outf2'
            m1.sav2fil( outf3, format='isotropic average' )
            print 'saving inpf - inpf2 as isotropic average in outf3'
            return False
        elif task=='filter out high frequency modes':
            """from the covariance matrix, keep only the modes with
            the lowest frequency giving rise to a certain fraction of
            the total mean square fluctuations """
            dot  = numpy.dot
            m1 = mwcovar( inpf )
            co = float( inpf2 ) #filter cutoff
            m1.findModes()
            L = m1.evals  #L: list eigenv, biggest to smallest.
            E = m1.evecs  #E: rotation matrix
            MSF = L.sum()
            msfco = co * MSF  #cut-off mean-square-fluctuations
            i = 0; msf = 0.0
            while msf < msfco: 
                msf += L[i]; i += 1
            L[i:] = 0.0 # remove rest of values
            m1.mwc = dot(E, dot( numpy.diag(L),E.transpose() ) )
            m1.sav2fil(outf,format='full precision')
            print 'MSF=%f msf=%f'%( MSF, msf)
            return False
        elif task=='difference approximation':
            """note that we renormalize the particle
            fluctuations of m2 by those of m1"""
            pinv = numpy.linalg.pinv #pseudoinverse
            dot  = numpy.dot
            #did we pass cutoffs?
            co_1=None
            if inpf3: co_1=float(inpf3)
            co_2=None
            if inpf4: co_2=float(inpf4)

            #Find the pseudoinverse of correlation matrix
            m1 = mwcovar( inpf )
            #D1 diagonalized m1, K1 is a 1D-array of the inverse of
            #the eigenvalues. Set all eigenvalues of m1 smaller than
            #co to zero.

            K1, D1 = m1.pseudoinverse(co=co_1)
            #buf= '\n'.join([str(x) for x in K1])+'\n&'; open('junk_K1.dat','w').write(buf)
            K1 = numpy.diag( K1 ) #make a diagonal matrix using the 1D-array list of eigenvalues
            deltaK  = dot( D1.transpose(), dot( K1, D1) ) #pseudoinverse of m1
            
            #Find the pseudoinverse of correlation matrix
            m2 = mwcovar( inpf2 )
            K2, D2 = m2.pseudoinverse(co=co_1)
            #buf= '\n'.join([str(x) for x in K2])+'\n&'; open('junk_K2.dat','w').write(buf)
            K2 = numpy.diag( K2 )
            deltaK -= dot( D2.transpose(), dot( K2, D2) )

            #substract the two pseudoinverse - Force constant matrix
            mK = mwcovar(deltaK)
            dK, E = numpy.linalg.eigh( deltaK ) #diagonalize deltaK
            #buf= '\n'.join([str(x) for x in dK])+'\n&'; open('junk_dK.dat','w').write(buf)

            #invert the force constant matrix - correlation changes matrix 

            dC, D3 = mK.pseudoinverse(co=co_2, tZero=False)

            buf= '\n'.join([str(x) for x in dC])+'\n&'; open('junk_dC.dat','w').write(buf)

            dC = dot( D3.transpose(), dot( numpy.diag(dC), D3) )
            m3 = mwcovar( dC )

            m3.sav2fil( outf, format='full precision' )
            m3.sav2fil( outf2, format='as 3 columns, isotropic average' )
            
            return False
        elif task=='force constant':
            """note that we renormalize the particle
            fluctuations of m2 by those of m1"""
            m1 = mwcovar( inpf )
            #find first contact map of native file
            co = numpy.zeros( m1.n * m1.n ).reshape( m1.n, m1.n )
            for l in open(inpf2).readlines():
                if l[0]=='#':continue
                resi, resj, x = l.split( )
                resi = int( resi ) -1 ; resj = int( resj ) -1
                if resi > 51 or resj > 51: continue
                if float( x ) < 0.75: continue #no native contact
                for i in range( 3 ):
                    for j in range( 3 ):
                        co[ 3*resi+i ][ 3*resj+j ] = 1.0
            mlist = [ 12.01 ] * (m1.n/3)
            m1.mwc = m1.force_constant_matrix(T=T, mlist=mlist, co_map=co)
            #save to file
            m1.sav2fil(outf,format='as 3 columns, isotropic average')
            print 'saving force constant matrix in outf3'
            return False
        elif task=='output mwcovar as 3 column - isotropic average':
            mwc = mwcovar( inpf )
            mwc.sav2fil( inpf2, format='as 3 columns, isotropic average' )
            return False
        elif task=='find modes from mwcovar':
            mwc = mwcovar( inpf )           
            natoms = mwc.n/3
            M = [ 12.01 ] * natoms
            buf = mwc.saveModes( outf, eigFormat='frequency',
                                 atomic_masses = M,
                                 format='ptraj evecs' )
            if not buf: return False
            buf = '#mode-number eigenvalue [ a.m.u * Angstrom^2 ]\n'
            for i in range( len( mwc.evals) ):
                buf += '%3d %8.3f\n' %( 1+i, mwc.evals[ i ] )
            sys.stdout.write(buf)
            return False
        elif task=='make pdbs':
            """create a pdb for each mode, where we assign a bfactor
            to each residue based on the absolute value of the mode
            component, averaged over the X,Y,Z directions"""
            if not isthere( outf ): os.system('mkdir -p '+outf)
            ev = evecs( inpf )
            masses = [12.01] * ev.natoms
            ev.massProcess( masses )
            ev.create_PDBs( inpf2, outdir=outf, byres=True )
            return False
        elif task=='make pdb movies':
            """create a pdb movie for each mode"""
            ev = evecs( inpf )
            masses = [12.01] * ev.natoms
            ev.massProcess( masses )
            ev.gen_movies( inpf2, outdir=outf, byres=True )
            return False

        
    def namd(self,conditions='T_340',task=None,group=None,
             inp=None, inp2=None, out=None):
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        if task=='create box.xsc':
            tpl=open( '%s/%s/namd/box.xsc'%(pjdout,self.id) ).read()
            l=os.popen('tail -1 %s/solv/first.rst'%sd).read()
            b=[float(x) for x in l.split()][0:3]
            d={'A_X':b[0], 'B_Y':b[1], 'C_Z':b[2],
               'O_X':b[0]/2, 'O_Y':b[1]/2, 'O_Z':b[2]/2}
            for (k,v) in d.items(): tpl=tpl.replace(k,str(v))
            open('%s/solv/namd/box.xsc'%sd,'w').write(tpl)
            return False
        elif task=='energy evaluation':
            wd=mkdtemp(dir='/projects/tmp')
            print 'Temporary directory: ',wd
            nd='%s/%s/namd'%(pjdout,self.id)
            #move input files, then run namd
            cmd='cd %s/solv && cp top %s && cp %s/solv/namd/box.xsc '\
                 '%s && cd %s && cp %s.namd %s.pdb %s'\
                 %(sd,wd,sd,wd,nd,group,group,wd)
            os.system(cmd)
            cmd='cd %s; ln -s %s/solv/dcd dcd; namd2 %s.namd > %s.log && '\
                 'mv %s.log %s/solv/namd'%(wd,sd,group,group,group,sd)
            os.system(cmd)
            os.system('rm -rf '+wd) #remove working directory
            return False
        elif task=='generate mturn input by residue':
            """create *.namd and *.pdb input files for each residue of
            the multi-turn region. We will evaluate the energy/force
            cause on these residue by water and by the protein
            residues not-belonging to the multi-turn region """
            start, end = [ int(x) for x in inp.split('-') ]
            partner = group.split('-')[-1]
            #create %02d-****.pdb and %02d-****.namd input files
            lines = open( '%s/%s.pdb'%(sd,group) ).readlines()
            buf=''
            for res in range(start,1+end):
                for l in lines:
                    if l[0:4] == 'ATOM':
                        rid = int(l[21:26]) #residue number
                        if rid >= start and rid <= end and rid != res:
                            l = l.replace('    1.00','    0.00')
                    buf += l
                open( '%s/%02d-%s.pdb'%(sd,res,partner), 'w' ).write( buf )
                buf = open('%s/%s.namd'%(sd,group) ).read()
                buf = buf.replace( group, '%02d-%s'%(res,partner) )
                open( '%s/%02d-%s.namd'%(sd,res,partner),'w').write( buf )
            return False
        elif task=='parse log file':
            buf='#    N  electr   VdW            f_electr(4-6)'\
                 '            f_VdW(7-9)\n'
            for l in open('%s/solv/namd/%s.log'%(sd,group)).readlines():
                if l[0:6]=='ENERGY':
                    items=l.split()
                    N=float( items[1] )
                    ele,vdw=[ float(x) for x in items[6:8] ]
                    buf+=' %5d %8.2f %6.2f'%(1+N,ele,vdw)
                elif l[0:16]=='PAIR INTERACTION':
                    items=l.split()
                    fvdw=[ float(x) for x in items[5:8] ]
                    fele=[ float(x) for x in items[9: ] ]
                    buf+=' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n'\
                          %tuple(fele+fvdw)
            open('%s/solv/namd/%s.dat'%(sd,group),'w').write(buf)
            return False
        elif task=='parse log file II':
            buf  = '#    N   E-bonded E-nonbonded       E      electr'
            buf += '       VdW\n'
            avEb=0.0; avEnb=0.0; avE=0.0
            avELECT =0.0; avVDW=0.0
            for l in open('%s/solv/namd/%s.log'%(sd,group)).readlines():
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
                    buf+=' %5d %10.2f %10.2f %10.2f %10.2f %10.2f\n'\
                          %(1+TS,eb,enb,e,ELECT,VDW)
            TS+=1; print TS
            avEb/=TS; avEnb/=TS; avE/=TS; avELECT/=TS; avVDW/=TS;
            buf2='#     avEb      avEnb       avE      avELECT     avVDW\n'
            buf2+='#%10.2f %10.2f %10.2f %10.2f %10.2f\n'\
                   %(avEb,avEnb,avE,avELECT,avVDW)
            open('%s/solv/namd/%s.II.dat'%(sd,group),'w').write(buf2+buf)
            return False
        elif task=='rotate forces':
            pin=open('%s/solv/sulfur_cage.rmsRot.dat'%sd)
            buf='#    N  electr   VdW            f_electr(4-6)            f_VdW(7-9)\n'
            for l in open('%s/solv/namd/%s.dat'%(sd,group)):
                if l[0]=='#': continue  #skip comment lines
                #retrieve electrostatic and VdW forces
                items=l.split()
                N=int( items[0] )
                ele,vdw=[ float(x) for x in items[1:3] ]
                fele=numpy.array([float(x) for x in items[3:6] ])
                fvdw=numpy.array([float(x) for x in items[6:9]])
                #retrieve rotation to native. We need the transpose to
                #rotate the coordinates, since ptraj outputs the
                #matrix to rotate the vectors.
                R=[ float(x) for x in pin.readline().split()[1:] ]
                R=numpy.array(R).reshape(3,3).transpose()
                #rotate forces
                fele=numpy.dot(R,fele)
                fvdw=numpy.dot(R,fvdw)
                #save to file
                buf+=' %5d %8.2f %6.2f'%(1+N,ele,vdw)
                buf+=' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n'%tuple(list(fele)+list(fvdw))
            open('%s/solv/namd/%s.rotd.dat'%(sd,group),'w').write(buf)
        elif task=='rotate forces by CA':
            pin=open('%s/solv/CA.rmsRot.dat'%sd)
            buf='#    N  electr   VdW            f_electr(4-6)            f_VdW(7-9)\n'
            for l in open('%s/solv/namd/%s.dat'%(sd,group)):
                if l[0]=='#': continue  #skip comment lines
                #retrieve electrostatic and VdW forces
                items=l.split()
                N=int( items[0] )
                ele,vdw=[ float(x) for x in items[1:3] ]
                fele=numpy.array([float(x) for x in items[3:6] ])
                fvdw=numpy.array([float(x) for x in items[6:9]])
                #retrieve rotation to native. We need the transpose to
                #rotate the coordinates, since ptraj outputs the
                #matrix to rotate the vectors.
                R=[ float(x) for x in pin.readline().split()[1:] ]
                R=numpy.array(R).reshape(3,3).transpose()
                #rotate forces
                fele=numpy.dot(R,fele)
                fvdw=numpy.dot(R,fvdw)
                #save to file
                buf+=' %5d %8.2f %6.2f'%(1+N,ele,vdw)
                buf+=' %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n'%tuple(list(fele)+list(fvdw))
            open('%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group),'w').write(buf)
        elif task=='compute average forces':
            """compute the average electrostatic and Van der Waals
            force for the different groups interacting with the Fe
            ion"""
            if group:
                buf='      GROUP      <f_ele>(2-4)         <f_VdW>(5-7)\n'
                buf+=' %10s'%group
                infile='%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group)
                N=0; sm=[0.0]*6
                for l in open(infile).readlines():
                    if l[0]=='#': continue
                    N+=1
                    items=[ float(x) for x in l.split()[3:] ]
                    for i in range(6): sm[i]+=items[i]                    
                for i in range(6): buf+=' %6.2f'%(sm[i]/N,)
                print buf
                return False
            buf='by superposition of sulfur cage\n'
            buf+='      GROUP      <f_ele>(2-4)         <f_VdW>(5-7)\n'
            for group in 'FE_all FE_ions FE_protein FE_S FE_WAT'.split():
                buf+=' %12s'%group
                infile='%s/solv/namd/%s.rotd.dat'%(sd,group)
                N=0; sm=[0.0]*6
                for l in open(infile).readlines():
                    if l[0]=='#': continue
                    N+=1
                    items=[ float(x) for x in l.split()[3:] ]
                    for i in range(6): sm[i]+=items[i]                    
                for i in range(6): buf+=' %6.2f'%(sm[i]/N,)
                buf+='\n'
            buf+='by superposition of CA atoms\n'
            buf+='      GROUP      <f_ele>(2-4)         <f_VdW>(5-7)\n'
            for group in 'FE_all FE_ions FE_protein FE_S FE_WAT'.split():
                buf+=' %10s'%group
                infile='%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group)
                N=0; sm=[0.0]*6
                for l in open(infile).readlines():
                    if l[0]=='#': continue
                    N+=1
                    items=[ float(x) for x in l.split()[3:] ]
                    for i in range(6): sm[i]+=items[i]                    
                for i in range(6): buf+=' %6.2f'%(sm[i]/N,)
                buf+='\n'
            print buf
            return False

    def forces(self,conditions='T_340',task=None,group=None,zoom=15.0,
               inp=None,out=None, out2=None):
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        if task=='scatter plot elec byCA':
            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt
            infile='%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group)
            fxs=[]; fys=[]; fzs=[]; skip=9; iskip=9
            for l in open(infile).readlines():
                if l[0]=='#': continue
                if iskip!=skip: iskip+=1; continue
                iskip=0
                fx,fy,fz=[ float(x) for x in l.split()[3:6] ]
                fxs.append(fx); fys.append(fy); fzs.append(fz)
            fig=plt.figure()
            ax=Axes3D(fig)
            ax.scatter(fxs,fys,fzs)
            ax.set_xlabel('fx')
            ax.set_ylabel('fy')
            ax.set_zlabel('fz')
            plt.show()
            return False
        elif task=='scatter plot elec byCA for pymol':
            """create pdb file storing the electrostatic forces. The
            goal is to plot this as a point cloud around the Fe ion if
            we load the native state of rubredoxin"""
            from math import sqrt
            buf=''
            infile='%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group)
            fxs=[]; fys=[]; fzs=[]; fs=[]; F=0.0; avf=0.0
            for l in open(infile).readlines():
                if l[0]=='#': continue
                fx,fy,fz=[float(x) for x in l.split()[3:6]]
                fxs.append(fx); fys.append(fy); fzs.append(fz);
                f=sqrt(fx*fx + fy*fy + fz*fz);#  print f
                fs.append(f)
                if f>F: F=f
                avf+=f
            avf/=len(fxs)
            zoom=float(zoom)/F #re-scale the forces.
            native='%s/%s/unsolv/native.pdb'%(pjdout,self.id)
            for l in open(native).readlines():
                if 'FE   FE3' in l: break
            x,y,z=[float(r) for r in l[30:54].split()]
            #pdb.set_trace()
            #calculate B-factors like for the force modulus
            for i in range(len(fxs)):
                bf=fs[i]/F
                fxs[i]*=zoom; fys[i]*=zoom; fzs[i]*=zoom;
                fx=x+fxs[i]; fy=y+fys[i]; fz=z+fzs[i];
                xyz='%8.3f%8.3f%8.3f'%(fx,fy,fz)
                buf+='HETATM                        %s      %6.2f\n'%(xyz,bf)
            print buf
            return False
        elif task=='scatter plot of forces byCA for pymol':
            """create pdb file storing the forces. The goal is to plot
            this as a point cloud around the center of the input
            file"""
            from math import sqrt
            #read electrostatic force
            infile='%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group)
            fxs=[]; fys=[]; fzs=[]; fs=[]; F=0.0
            for l in open(infile).readlines():
                if l[0]=='#': continue
                fx,fy,fz=[float(x) for x in l.split()[3:6]]
                fxs.append(fx); fys.append(fy); fzs.append(fz);
                f=sqrt(fx*fx + fy*fy + fz*fz);  #print f
                fs.append(f)
                if f>F: F=f
            zoom2=float(zoom)/F #re-scale the forces.
            #calculate B-factors like for the force modulus
            buf=''
            for i in range(len(fxs)):
                bf=fs[i]/F
                fxs[i]*=zoom2; fys[i]*=zoom2; fzs[i]*=zoom2;
                fx=fxs[i]; fy=fys[i]; fz=fzs[i];
                xyz='%8.3f%8.3f%8.3f'%(fx,fy,fz)
                buf+='HETATM                        %s      %6.2f\n'%(xyz,bf)
            open(out,'w').write(buf)
            #read VdW force
            infile='%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group)
            fxs=[]; fys=[]; fzs=[]; fs=[]; F=0.0
            for l in open(infile).readlines():
                if l[0]=='#': continue
                fx,fy,fz=[float(x) for x in l.split()[6:9]]
                fxs.append(fx); fys.append(fy); fzs.append(fz);
                f=sqrt(fx*fx + fy*fy + fz*fz);  #print f
                fs.append(f)
                if f>F: F=f
            zoom2=float(zoom)/F #re-scale the forces.
            #calculate B-factors like for the force modulus
            buf=''
            for i in range(len(fxs)):
                bf=fs[i]/F
                fxs[i]*=zoom2; fys[i]*=zoom2; fzs[i]*=zoom2;
                fx=fxs[i]; fy=fys[i]; fz=fzs[i];
                xyz='%8.3f%8.3f%8.3f'%(fx,fy,fz)
                buf+='HETATM                        %s      %6.2f\n'%(xyz,bf)
            open(out2,'w').write(buf)
            return False
        elif task=='elec force modulus':
            from math import sqrt
            print '#modulus of electrostatic force'
            infile='%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group)
            for l in open(infile).readlines():
                if l[0]=='#': continue
                fx,fy,fz=[float(x) for x in l.split()[3:6]]
                #pdb.set_trace()
                f=sqrt(fx*fx + fy*fy + fz*fz)
                print '%6.2f'%f
            return False
        elif task=='force moduli':
            from math import sqrt
            infile='%s/solv/namd/%s.by_CA.rotd.dat'%(sd,group)
            buf='#f_ele  f_vdw\n'
            for l in open(infile).readlines():
                if l[0]=='#': continue
                fx,fy,fz=[float(x) for x in l.split()[3:6]]
                f=sqrt(fx*fx + fy*fy + fz*fz)
                buf+=' %6.2f'%f
                fx,fy,fz=[float(x) for x in l.split()[6:9]]
                f=sqrt(fx*fx + fy*fy + fz*fz)
                buf+=' %6.2f'%f
                buf+='\n'
            if out: open(out,'w').write(buf)
            else: print buf
            return False
        elif task=='residue-based average force projections':
            #load average force on whole multi-turn region
            for l in open( 'average_forces.dat' ).readlines():
                if 'mturn-'+inp in l:
                    items = [ float(x) for x in l.split()[1:] ]
                    Fele = numpy.array( items[0:3] )
                    Fvdw = numpy.array( items[3:] )
            #load average force on each residue and find projection
            buf='#GROUP    p_ELE  p_VDW\n'
            for l in open( 'xx-%s_avg_F.dat'%inp ).readlines():
                if l[0]=='#': continue                
                res = int( l.split('-')[0] )
                items = [ float(x) for x in l.split()[1:] ]
                fele = numpy.array( items[0:3] )
                p = numpy.dot(fele,Fele) / numpy.sum( (Fele*Fele) )
                buf+='%2d-%s %6.3f'%(res,inp,p)
                fvdw = numpy.array( items[3:] )
                p = numpy.dot(fvdw,Fvdw) / numpy.sum( (Fvdw*Fvdw) )
                buf+=' %6.3f\n'%p
            open(out,'w').write(buf)
            return False

    def phi_psi(self,conditions='T_340 solv phi_psi.00001_00843',
                task='average'):
        """go over phi and psi data files and find average. Care when
        angle jumps from 180 to -180"""
        from math import sqrt
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        buf='# resi <phi> std\n'
        for iphi in range(1,53):
            avg=0.00; std=0.00; prev_a=0.00; N=0
            for l in open('%s/phi%02d.dat'%(sd,iphi)).readlines():
                t,a=[float(x) for x in l.split()]
                if abs(prev_a-a)>180:
                    if a<0: a=a+360
                    else: a=a-360
                prev_a=a; avg+=a; std+=a*a; N+=1
            avg/=N; std=sqrt(std/N-avg*avg)
            buf+='%02d %6.1f %6.1f\n'%(iphi,avg,std)
        open('%s/avPhi.dat'%sd,'w').write(buf)
        buf='# resi <psi> std\n'
        for ipsi in range(2,54):
            avg=0.00; std=0.00; prev_a=0.00; N=0
            for l in open('%s/psi%02d.dat'%(sd,ipsi)).readlines():
                t,a=[float(x) for x in l.split()]
                if abs(prev_a-a)>180:
                    if a<0: a=a+360
                    else: a=a-360
                prev_a=a; avg+=a; std+=a*a; N+=1
            avg/=N; std=sqrt(std/N-avg*avg)
            buf+='%02d %6.1f %6.1f\n'%(ipsi,avg,std)
        open('%s/avPsi.dat'%sd,'w').write(buf)
        return False
    
    def rotamers(self,conditions='T_340 solv first_rotamer.00001_00843',
                resis=''):
        """go over ther rotamers average. Care when angle jumps from
        180 to -180"""
        from math import sqrt
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        buf='# resi <rot> std\n'
        for resi in resis.split():
            avg=0.00; std=0.00; prev_a=0.00; N=0
            for l in open('%s/%s.dat'%(sd,resi)).readlines():
                t,a=[float(x) for x in l.split()]
                if abs(prev_a-a)>180:
                    if a<0: a=a+360
                    else: a=a-360
                prev_a=a; avg+=a; std+=a*a; N+=1
            avg/=N; std=sqrt(std/N-avg*avg)
            buf+='%s %6.1f %6.1f\n'%(resi,avg,std)
        open('%s/avg.dat'%sd,'w').write(buf)
        return False

    def interactions(self,conditions='T_340 unsolv',task=None,
                    inpf=None,inpf2=None,outf=None,outf2=None):
        if task=='substract energyMaps':
            """read first ref15-shel map, then 1fhh map, then
            substract. Output numbers will be those of ref15-shel
            sequence. """
            from copy import copy
            comments=''
            #read ref15-shel
            ref15=['',]*100;
            for i in range(100): ref15[i]=['',]*100;
            for l in open(inpf).readlines():
                if l[0]=='#': comments+=l; continue
                xx = [ float(x) for x in l.split() ]
                i = int(xx[0]); j = int(xx[1])
                ref15[i][j]=copy( xx[2:] )
            #read 1fhh                
            fhh=['',]*100;
            for i in range(100): fhh[i]=['',]*100;
            for l in open(inpf2).readlines():
                if l[0]=='#': continue
                xx = [ float(x) for x in l.split() ]
                i = int(xx[0]); j = int(xx[1])
                fhh[i][j]=copy( xx[2:] )
            #substract and output to file
            buf=comments
            for i in range(100):
                for j in range(100):
                    if ref15[i][j] and fhh[i+1][j+1]:
                        buf+='%5d %5d'%(i,j)
                        a=ref15[i][j]; b=fhh[i+1][j+1]; dd=[i,j]
                        for k in range( len(a) ):
                            buf+=' %10.4f'%( a[k]-b[k], )
                        buf+='\n'
            open(outf,'w').write(buf)
            return False
        if task=='substract forceMaps':
            """read first ref15-shel map, then 1fhh map, then
            substract. Output numbers will be those of ref15-shel
            sequence. """
            from copy import copy
            comments=''
            #read ref15-shel
            ref15=['',]*100;
            for i in range(100): ref15[i]=['',]*100;
            for l in open(inpf).readlines():
                if l[0]=='#': comments+=l; continue
                xx = [ float(x) for x in l.split() ]
                i = int(xx[0]); j = int(xx[1])
                ref15[i][j]=copy( xx[2:] )
            #read 1fhh                
            fhh=['',]*100;
            for i in range(100): fhh[i]=['',]*100;
            for l in open(inpf2).readlines():
                if l[0]=='#': continue
                xx = [ float(x) for x in l.split() ]
                i = int(xx[0]); j = int(xx[1])
                fhh[i][j]=copy( xx[2:] )
            #substract vales and output to file
            buf=comments
            for i in range(100):
                for j in range(100):
                    #notice the shift [+1][+1] in fhh array
                    if ref15[i][j] and fhh[i+1][j+1]:
                        buf+='%5d %5d'%(i,j)
                        a=ref15[i][j]; b=fhh[i+1][j+1]; dd=[i,j]
                        for k in range( len(a) ):
                            buf+=' %10.4f'%( a[k]-b[k], )
                        buf+='\n'
            open(outf,'w').write(buf)
            #substract absolute vales and output to file
            buf=comments
            for i in range(100):
                for j in range(100):
                    #notice the shift [+1][+1] in fhh array
                    if ref15[i][j] and fhh[i+1][j+1]:
                        buf+='%5d %5d'%(i,j)
                        a=ref15[i][j]; b=fhh[i+1][j+1]; dd=[i,j]
                        for k in range( len(a) ):
                            buf+=' %10.4f'%( abs(a[k])-abs(b[k]), )
                        buf+='\n'
            open(outf2,'w').write(buf)
            return False
        if task=='compute SASA':
            """pass a list of PDB files to run dssp"""
            from secondary.dssp import rundssp,dssp
            from math import sqrt
            sasa_av=[0.0]*100; sasa_dev=[0.0]*54; N=0
            for file in open(inpf).readlines()[0:10]:
                rundssp(file.strip(),'junk.dat')
                obj=dssp('junk.dat')
                for i in range(len(obj.resi)):
                    ires=int(obj.resi[i]['RESIDUE'].split()[0])
                    acc= float( obj.resi[i]['ACC'] )
                    sasa_av[ires] += acc
                    sasa_dev[ires] += acc*acc
                N +=1
            buf='#RES <SASA> dev(SASA)\n'
            for i in range(1,54):
                x=sasa_av[i]/N
                y=sqrt( sasa_dev[i]/N - x*x)
                buf+='%3d %6.2f %6.2f\n'%(i,x,y)
            open(outf,'w').write(buf)
        return False

    def fluctuations(self,conditions='T_340 unsolv',task=None,
                 inp=None,inp2=None,inp3=None,inp4=None,
                 out=None):
        from math import sqrt
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        if task=='average RMSF for heavy atoms':
            MSF=0.0; N=0
            for l in open( '%s/%s'%(sd,inp) ).readlines():
                if l[0]=='#':continue
                x = float( l.split()[1] )
                MSF += x*x; N += 1
            MSF = sqrt(MSF/N)
            print '%s heavy atoms <RMSF> = %5.2f'%(conditions,MSF)
        elif task=='correlation coeff for CA fluctuations':
            """find the correlation coefficient for two
            atomicfluct_CA_1-52.dat files"""
            from scipy.stats import linregress
            from utilities.readingWritingFiles import read_column
            fluct1 = read_column(inp,2,isFloat=True)
            fluct2 = read_column(inp2,2,isFloat=True)
            (a_s,b_s,r,tt,stderr) = linregress(fluct1,fluct2)
            print r
        return False
    
    def contacts(self,conditions='T_340 unsolv',task=None,
                 inp=None,inp2=None,inp3=None,inp4=None,
                 out=None):
        pt=self.conditions2dict(conditions); sd=pt['simdir']
        if task=='ratio of contacts':
            ij = []; comments = ''; rf = []; frM = 1.00
            i1=open('residue_map.dat')
            i2=open(inp)
            for l1 in i1.readlines():
                l2=i2.readline()
                if l1[0]=='#':
                    comments += l1
                    continue
                i,j,f1 = l1.split(); f1 = float(f1)
                i,j,f2 = l2.split(); f2 = float(f2)
                ij.append( '%4s %4s'%(i,j) )
                if f1 < 1.00:
                    if f2 < 1.00: rf.append(1.00) # small frequencies
                    else:
                        if f2 > frM: frM = f2
                        rf.append(f2)           # non-native contacts
                else:
                    if f2/f1 > frM: frM = f2/f1
                    rf.append(f2/f1)
            print 'max ratio = ',str(frM)
            #trace non-native contacts and save to buffer
            buf  = '# Ratio of residue pair contact number\n'
            buf += comments                
            for i in range( len(rf) ):
                buf += '%s %6.3f\n'%(ij[i],rf[i])
            open(out,'w').write(buf)
            return False
        if task=='difference of contact ratios':
            """n_thermo(T)/n_thermo(native) - n_meso(T)/n_meso(native)"""
            def populate_list(inpx):
                lst = [ '' ]*55
                for i in range(55): lst[i] = [ '' ]*55
                for l in open(inpx).readlines():
                    if l[0]=='#': continue
                    i,j,f = l.split()
                    lst[ int(i) ][ int(j) ] = float( f )
                return lst

            lmN = populate_list(inp)  #native mesophilic
            lmT = populate_list(inp2) #mesophilic at some Temperature
            ltN = populate_list(inp3) #native thermophilic
            ltT = populate_list(inp4) #thermophilic at some Temperature

            xm = 0.00; xM = 0.00;
            buf = '#difference of ratio contact number\n'
            for i in range(1,53):
                for j in range(1,53):
                    if i == j: continue
                    x = 0.00
                    f0 = ltN[i][j] * lmN[i+1][j+1] #align meso and thermo
                    f  = ltT[i][j] * lmN[i+1][j+1]
                    f -= ltN[i][j] * lmT[i+1][j+1]
                    if f0 < 1.00: #definition of non-native
                        if f > 1.00: x = f
                    else:
                        x = f / f0
                    if x > xM: xM = x
                    if x < xm: xm = x
                    buf += '%3d %3d %6.3f\n' % (i,j,x)
            print 'maximum difference = ',str(xM)
            print 'minimum difference = ',str(xm)            
            open(out,'w').write(buf)
            return False
        if task=='ratio of contact ratios':
            """n_thermo(T)/n_thermo(native) / n_meso(T)/n_meso(native)"""
            def populate_list(inpx):
                lst = [ '' ]*55
                for i in range(55): lst[i] = [ '' ]*55
                for l in open(inpx).readlines():
                    if l[0]=='#': continue
                    i,j,f = l.split()
                    lst[ int(i) ][ int(j) ] = float( f )
                return lst

            lmN = populate_list(inp)  #native mesophilic
            lmT = populate_list(inp2) #mesophilic at some Temperature
            ltN = populate_list(inp3) #native thermophilic
            ltT = populate_list(inp4) #thermophilic at some Temperature

            xm = 0.00; xM = 0.00;
            buf = '#ratio of ratio contact number\n'
            for i in range(1,53):
                for j in range(1,53):
                    if i == j: continue
                    x = 0.00
                    f  = ltN[i][j] * lmT[i+1][j+1] #align meso and thermo
                    f0 = ltT[i][j] * lmN[i+1][j+1]
                    if f0 < 1.00:
                        if f < 1.00: x = 1.00
                    else:
                        x = f / f0
                    if x > xM: xM = x
                    if x < xm: xm = x
                    buf += '%3d %3d %6.3f\n' % (i,j,x)
            print 'maximum ratio = ',str(xM)
            print 'minimum ratio = ',str(xm)            
            open(out,'w').write(buf)
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
            
def forces(*kargs,**kwargs):
    cl.cycleOverList('forces',*kargs,**kwargs)
            
def phi_psi(*kargs,**kwargs):
    cl.cycleOverList('phi_psi',*kargs,**kwargs)

def rotamers(*kargs,**kwargs):
    cl.cycleOverList('rotamers',*kargs,**kwargs)

def interactions(*kargs,**kwargs):
    cl.cycleOverList('interactions',*kargs,**kwargs)

def fluctuations(*kargs,**kwargs):
    cl.cycleOverList('fluctuations',*kargs,**kwargs)

def contacts(*kargs,**kwargs):
    cl.cycleOverList('contacts',*kargs,**kwargs)

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

def mdanalysis(*kargs,**kwargs):
    """make use of MDAnalysis toolkit to analyze trajectories"""
    import MDAnalysis
    jobList = ('fluctuations',
               )
    job = kwargs['job']
    if job not in jobList: raise RuntimeError('Job '+job+' not in jobList')
    if job=='fluctuations':
        from copy import deepcopy
        from myMDAnalysis.analysis.align import rms_fit_trj
        from myMDAnalysis.analysis.fluctuations import avg_conf,msf_trj
        from myMDAnalysis.analysis.listCmd import listCmd
        refFile=kwargs['ref'] ; trajFile=kwargs['traj'] ; topFile=kwargs['top'] ; outFile=kwargs['outFile']
        reference = MDAnalysis.Universe(topFile,refFile,topology_format='TOP')
        mobile = MDAnalysis.Universe(topFile,trajFile,topology_format='TOP')
        commands=listCmd(mobile)
        commands.append(rms_fit_trj,mobile,reference)
        average = deepcopy( reference )
        commands.append(avg_conf,mobile,average)
        commands.execute()
        commands.clearCommandList()
        commands.append(msf_trj,mobile,filename=outFile,average=average)
        commands.execute()
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='service provider for RdCp_RdPf project')
    parser.add_argument('service',help='requested service, the name of a function defined in this module')
    parser.add_argument('--kargs',help='required arguments of service. Ex: "arg1,arg2,arg3"')
    parser.add_argument('--kwargs',help='optional arguments of service. Ex: "arg1=val1,arg2=val2"')
    args = parser.parse_args()
    service=args.service
    reqargs=[]
    if args.kargs: reqargs=args.kargs.split(',')
    optargs={}
    if args.kwargs: optargs=dict( map( lambda x: x.split('='), args.kwargs.split(',') ) )
    exitCode=0
    try:
        locals()[service](*reqargs,**optargs)
    except:
        sys.stderr.write('Error: service '+service+' not found or error in service\n')
        exitCode=1
    sys.exit(exitCode)

