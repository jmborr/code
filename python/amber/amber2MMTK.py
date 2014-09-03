# Written by Jose Borreguero
# last revision:
#
import sys,os,re,numpy,Numeric
import pdb
from copy import copy
from Scientific.IO.NetCDF import *
from MMTK import PDB, Trajectory, Units
from MMTK.ParticleProperties import Configuration,ParticleVector
from MMTK.Trajectory import TrajectoryOutput,SnapshotGenerator
from utilities.small_utilities import Bye

basename=os.path.basename
#alternative names in Amber and MMTK
alt={'CYM':'CYX',}

def fixPDB(prmpdb,pdbfixed=None):
    """try to convert prmpdb to comply with MMTK requirements"""
    fixes=[]
    if not pdbfixed: pdbfixed=prmpdb+'.fixed'
    pat=re.compile('ATOM\s|HETATM\s')
    pdborig=open(prmpdb).read() ; pdb=pdborig
    for a,m in alt.items():
        if pdb.find(a)>=0:
            fixes.append(a+' switched to '+m)
            pdb=pdb.replace(a,m)
    natoms=len(pat.findall(pdborig))
    if natoms!=len(pat.findall(pdb)):
        buf='Could not convertfix %s to MMTK format'+basename(prmpdb)
        sys.stderr.write(buf+'\n')
        raise systemExit,1
    open(pdbfixed,'w').write(pdb)
    return pdbfixed,fixes,natoms

def getBox(prmcrd,natoms,format='crd'):
     """guess Box size from trajectory"""
     crdobj=crd(prmcrd,natoms)    #create trajectory object
     crdobj.readframe()
     box=crdobj.frame.box
     if box.any(): return box
     return None       
    
class frame:

    def __init__(self):
        self.xyz=None  #coordinates
        self.box=numpy.array([0.0,0.0,0.0])  #box size
        self.fn=0      #frame number
        self.natoms=None
        
    def readFrame(self,fp,natoms):
        """read certain number of atoms"""
        nread=0 ; xyz=[] ; box=[] ; nmax=natoms*3
        pos=fp.tell()
        if len(fp.readline())<2: return None #reached end of file
        fp.seek(pos)
        self.fn+=1
        while nread<nmax:
            line=fp.readline()
            if not line:
                sys.stderr.write('ERROR while reading frame %d\n'%self.fn)
                sys.exit(1)
            row=line.split() ; nread+=len(row)
            xyz+=[float(x) for x in row]
        if nread>nmax: box=xyz[nmax:]
        pos=fp.tell()
        row=fp.readline().split()
        if len(row)==3:
            box+=[float(x) for x in row]
        else:
            fp.seek(pos) #rewind
        self.xyz=numpy.array(xyz[:nmax]).reshape(natoms,3)
        self.xyz.reshape(natoms,3)
        if box: self.box=numpy.array(box)
        self.natoms=natoms
        return self

    def write(self):
        buf='generated with amber2MMTK.frame.write\n'
        flat=list(self.xyz.flatten()) ; flat.reverse() ; n=0
        while flat:
            buf+='%8.3f'%flat.pop() ; n+=1
            if not n%10: buf+='\n'
        if buf[-1]!='\n': buf+='\n'
        if self.box.any():
            flat=list(self.box)
            while flat: buf+='%8.3f'%flat.pop()
            buf+='\n'
        return buf
        
class crd:
    def __init__(self,crdfile,natoms,format='crd'):
        self.idx=[]         #index the frames
        self.fp=None
        self.natoms=natoms
        try:
            if crdfile[-3:]=='.gz':
                import gzip
                self.fp=gzip.open(crdfile)
            else:
                self.fp=open(crdfile)
        except:
            sys.stderr.write('Could not open '+crdfile+'\n')
            raise SystemExit,1
        self.frame=frame()
        
    def readframe(self):
        if self.fp.tell()==0: self.fp.readline() #first line is comment
        self.idx.append(self.fp.tell())
        if not self.frame.readFrame(self.fp,self.natoms):
            self.idx.pop() #remove last index
            #print 'self.idx.pop()'
            return None
        return self
        
    def __iter__(self):
        while self.readframe(): yield self.frame

               
class prmcrd2prmcdf:
    """create amber NetCDF file from a crd file

    ptraj is quite lame at prmcrd to prmcdf converting """

    def __init__(self,outcdf,natoms,prmcrd,
                 prmvel=None,box=None,dtime=0.0):
        
        f=NetCDFFile(outcdf, 'w',NC_64BIT_OFFSET)

        glob_attr={'title':'','application':'AMBER',
                   'program':'sander','programVersion':'9.0',
                   'Conventions':'AMBER','ConventionVersion':'1.0'}
        for key,val in glob_attr: setattr(f,key,val)

        dims={'frame':'UNLIMITED', 'spatial':3, 'atom':natoms,
              'label':5, 'cell_spatial':3, 'cell_angular':3}
        for key,val in dims: f.createDimension(key,val)

        f.createVariable('spatial', 'char',('spatial',))
        f.createVariable('time','float',('frame',))
        f.createVariable('coordinates','float',('frame','atom','spatial'))
        f.createVariable('velocities','float',('frame','atom','spatial'))
        f.createVariable('cell_spatial','char',('cell_spatial',))
        f.createVariable('cell_angular','char',('cell_angular','label'))
        f.createVariable('cell_lenghts','double',('frame','cell_spatial'))
        f.createVariable('cell_angles','double',('frame','cell_angles'))

        settatr(f.variables['time'],'units','picosecond')
        settatr(f.variables['coordinates'],'units','angstrom')
        settatr(f.variables['velocities'],'units','angstrom/picosecond')
        settatr(f.variables['cell_lenghts'],'units','angstrom')
        settatr(f.variables['cell_angles'],'units','degree')

        f.variables['spatial']='xyz'
        f.variables['cell_spatial']='abc'
        f.variables['cell_angular']='alpha','beta','gamma'
        f.variables['cell_angles'][0:3]=90.0,90.0,90.0 #orthorombic
        
        crdobj=crd(prmcrd,natoms,format='crd')
        velobj=None
        if prmvel: velobj=crd(prmve,natoms,format='crd')
        nfr=0 ; time=0.0
        for cframe in crdobj:
            f.variables['time'][nfr]=time
            f.variables['coordinates'][nfr][:]=cframe.xyz
            if box:
                for j in range[3]:
                    f.variables['cell_lenghts'][j]=box[j]
            elif cframe.box.any():
                f.variables['cell_lenghts'][:]=cframe.box[j]
            vframe=None
            if velobj: vframe=velobj.readframe()
            f.variables['coordinates'][nfr][:]=vframe.xyz
            nframe+=1 ; time+=dtime
        f.close()
            
class amberReader:

    def __init__(self,traj,universe,natoms,prmcrd,**options):
        self.traj=traj
        self.universe=universe
        self.natoms=natoms
        self.prmcrd=prmcrd
        self.prmvel=None
        self.__dict__.update(options)

    def __call__(self):
        actions=[TrajectoryOutput(self.traj, ['all'], 0, None, 1)]
        snapshot=SnapshotGenerator(self.universe,actions=actions)
        crdobj=crd(self.prmcrd,self.natoms)
        velobj=None
        if self.prmvel: velobj=crd(self.prmvel,self.natoms)
        time=0.0 ;
        nread=0  #for debugging purposes only
        array=Numeric.array; f=Numeric.Float; reshape=Numeric.reshape
        for cframe in crdobj:
            #print nread
            open('junk.crd','w').write(cframe.write())
            #print 'write(cframe.write())'
            if velobj:
                vframe=velobj.readframe().frame
                #print 'velobj.readframe().frame'
                #numpy->Numeric trick to avoid memory leaks
                da=reshape(array(list(vframe.xyz.flatten()),f)*Units.Ang,
                           vframe.xyz.shape)
                #print 'reshape(array(list(vframe.xyz.flatten()),f)'
                velocities=ParticleVector(self.universe,data_array=da)
                #print 'ParticleVector(self.universe,data_array=da)'
                self.universe.setVelocities(velocities)
                #print 'self.universe.setVelocities(velocities)'
                open('junk.vel','w').write(vframe.write())
                #print 'write(vframe.write())'
            da=reshape(array(list(cframe.xyz.flatten()),f)*Units.Ang,
                       cframe.xyz.shape)
            #print 'reshape(array(list(cframe.xyz.flatten()),f)'
            if cframe.box.any():
                self.box=array(list(cframe.box.flatten()),f)*Units.Ang
            #print 'array(list(cframe.box.flatten()),f)'
            config=Configuration(self.universe,data_array=da,cell=self.box)
            #print 'Configuration(self.universe,data_array=da,cell=box)'
            conf=self.universe.contiguousObjectConfiguration(None,config)
            #print 'universe.contiguousObjectConfiguration(None,config)'
            self.universe.setConfiguration(conf)
            #print 'self.universe.setConfiguration(conf)'
            snapshot(data={'time':time})
            #print'snapshot(data='
            time+=self.dtime
            #for debuggin purposes only
            nread+=1
        #print 'finished looping'



