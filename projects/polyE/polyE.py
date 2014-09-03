'''
Created on Jun 3, 2013

@author: jbq
'''
#!/usr/bin/python

import argparse, sys, os, numpy, random
from math import exp, erfc, sqrt
from pdb import set_trace as tr
from tempfile import mkstemp

projd='/projects/research/polyE'
coded='/home/jbq/code/projects/polyE'

def vmdexec(cmds):
  """execute a set of commands by calling vmd"""
  handle,filename=mkstemp(dir='/tmp')
  open(filename,'w').write(cmds)
  os.system('vmd -dispdev text -e %s'%filename) # run vmd in the terminal
  os.system('/bin/rm %s'%filename) # clean-up


def multi_chains(*kargs,**kwargs):
  job = kwargs['job']
  if job=='create chains':
    """generate PDB's for polyE chains"""
    avg=331.*1.5;  std=331./4;  nchain=int(kwargs['nchain'])
    for ichain in range(1,nchain+1):
      seq='e'*int(random.gauss(avg,std)) # extended structure with number of residues according to gaussian distribution
      #seq='E'*int(random.gauss(avg,std)) # helical structure with number of residues according to gaussian distribution
      os.system('crankite-peptide -r 10x1 %s -o %s_%02d.pdb &'%(seq,kwargs['prefix'],ichain))

  else:
    print 'Job not recognized. Nothing to do'

def scattering(*kargs,**kwargs):
  """calculation of structure factors"""
  from scatter.scatter import SIIScF,storeSIIScF,FourierSIIScf,fourierRealSymmetric
  from amber.amber10 import crd
  job=kwargs['job']
  if job=='isf':
    """intermediate scattering function of all the hydrogens"""
    qq=[1.9,]
    if 'qlist' in kwargs.keys(): qq=[float(q) for q in kwargs['qlist'].split()] #Q-values
    dt=float(kwargs['dt']) #time in between frames, in picoseconds
    tl=float(kwargs['tl']) #limit t (in psec). Compute SIIScF(q,t) from t=0 up to t=tl
    nt=int(tl/dt)          #number of time points in which to divide tl
    Nhydr=int(kwargs['Nhydr']) #number of hydrogen atoms
    nb=Nhydr   #consider only these number of hydrogens
    tt=[]
    for i in range(nt): tt.append(i*dt)
    crdf=kwargs['crdfile']
    crdobj=crd(crdf,Nhydr)
    crdobj.loadframes(skip=0)
    scat=SIIScF(qq,crdobj.frames[:,0:nb],dt,tl,nt)
    #store all intermediate factors
    if 'storagedir' in kwargs.keys():
      os.system('mkdir -p %s'%kwargs['storagedir'])
      for i in range( nb ):
        storeSIIScF( '%s/%04d.dat'%(kwargs['storagedir'],i+1), scat[i], qq, tt )
    avscat=scat.sum(axis=0)/len(scat)
    storeSIIScF( '%s'%kwargs['isffile'], avscat, qq, tt ) #average over all hydrogens
    if 'sffile' in kwargs.keys():
      sf=FourierSIIScf(avscat,transform=fourierRealSymmetric)
      ww=[]; dw=4.136/(nt*dt)
      for i in range(nt): ww.append(i*dw)
      storeSIIScF(kwargs['sffile'],sf,qq,ww)


def single_chain(*kargs,**kwargs):
  job = kwargs['job']
  if job=='translate AMBER PDB':
    """generate a PDB with atom type understood by CHARMM"""
    buf=open(kwargs['inpdb']).read()
    for line in open(projd+'/nchain_1/charmm/atom_type_conversion_table.dat'):
      if line[0]=='#': continue
      amber_type=line[1:8]; charmm_type=line[11:18]
      if charmm_type!='OT1 GLU': buf=buf.replace(amber_type,charmm_type) #replace amber atom with charmm atom
    open(kwargs['outpdb'],'w').write(buf)
    print "NOTICE: Manually edit output pdb and:\n1. replace last 'O   GLU' with 'OT1 GLU'\n2. remove the Na+ atoms"

  elif job=='neutralize system':
    """use cionize pluging of VMD to neutralize the system
    type cionize in the vmd prompt to inquire about argument options"""
    cmds='''mol load psf %s pdb %s
    package require cionize
    namespace import ::cionize::*
    cionize -mol 0 -prefix %s -ris 3.0 -rii 5.0 -border 4.0 -gridspacing 1.0 -ions "{SOD 496 1}"
    quit'''%(kwargs['psf'],kwargs['pdb'],kwargs['prefix'])
    vmdexec(cmds)

  elif job=='add charges to PDB':
    """add the partial atomic charges to a charmm PDB
    The partial charges below were retrieved from top_all27_prot_lipid.rtf"""
    charge={'CA  GLU':0.07, 'CB  GLU':-0.18, 'CD  GLU':0.62, 'CG  GLU':-0.28, 'C   GLU':0.51, 'HA  GLU':0.09,
            'HB1 GLU':0.09, 'HB2 GLU':0.09, 'HG1 GLU':0.09, 'HG2 GLU':0.09, 'HN  GLU':0.31, 'HT1 GLU':0.33,
            'HT2 GLU':0.33, 'HT3 GLU':0.33, 'N   GLU':-0.47, 'OE1 GLU':-0.76, 'OE2 GLU':-0.76, 'O   GLU':-0.51,
            'OT1 GLU':-0.67, 'OT2 GLU':-0.67,
             }
    buf=''
    for line in open(kwargs['inpdb']).readlines():
      atom_type=line[13:20]
      Q='%6.2f'%charge[atom_type]
      buf+=line[0:60]+Q+line[66:] # place the charge in the field for occupancy
    open(kwargs['outpdb'],'w').write(buf)

  elif job=='set CoG at origin':
    """set the geometry center at the origin"""
    from tempfile import mkstemp
    cmds='''source %s/vmd_procs.tcl
    mol new %s
    set sel [atomselect top "all"]
    set cog [geom_center $sel]
    $sel moveby -$cog
    $sel writepdb %s 
    quit
    '''%(coded,kwargs['inpdb'],kwargs['outpdb'])
    vmdexec(cmds)

  elif job=='create CA harmonic restraints':
    """Create harmonic restraints on the CA atoms that will shrink the polymer chain.
    Given a certain CA, we select a point of space as the center of attraction, and provide a force constant.
    The force constant should grow weaker as we move away from the center of the chain. A C-C bond has a k=600Kcal/(mol*A^2), and a C-C bond has
    length fluctuations of the order of 0.1A, thus, the typical force is ~600*0.1=60Kcal/(mol*A).
    The CA in the tip of the polymer are about 750A from the center. Furthermore, assume they are attracted to a point that lies in the surface
    of a sphere of radius 150A from the center of the polymer chain. Thus, let's assign a force constant such that 60 = k * (750+150), or k=0.067.
    This is the equation we will use to set the spring constant for all CA atoms, with 10 instead of 60.
    For a given CA lying a distance r from the center of the chain, we will select as center of attraction a point in the surface of a sphere of radius
    r/4 and with center in the center of the chain.
    """
    def surface_point(radius):
      """random point on the surface of a sphere of radius 'radius'"""
      z=random.uniform(-1,1)
      chi=random.uniform(0,2*numpy.pi)
      x=numpy.sqrt(1-z*z)*numpy.cos(chi)
      y=numpy.sqrt(1-z*z)*numpy.sin(chi)
      return radius*numpy.array([x,y,z])
    force=10. #typical force felt by each CA towards the center of attraction when the chain is stretched, in Kcal/(mol*A)
    ratio=4.
    if 'force' in kwargs.keys(): force=float(kwargs['force'])
    caxyz=[] # coordinates of CA atoms
    for line in open(kwargs['pdbref']).readlines():
      if 'CA  GLU' in line:
        xyz=[float(x) for x in line[30:54].split()]
        caxyz.append(xyz)
    caxyz=numpy.array(caxyz)
    middle=len(caxyz)/2 # middle index
    offset=0
    if len(caxyz)%2: offset=1 # odd number of CA atoms
    cog=(caxyz[middle]+caxyz[middle+1+offset])/2. # assumed center of the chain
    leftr=[]; rightr=[]; leftk=[]; rightk=[]
    for i in range(middle):
      index=i # index of CA 'to the left' of cog
      rel=caxyz[index]-cog # vector from CA to center of the chain
      r=numpy.linalg.norm(rel) # distance to center of chain
      leftk.append(force/(r+r/ratio)) # force constant such that typical forces are 10Kcal/(mol*A)
      leftr.append(cog+surface_point(r/ratio)) # center of attraction for the CA atom
      index=middle+i+offset # index of CA 'to the right' of cog
      rel=caxyz[index]-cog
      r=numpy.linalg.norm(rel)
      rightk.append(force/(r+r/ratio))
      rightr.append(cog+surface_point(r/ratio))
    k=numpy.array(leftk+rightk) # list of force constants
    caref=numpy.array(leftr+rightr) # list of center of attractions
    if offset:
      k=numpy.array(leftk+[0,]+rightk) # don't set harmonic restraint for middle CA atom
      caref=numpy.array(leftr+caxyz[middle]+rightr)
    bufk='' # buffer for PDB file containing force constant values
    bufr='' # buffer for PDB file containing constraint reference positions
    index=0 # CA index
    for line in open(kwargs['pdbref']).readlines():
      if 'CA  GLU' in line:
        bufk+=line[0:46] + '%8.3f'%k[index] + line[54:] # substitute Z-coordinate with force constant
        bufr+=line[0:30] + '%8.3f%8.3f%8.3f'%tuple(caref[index]) + line[54:] # substitute for center of attraction
        index+=1
      else:
        if line[0:5]=='ATOM ':
          bufk+=line[0:46] + '%8.3f'%(0.,) + line[54:] # make sure force constant is zero
        else:
          bufk+=line
        bufr+=line
    open(kwargs['conskfile'],'w').write(bufk)
    open(kwargs['consref'],'w').write(bufr)

  else:
    print 'Job not recognized. Nothing to do'

if __name__=='__main__':
  parser = argparse.ArgumentParser(description='service provider for dhfr_solv project')
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

