'''
Created on Jun 3, 2013

@author: jbq
'''
#!/usr/bin/python

import argparse, sys, os, numpy, random, re
from math import exp, erfc, sqrt
from pdb import set_trace as tr
from tempfile import mkstemp

projd='/projects/research/POSS/me8t8'
coded='/home/jbq/code/projects/POSS/me8t8'

def vmdexec(cmds,redirect=False):
  """execute a set of commands by calling vmd"""
  handle,filename=mkstemp(dir='/tmp')
  open(filename,'w').write(cmds)
  cmd='vmd -dispdev text -e %s'%filename
  if redirect:
    return os.popen(cmd)
  else:
    os.system(cmd)
  return None

def mantidjobs(*kargs, **kwargs):
  """Jobs carried out in Mantid"""
  from mantid.simpleapi import (LoadNexus, SaveNexus, Scale, ScaleX)
  job=kwargs['job']

  if job=='create experimental Nexus files':
    rootd='/projects/research/POSS/me8t8/expdata'
    for T in (100,150,175,200,225,250,300,350):
      LoadDaveGrp(Filename='%s/exp%dK.dat'%(rootd,T),XAxisUnits='DeltaE',YAxisUnits='Empty',IsMicroEV=1,OutputWorkspace='exp%dK'%T)
      Rebin(InputWorkspace='exp%dK' %T,Params=[-0.14,0.0004,0.14],OutputWorkspace='exp%dK_2' %T)
      SaveNexus(InputWorkspace='exp%dK_2' %T,Filename='/projects/research/POSS/me8t8/expdata/exp%dK.nxs'%T)
    T=13
    LoadDaveGrp(Filename='%s/exp%dK.dat'%(rootd,T),XAxisUnits='DeltaE',YAxisUnits='Empty',IsMicroEV=1,OutputWorkspace='elastic')
    Rebin(InputWorkspace='elastic',Params=[-0.14,0.0004,0.14],OutputWorkspace='elastic')
    SaveNexus(InputWorkspace='elastic',Filename='/projects/research/POSS/me8t8/expdata/elastic.nxs')
    Scale(InputWorkspace='elastic',Factor=-1,Operation='Multiply',OutputWorkspace='resolution')
    SaveNexus(InputWorkspace='resolution',Filename='/projects/research/POSS/me8t8/expdata/resolution.nxs')

  else:
    print 'job not recognized, nothing to calculate'

def preparing(*kargs,**kwargs):
  """tasks related to preparation of simulation system"""

  job=kwargs['job']
  if job=='neutralize mol2 file':
    """the mol2 file has a residual charge, remove it"""
    inmol2=open(kwargs['inmol2'],'r');  line=inmol2.readline()
    while '@<TRIPOS>ATOM' not in line: line=inmol2.readline()
    line=inmol2.readline();  Q=0.;  natoms=0
    while '@<TRIPOS>BOND' not in line:
      Q+=float(line.split()[-1]);  natoms+=1;  line=inmol2.readline()
    Q=Q/natoms;  inmol2.seek(0);  line=inmol2.readline();  buf=line
    while '@<TRIPOS>ATOM' not in line:
      line=inmol2.readline();  buf+=line
    line=inmol2.readline()
    while '@<TRIPOS>BOND' not in line:
      q=float(line.split()[-1])-Q;  line=line[:67]+'%14.10f\n'%q;  buf+=line;  line=inmol2.readline()
    while line:
      buf+=line;  line=inmol2.readline()
    open(kwargs['outmol2'],'w').write(buf)

  elif job=='translate PDB from CSD to AMBER':
    """translate PDB from CSD to AMBER"""
    csd2amber={'Si1 ':' Si1', 'Si2 ':' Si2', 'Si2K':' Si3', 'Si2A':' Si4', 'Si2I':' Si5',
               'Si2B':' Si6', 'Si2J':' Si7', 'Si1J':' Si8', ' O1 ':' O2 ', ' O2K':' O3 ',
               ' O2A':' O4 ', ' O2I':' O10', ' O1A':' O1 ', ' O2B':' O9 ', ' O2J':' O7 ',
               ' O2 ':' O12', ' O1J':' O8 ', ' O1K':' O5 ', ' O1B':' O11', ' O1I':' O6 ',
               ' C2 ':' C2 ', ' C2K':' C7 ', ' C2A':' C6 ', ' C2I':' C3 ', ' C1 ':' C4 ',
               ' C2B':' C1 ', ' C2J':' C5 ', ' C1J':' C8 ', ' H3 ':' H6 ', ' H4 ':' H4 ',
               ' H2 ':' H5 ', ' H3K':' H19', ' H4K':' H21', ' H2K':' H20', ' H2A':' H18',
               ' H4A':' H16', ' H3A':' H17', ' H2I':' H8 ', ' H4I':' H7 ', ' H3I':' H9 ',
               ' H1A':' H10', ' H1 ':' H11', ' H1B':' H12', ' H2B':' H3 ', ' H3B':' H1 ',
               ' H4B':' H2 ', ' H3J':' H14', ' H4J':' H13', ' H2J':' H15', ' H1I':' H24',
               ' H1K':' H23', ' H1J':' H22',
               }
    amberorder=[' C1 ',' H1 ',' H2 ',' H3 ',' C2 ',' H4 ',' H5 ',' H6 ',' C3 ',' H7 ',' H8 ',
                ' H9 ',' C4 ',' H10',' H11',' H12',' C5 ',' H13',' H14',' H15',' C6 ',' H16',
                ' H17',' H18',' C7 ',' H19',' H20',' H21',' C8 ',' H22',' H23',' H24',' O1 ',
                ' Si1',' O2 ',' Si2',' O3 ',' Si3',' O4 ',' Si4',' O5 ',' O6 ',' Si5',' Si6',
                ' O7 ',' Si7',' O8 ',' Si8',' O9 ',' O10',' O11',' O12'] # order of the atoms in residue N11
    linespermolecule=52;  extralinespermolecule=8;  iatom_global=0;  iatom_molecule=0
    iresidue=1;  lines=[];  linespm={}; inpdb=open(kwargs['inpdb'],'r');  line=inpdb.readline()
    while line:
      if 'HETATM' in line:
        iatom_molecule+=1
        if iatom_molecule>linespermolecule:
          for i in range(extralinespermolecule-1): inpdb.readline()
          iatom_molecule=0
        else:
          iatom_global+=1;  line=line.replace('HETATM','ATOM  ');  line=line.replace('UNK','N11')
          amberatom=csd2amber[line[12:16]]
          line=line[:12]+amberatom+line[16:];  line=line[:6]+'%5d'%iatom_global+line[11:]
          line=line[:23]+'%3d'%iresidue+line[26:]; 
          linespm[amberatom]=line
          if not iatom_global%linespermolecule:
            for amberatom in amberorder: lines+=linespm[amberatom]
            linespm={}; iresidue+=1;  lines+='TER\n'
      line=inpdb.readline()
    lines+='END\n'
    open(kwargs['outpdb'],'w').write(''.join(lines))

  else:
    print 'job not recognized, nothing to calculate'


def single_molecule(*kargs,**kwargs):
  from utilities.readingWritingFiles import write_from_numpy
  job = kwargs['job']
  if job=='remove residual total charge':
    """Remove residual total charge"""
    pt=open(kwargs['prmtop'])
    buf=''
    line=pt.readline()
    while '%FLAG CHARGE' not in line:
      buf+=line
      line=pt.readline()
    buf+=line+pt.readline()
    line=pt.readline()
    charges=[]
    while '%FLAG MASS' not in line:
      charges+=[float(q) for q in line.split()]
      line=pt.readline()
    charges=numpy.array(charges)
    charges=charges-numpy.sum(charges)/len(charges)
    buf+=write_from_numpy(None,charges,format='%16.8e',ncols=5)
    buf+=line
    buf+=''.join(pt.readlines())
    open(kwargs['outprmtop'],'w').write(buf)

  elif job=='convert ptraj DCD to vmd DCD':
    """sassena cannot read a DCD trajectory output by ptraj. We must load it into VMD and then
    save it as DCD. Now sassena can read the new DCD trajectory"""
    cmds="""mol new %s waitfor all
mol addfile %s waitfor all
animate write dcd %s
quit"""%(kwargs['prmtop'],kwargs['indcd'],kwargs['outdcd'])
    vmdexec(cmds)
  else:
    print 'job not recognized, nothing to calculate'

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

