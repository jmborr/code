#!/usr/bin/python

import os,re,sys,pickle,gzip,pdb
import numpy
from tempfile import mkstemp,mkdtemp
from math import sqrt,pi,cos
from copy import copy,deepcopy

joink=os.path.join
basename=os.path.basename
dirname=os.path.dirname
isthere=os.path.exists
codedir=os.environ['CODED']
tmpdir=os.environ['TMPD']

try:
    AMBERHOME=os.environ['AMBERHOME10']
except:
    sys.stderr.write("WARNING: environment variable AMBERHOME10 not set\n")

ptrajexe='ptraj10' #link to appropriate executable

def extricateBC(fac,box):
    """
    remove periodic boundary conditions from a trajectory

    fac.shape = (#frames, #atoms, 3 ) is the trajectory

    Overwrite fac by taking the first snapshot as reference, then
    removing the periodic boundary conditions for the following
    snapshots
    
    """

    prev = fac[0]
    for ifr in range(1, len( fac ) ):
        dfr = (fac[ ifr ] - prev).transpose() #dfr.shape=(3,#atoms)
        prev = deepcopy(fac[ ifr ]) #store fac[ifr] for next round
        #extricate periodic boundary conditions
        for ix in range( 3 ):
            d = box[ix]; dh = d/2.0
            Y = dfr[ix]  # Y is just a reference to dfr[ix]
            Y[numpy.where(Y > dh )] -= d
            Y[numpy.where(Y < -dh)] += d
        fac[ ifr ] = fac[ ifr-1 ] + dfr.transpose()
        

def genRST(inpdb,leeway=0.2,rk2=10.0,rk3=10.0,task='by CB map',co=7.5):
    """generate NMR restraints

R1, R2, R3, R4 define a flat-welled parabola which becomes linear
beyond a specified distance. I.e.

Example
&rst
    iat=5,63,
    iresid=0,
        r1=0.0, r2=5.2, r3=5.6, r4=999.0,
         rk2=25.0, rk3=25.0,
 &end

         \                       /
          \                     /
           \                   /
            .                 .
              .             .
                 ._______.

           r1    r2      r3   r4   
           
      '\' = lower bound linear response region 
      '/' = lower bound linear response region 
      '.' = parobola
      '_' = flat region
      rk2 (rk3) is force constant of left (right) parabola.2
      
If you have determined lower and upper bounds from an NMR experiment,
those would typically correspond to R2 and R3. Note that the flat well
means that any value R2 <= value <= R3 is equally acceptable. R1 and
R4 define linear response regions. These are sometimes used so that
restraints that severely violate the lower and upper bounds don't tear
the structure up. A typical value of R1 is R2-2.0. And a typical value
of R4 is R3+2.0 (angstroms)."""
    if   isinstance(inpdb,str) : pin=open(inpdb)
    elif isinstance(inpdb,file): pin=inpdb

    template="""&rst
    iat=_ATOM1_,_ATOM2_,
    iresid=0,
        r1=0.0, r2=_R2_, r3=_R3_, r4=999.0,
         rk2=_RK2_, rk3=_RK3_,
 &end\n"""

    if task=='by CB map':
        from mypdb.smallThings import CBmap
        cbm=CBmap(pin,co=co)
        neigh=cbm['neighbors']; d=cbm['distances']; L=len(neigh)
        #find atomic number from pdblines
        atn=[]
        for l in cbm['pdblines']: atn.append(l[6:11])
        #create restraints
        buf=''
        for iat in range(L):
            nn=neigh[iat]; k=0
            for jat in nn:
                if jat<iat: continue #avoid duplicate the restraint
                restr=template
                center=d[iat][k]; wiggle=center*leeway/2
                repl={'_ATOM1_':atn[iat], '_ATOM2_':atn[jat],
                      '_R2_':center-wiggle, '_R3_':center+wiggle,
                      '_RK2_':rk2,'_RK3_':rk3}
                for key,val in repl.items():
                    restr=restr.replace(key,str(val))
                buf+=restr; k+=1
        return buf
    
class leap:

    """dealing with LEap"""
    
    def __init__(self):

        self.AAA=('ALA', 'ARG', 'ASN', 'ASP', 'ASH',\
                  'CYS', 'CYM', 'CYX', 'GLU', 'GLH',\
                  'GLN', 'GLY', 'HID', 'HIE', 'HIP',\
                  'ILE', 'LEU', 'LYS', 'MET', 'PHE',\
                  'PRO', 'SER', 'THR', 'TRP', 'TYR',\
                  'VAL')

        self.basenameAAA={
            'ALA':'ALA', 'ARG':'ARG', 'ASN':'ASN', 'ASP':'ASP',\
            'ASH':'ASP', 'CYS':'CYS', 'CYM':'CYS', 'CYX':'CYS',\
            'GLU':'GLU', 'GLH':'GLU', 'GLN':'GLN', 'GLY':'GLY',\
            'HID':'HIS', 'HIE':'HIS', 'HIP':'HIS', 'ILE':'ILE',\
            'LEU':'LEU', 'LYS':'LYS', 'MET':'MET', 'PHE':'PHE',\
            'PRO':'PRO', 'SER':'SER', 'THR':'THR', 'TRP':'TRP',\
            'TYR':'TYR', 'VAL':'VAL'}
        self.isstandarAAAname={
            'ALA':True, 'ARG':True, 'ASN':True, 'ASP':True, \
            'ASH':False,'CYS':True, 'CYM':False,'CYX':False,\
            'GLU':True ,'GLH':False,'GLN':True, 'GLY':True, \
            'HID':False,'HIE':False,'HIP':False,'ILE':True, \
            'LEU':True, 'LYS':True, 'MET':True, 'PHE':True, \
            'PRO':True, 'SER':True, 'THR':True, 'TRP':True, \
            'TYR':True, 'VAL':True}

class tleap(leap):

    """derived class from LEaP"""

    def __init__(self,cmd=[],seq=''):
        self.cmd=cmd #list of commands
        if seq: self.cmd+=seq.split('\n')
        
    def go(self,script=None,fake=None):
        """execute the stored commands"""
        if not script: xxx,script=mkstemp(dir=tmpdir)
        open(script,'w').write('\n'.join(self.cmd)+'\nquit')
        cmd='tleap -f  '+script
        if fake:
            print cmd
        else:
            os.system(cmd)
        os.system('/bin/rm '+script)

class top:
    """manage topology file"""

    def __init__(self,topf):
        self.top=topf
        self.FLAG={}
        self.parse()

    def parse(self):

        """http://ambermd.org/formats.html"""

        """read in 'FLAG POINTERS'
FORMAT(12i6)  NATOM,  NTYPES, NBONH,  MBONA,  NTHETH, MTHETA,
              NPHIH,  MPHIA,  NHPARM, NPARM,  NEXT,   NRES,
              NBONA,  NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA,
              NATYP,  NPHB,   IFPERT, NBPER,  NGPER,  NDPER,
              MBPER,  MGPER,  MDPER,  IFBOX,  NMXRS,  IFCAP
  NATOM  : total number of atoms 
  NTYPES : total number of distinct atom types
  NBONH  : number of bonds containing hydrogen
  MBONA  : number of bonds not containing hydrogen
  NTHETH : number of angles containing hydrogen
  MTHETA : number of angles not containing hydrogen
  NPHIH  : number of dihedrals containing hydrogen
  MPHIA  : number of dihedrals not containing hydrogen
  NHPARM : currently not used
  NPARM  : currently not used
  NEXT   : number of excluded atoms
  NRES   : number of residues
  NBONA  : MBONA + number of constraint bonds
  NTHETA : MTHETA + number of constraint angles
  NPHIA  : MPHIA + number of constraint dihedrals
  NUMBND : number of unique bond types
  NUMANG : number of unique angle types
  NPTRA  : number of unique dihedral types
  NATYP  : number of atom types in parameter file, see SOLTY below
  NPHB   : number of distinct 10-12 hydrogen bond pair types
  IFPERT : set to 1 if perturbation info is to be read in
  NBPER  : number of bonds to be perturbed
  NGPER  : number of angles to be perturbed
  NDPER  : number of dihedrals to be perturbed
  MBPER  : number of bonds with atoms completely in perturbed group
  MGPER  : number of angles with atoms completely in perturbed group
  MDPER  : number of dihedrals with atoms completely in perturbed groups
  IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
  NMXRS  : number of atoms in the largest residue
  IFCAP  : set to 1 if the CAP option from edit was specified
"""
        regexp='FORMAT\((\d+)\w(\d+)'
        all=open(self.top)
        
        keys='NATOM  NTYPES NBONH  MBONA  NTHETH MTHETA '+\
             'NPHIH  MPHIA  NHPARM NPARM  NEXT   NRES '  +\
             'NBONA  NTHETA NPHIA  NUMBND NUMANG NPTRA ' +\
             'NATYP  NPHB   IFPERT NBPER  NGPER  NDPER ' +\
             'MBPER  MGPER  MDPER  IFBOX  NMXRS  IFCAP '
        keys=keys.split()
        all.seek(0); line=''
        while('FLAG POINTERS') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values)<len(keys):
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.ITITL={}
        for i in range(len(keys)): self.ITITL[ keys[i] ] = values[i]
        self.FLAG['ITITL']=self.ITITL

        """read in FLAG ATOM_NAME"""
        all.seek(0); line=''
        while('FLAG ATOM_NAME') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[line[i:i+D] for i in range(0,L,D)]
        self.IGRAPH=values[:self.ITITL['NATOM']]        
        self.FLAG['IGRAPH']=self.IGRAPH

        """read in FLAG CHARGE"""
        all.seek(0); line=''
        while('FLAG CHARGE') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.CHRG=values[:self.ITITL['NATOM']]        
        self.FLAG['CHRG']=self.CHRG
        
        """read in 'FLAG MASS'"""
        all.seek(0); line=''
        while('FLAG MASS') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.AMASS=values[:self.ITITL['NATOM']]
        self.FLAG['AMASS']=self.AMASS


        """read in FLAG ATOM_TYPE_INDEX
          IAC : index for the atom types involved in Lennard Jones
           (6-12) interactions.  See ICO below.  """
        all.seek(0); line=''
        while('FLAG ATOM_TYPE_INDEX') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.IAC=values[:self.ITITL['NATOM']]
        self.FLAG['IAC']=self.IAC


        """read in FLAG NUMBER_EXCLUDED_ATOMS
          NUMEX : total number of excluded atoms for atom "i".  See
           NATEX below."""
        all.seek(0); line=''
        while('FLAG NUMBER_EXCLUDED_ATOMS') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[line[i:i+D] for i in range(0,L,D)]
        self.NUMEX=values[:self.ITITL['NATOM']]        
        self.FLAG['NUMEX']=self.NUMEX


        """read in FLAG NONBONDED_PARM_INDEX
        ICO(i), i=1,NTYPES*NTYPES        
        ICO : provides the index to the nonbon parameter arrays CN1,
           CN2 and ASOL, BSOL.  All possible 6-12 or 10-12 atoms type
           interactions are represented.  NOTE: A particular atom type
           can have either a 10-12 or a 6-12 interaction, but not
           both.  The index is calculated as follows: index =
           ICO(NTYPES*(ICO(i)-1)+ICO(j)) If index is positive, this is
           an index into the 6-12 parameter arrays (CN1 and CN2)
           otherwise it is an index into the 10-12 parameter arrays
           (ASOL and BSOL)."""
        nread=self.ITITL['NTYPES']*self.ITITL['NTYPES']
        all.seek(0); line=''
        while('FLAG NONBONDED_PARM_INDEX') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.ICO=values[:nread]
        self.FLAG['ICO']=self.ICO

        
        """read in FLAG RESIDUE_LABEL
        LABRES : the residue labels"""
        all.seek(0); line=''
        while('FLAG RESIDUE_LABEL') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NRES']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[line[i:i+D] for i in range(0,L,D)]
        self.LABRES=values[:self.ITITL['NRES']]        
        self.FLAG['LABRES']=self.LABRES

        
        """read in FLAG RESIDUE_POINTER"""
        all.seek(0); line=''
        while('FLAG RESIDUE_POINTER') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NRES']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.IPRES=values[:self.ITITL['NRES']]        
        self.FLAG['IPRES']=self.IPRES

        
        """read in FLAG BOND_FORCE_CONSTANT
        RK: force constant for the bonds of each type, kcal/mol"""
        all.seek(0); line=''
        while('FLAG BOND_FORCE_CONSTANT') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NUMBND']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.RK=values[:self.ITITL['NUMBND']]
        self.FLAG['NUMBNK']=self.RK

        
        """read in FLAG BOND_EQUIL_VALUE
        REQ: the equilibrium bond length for the bonds of each type,
        angstroms """
        all.seek(0); line=''
        while('FLAG BOND_EQUIL_VALUE') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NUMBND']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.REQ=values[:self.ITITL['NUMBND']]
        self.FLAG['REQ']=self.REQ

        
        """read in FLAG ANGLE_FORCE_CONSTANT
        TK: force constant for the angles of each type, kcal/mol A**2"""
        all.seek(0); line=''
        while('FLAG ANGLE_FORCE_CONSTANT') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NUMANG']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.TK=values[:self.ITITL['NUMANG']]
        self.FLAG['TK']=self.TK

        
        """read in FLAG ANGLE_EQUIL_VALUE
        TEQ: the equilibrium angle for the angles of each type, radians"""
        all.seek(0); line=''
        while('FLAG ANGLE_EQUIL_VALUE') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NUMANG']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.TEQ=values[:self.ITITL['NUMANG']]
        self.FLAG['TEQ']=self.TEQ

        
        """read in FLAG DIHEDRAL_FORCE_CONSTANT
        PK: force constant for the dihedrals of each type, kcal/mol"""
        all.seek(0); line=''
        while('FLAG DIHEDRAL_FORCE_CONSTANT') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NPTRA']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.PK=values[:self.ITITL['NPTRA']]
        self.FLAG['PK']=self.PK

        
        """read in FLAG DIHEDRAL_PERIODICITY
        PN: periodicity of the dihedral of a given type"""
        all.seek(0); line=''
        while('FLAG DIHEDRAL_PERIODICITY') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NPTRA']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.PN=values[:self.ITITL['NPTRA']]
        self.FLAG['PN']=self.PN

        
        """read in FLAG DIHEDRAL_PHASE
        PHASE: phase of the dihedral of a given type, radians"""
        all.seek(0); line=''
        while('FLAG DIHEDRAL_PHASE') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NPTRA']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.PHASE=values[:self.ITITL['NPTRA']]
        self.FLAG['PHASE']=self.PHASE

        
        """read in FLAG SOLTY
        SOLTY  : currently unused (reserved for future use)"""
        all.seek(0); line=''
        while('FLAG SOLTY') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATYP']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.SOLTY=values[:self.ITITL['NATYP']]
        self.FLAG['SOLTY']=self.SOLTY

        
        """read in FLAG LENNARD_JONES_ACOEF

          LENNARD_JONES_ACOEF: Lennard Jones r**12 terms for all
           possible atom type interactions, indexed by ICO and IAC;
           for atom i and j where i < j, the index into this array is
           as follows (assuming the value of ICO(index) is positive):
           CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j)))."""
        nread=self.ITITL['NTYPES']*(self.ITITL['NTYPES']+1)/2
        all.seek(0); line=''
        while('FLAG LENNARD_JONES_ACOEF') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.CN1=values[:nread]
        self.FLAG['CN1']=self.CN1

        
        """read in FLAG LENNARD_JONES_BCOEF

          LENNARD_JONES_BCOEF: Lennard Jones r**6 terms for all
           possible atom type interactions.  Indexed like CN1
           above."""
        nread=self.ITITL['NTYPES']*(self.ITITL['NTYPES']+1)/2
        all.seek(0); line=''
        while('FLAG LENNARD_JONES_BCOEF') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.CN2=values[:nread]
        self.FLAG['CN2']=self.CN2

        
        """NOTE: the atom numbers in the following arrays that
        describe bonds, angles, and dihedrals are coordinate array
        indexes for runtime speed. The true atom number equals the
        absolute value of the number divided by three, plus one. In
        the case of the dihedrals, if the fourth atom is negative,
        this implies that the dihedral is an improper. If the third
        atom is negative, this implies that the end group interations
        are to be ignored. End group interactions are ignored, for
        example, in dihedrals of various ring systems (to prevent
        double counting of 1-4 interactions) and in multiterm
        dihedrals."""

        """read in FLAG BONDS_INC_HYDROGEN
        (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
        IBH    : atom involved in bond "i", bond contains hydrogen
        JBH    : atom involved in bond "i", bond contains hydrogen
        ICBH   : index into parameter arrays RK and REQ"""
        mread=self.ITITL['NBONH']; nread=3*mread
        all.seek(0); line=''
        while('FLAG BONDS_INC_HYDROGEN') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.IBH =values[0::3];  self.FLAG['IBH'] =self.IBH
        self.JBH =values[1::3];  self.FLAG['JBH'] =self.JBH
        self.ICBH=values[2::3];  self.FLAG['ICBH']=self.ICBH

        
        """read in FLAG BONDS_WITHOUT_HYDROGEN
        (IB(i),JB(i),ICB(i), i=1,NBONA)
        IB : atom involved in bond "i", bond does not contain hydrogen
        JB : atom involved in bond "i", bond does not contain hydrogen
        ICB: index into parameter arrays RK and REQ"""
        mread=self.ITITL['NBONA']; nread=3*mread
        all.seek(0); line=''
        while('FLAG BONDS_WITHOUT_HYDROGEN') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.IB =values[0::3];  self.FLAG['IB'] =self.IB
        self.JB =values[1::3];  self.FLAG['JB'] =self.JB
        self.ICB=values[2::3];  self.FLAG['ICB']=self.ICB

        
        """read in FLAG ANGLES_INC_HYDROGEN
        (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)
        ITH : atom involved in angle "i", angle contains hydrogen
        JTH : atom involved in angle "i", angle contains hydrogen
        KTH : atom involved in angle "i", angle contains hydrogen
        ICTH: index into parameter arrays TK and TEQ for angle
           ITH(i)-JTH(i)-KTH(i)"""
        mread=self.ITITL['NTHETH']; nread=4*mread
        all.seek(0); line=''
        while('FLAG ANGLES_INC_HYDROGEN') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.ITH =values[0::4];  self.FLAG['ITH'] =self.ITH
        self.JTH =values[1::4];  self.FLAG['JTH'] =self.JTH
        self.KTH =values[2::4];  self.FLAG['KTH'] =self.KTH
        self.ICTH=values[3::4];  self.FLAG['ICTH']=self.ICTH

        
        """read in FLAG ANGLES_WITHOUT_HYDROGEN
        (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)
        IT : atom involved in angle "i", angle does not contain hydrogen
        JT : atom involved in angle "i", angle does not contain hydrogen
        KT : atom involved in angle "i", angle does not contain hydrogen
        ICT: index into parameter arrays TK and TEQ for angle
           IT(i)-JT(i)-KT(i)"""
        mread=self.ITITL['NTHETA']; nread=4*mread
        all.seek(0); line=''
        while('FLAG ANGLES_WITHOUT_HYDROGEN') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.IT =values[0::4];  self.FLAG['IT'] =self.IT
        self.JT =values[1::4];  self.FLAG['JT'] =self.JT
        self.KT =values[2::4];  self.FLAG['KT'] =self.KT
        self.ICT=values[3::4];  self.FLAG['ICT']=self.ICT

        
        """read in FLAG DIHEDRALS_INC_HYDROGEN
        (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i), i=1,NPHIH)
        IPH : atom involved in dihedral "i", dihedral contains hydrogen
        JPH : atom involved in dihedral "i", dihedral contains hydrogen
        KPH : atom involved in dihedral "i", dihedral contains hydrogen
        LPH : atom involved in dihedral "i", dihedral contains hydrogen
        ICPH: index into parameter arrays PK, PN, and PHASE for
           dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i)"""
        mread=self.ITITL['NPHIH']; nread=5*mread
        all.seek(0); line=''
        while('FLAG DIHEDRALS_INC_HYDROGEN') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.IPH =values[0::5];  self.FLAG['IPH'] =self.IPH
        self.JPH =values[1::5];  self.FLAG['JPH'] =self.JPH
        self.KPH =values[2::5];  self.FLAG['KPH'] =self.KPH
        self.LPH =values[3::5];  self.FLAG['LPH'] =self.LPH
        self.ICPH=values[4::5];  self.FLAG['ICPH']=self.ICPH

        
        """read in FLAG DIHEDRALS_WITHOUT_HYDROGEN
        (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)
        IP : atom involved in dihedral "i", dihedral does not contain hydr
        JP : atom involved in dihedral "i", dihedral does not contain hydr
        KP : atom involved in dihedral "i", dihedral does not contain hydr
        LP : atom involved in dihedral "i", dihedral does not contain hydr
        ICP: index into parameter arrays PK, PN, and PHASE for
           dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i).  Note, if the
           periodicity is negative, this implies the following entry
           in the PK, PN, and PHASE arrays is another term in a
           multitermed dihedral."""
        mread=self.ITITL['NPHIA']; nread=5*mread
        all.seek(0); line=''
        while('FLAG DIHEDRALS_WITHOUT_HYDROGEN') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < nread:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.IP =values[0::5];  self.FLAG['IP'] =self.IP
        self.JP =values[1::5];  self.FLAG['JP'] =self.JP
        self.KP =values[2::5];  self.FLAG['KP'] =self.KP
        self.LP =values[3::5];  self.FLAG['LP'] =self.LP
        self.ICP=values[4::5];  self.FLAG['ICP']=self.ICP

        
        """read in FLAG EXCLUDED_ATOMS_LIST
        NATEX : the excluded atom list.  To get the excluded list for
           atom "i" you need to traverse the NUMEX list, adding up all
           the previous NUMEX values, since NUMEX(i) holds the number
           of excluded atoms for atom "i", not the index into the
           NATEX list.  Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
           excluded atoms are NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i))."""
        all.seek(0); line=''
        while('FLAG EXCLUDED_ATOMS_LIST') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NEXT']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.NATEX=values[:self.ITITL['NEXT']]
        self.FLAG['NATEX']=self.NATEX

        
        """read in FLAG HBOND_ACOEF
        ASOL: the value for the r**12 term for hydrogen bonds of all
           possible types.  Index into these arrays is equivalent to
           the CN1 and CN2 arrays, however the index is negative.  For
           example, for atoms i and j, with i < j, the index is
           -ICO(NTYPES*(IAC(i)-1+IAC(j))."""
        all.seek(0); line=''
        while('FLAG HBOND_ACOEF') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NPHB']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.ASOL=values[:self.ITITL['NPHB']]
        self.FLAG['ASOL']=self.ASOL

        
        """read in FLAG HBOND_BCOEF
        BSOL: the value for the r**10 term for hydrogen bonds of all
        possible types.  Indexed like ASOL."""
        all.seek(0); line=''
        while('FLAG HBOND_BCOEF') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NPHB']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.BSOL=values[:self.ITITL['NPHB']]
        self.FLAG['BSOL']=self.BSOL

        
        """read in FLAG HBCUT
        HBCUT: no longer in use"""
        all.seek(0); line=''
        while('FLAG HBCUT') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NPHB']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[float(line[i:i+D]) for i in range(0,L,D)]
        self.HBCUT=values[:self.ITITL['NPHB']]
        self.FLAG['HBCUT']=self.HBCUT

        
        """read in FLAG AMBER_ATOM_TYPE
        ISYMBL : the AMBER atom types for each atom"""
        all.seek(0); line=''
        while('FLAG AMBER_ATOM_TYPE') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[line[i:i+D] for i in range(0,L,D)]
        self.ISYMBL=values[:self.ITITL['NATOM']]
        self.FLAG['ISYMBL']=self.ISYMBL

        
        """read in FLAG TREE_CHAIN_CLASSIFICATION
        ITREE: the list of tree joining information, classified into
           five types.  M -- main chain, S -- side chain, B -- branch
           point, 3 -- branch into three chains, E -- end of the chain"""
        all.seek(0); line=''
        while('FLAG TREE_CHAIN_CLASSIFICATION') not in line:
            line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[line[i:i+D] for i in range(0,L,D)]
        self.ITREE=values[:self.ITITL['NATOM']]
        self.FLAG['ITREE']=self.ITREE

        
        """read in FLAG JOIN_ARRAY
        JOIN : tree joining information, potentially used in ancient
           analysis programs.  Currently unused in sander or gibbs."""
        all.seek(0); line=''
        while('FLAG JOIN_ARRAY') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.JOIN=values[:self.ITITL['NATOM']]
        self.FLAG['JOIN']=self.JOIN

        
        """read in FLAG IROTAT
        IROTAT : apparently the last atom that would move if atom i
           was rotated, however the meaning has been lost over time.
           Currently unused in sander or gibbs."""
        all.seek(0); line=''
        while('FLAG IROTAT') not in line: line=all.readline()
        values=[] ; line=all.readline() #line specifying format
        N,D=[int(x) for x in re.compile(regexp).search(line).groups()]
        while len(values) < self.ITITL['NATOM']:
            line=all.readline().strip('\n'); L=D*(len(line)/D)
            values+=[int(line[i:i+D]) for i in range(0,L,D)]
        self.IROTAT=values[:self.ITITL['NATOM']]
        self.FLAG['IROTAT']=self.IROTAT

        
        """more to follow!"""
        
        all.close()


    def numberAtoms(self,mask='*'):
        """number of Atoms with a particular mask"""
        N=-1
        buf='strip '+mask+'\ngo'
        xxx,scriptf=mkstemp(dir=tmpdir)
        open(scriptf,'w').write(buf)
        xxx,outf=mkstemp(dir=tmpdir)
        cmd=ptrajexe+' '+self.top+' < '+scriptf+' > '+outf+' 2>/dev/null'
        os.system(cmd)
        buf=''.join(open(outf,'r').readlines())
        os.system('/bin/rm '+scriptf+' '+outf) #cleanup
        match=re.compile('Mask.*represents (\d+) atoms').search(buf)
        if match: N=int(match.group(1))
        return N

    def indexes(self,regexp='.'):
        """return indexes of atoms for which their name mathes the
        regular expression regexp"""
        
        
class crd:
    """manage output trajectory coordinate files"""
    def __init__(self,crdf,nat,cformat=''):
        """index the frames in the trajectory
        cformat: field for every coordinate, in printf format"""
        from inputArgs.inputArgs import addAbsPath
        crdf=addAbsPath(crdf)
        if not isthere(crdf):
            sys.stderr.write(crdf+': no such file\n')
            return None
        self.cformat=cformat
        dumpf=joink(crdf+'.dump')
        if isthere(dumpf): #create object from the dump file
            obj=pickle.load(open(dumpf,'r'))
            self.__dict__=deepcopy(obj.__dict__) #no need to initialize
            return None
        self.dumpf=dumpf
        self.crdf=crdf
        self.zip=False
        self.junk=None
        if crdf[-3:]=='.gz': self.zip=True
        self.nat=nat    #number of atoms
        self.comment=None #comment line, if present
        self.hascomment()
        self.deltaidx=0 #jump in file-pointer to next frame
        self.box=None   #bounding box line, if present
        self.get_nlpf_box()
        self.nframe=0
        self.idx=[]      #frame indexing of crd
        self.fillidx()   #index trajectory
        self.frame=None  #current frame read
        self.frames=None
        self.ptf=None
        self.currframe=0
        self.dump()
        
    def dump(self):
        """store object in the production directory"""
        pickle.dump( self,open(self.dumpf,'w'), True )

    def open(self):
        """open self.crd for reading
        self.junk will store uncompress trajectory"""
        if not self.junk or not isthere(self.junk):
            if self.zip:
                xxx,self.junk=mkstemp(dir=tmpdir)
                cmd='/bin/cp '+self.crdf+' '+self.junk+'.gz && '+\
                     'gunzip -f '+self.junk+'.gz'
                os.system(cmd)
                crdf=self.junk  #;  print 'self.junk=',self.junk
            else:
                crdf=self.crdf
        else:
            crdf=self.junk
        return open(crdf,'r')

    def close(self,task=None):
        """close self.crd"""
        if self.junk:
            os.system('/bin/rm '+self.junk)
            self.junk=None
        if task:
            if task=='remove dump file':
                os.system('/bin/rm '+self.dumpf)
        
    def hascomment(self):
        """inquire if first line in file is a comment"""
        line=self.open().readline()
        flag=False
        if re.compile('[a-z]').search(line):
            self.comment=line
            flag=True
        self.close()
        return flag

    def get_nlpf_box(self):
        """find number of lines needed to advance one frame"""
        from utilities.readingWritingFiles import read_to_numpy
        #find how much to advance
        ptf=self.open()
        if self.comment: ptf.readline()
        startpos=ptf.tell()
        arr=read_to_numpy(ptf,3*self.nat,cformat=self.cformat)
        endpos=ptf.tell()
        #is there a box info?
        line=ptf.readline()
        if len(line.split())<7: #assume no more than six columns
            self.box=line
            endpos=ptf.tell()
        #how many positions in between two consecutive frames?
        self.deltaidx=endpos-startpos-1
        
        self.close()

    def fillidx(self):
        """index crd file and count number of frames"""
        ptf=self.open(); ptf.seek(0,0)
        if self.comment: ptf.readline()
        pos=0
        #count number of lines per snapshot
        nline=0; n=0; N=3*self.nat; line=''
        while n<N:
            line=ptf.readline()
            n+=len(line.split())
            nline+=1
        if self.box:
            ptf.readline(); nline+=1
        self.nframe+=1; self.idx.append(pos); pos=ptf.tell()
        #count number of frames
        while line:
            iline=0
            while iline<nline:
                line=ptf.readline(); iline+=1
            if line: self.nframe+=1; self.idx.append(pos); pos=ptf.tell()
        self.close()
        
    def opencrd(self,startframe=None):
        """open trajectory file"""
        self.ptf=self.open()
        if not startframe: startframe=0
        self.currframe=startframe       #current frame index
        self.curr=self.idx[startframe]  #current file-frame index
        self.ptf.seek(self.curr)        #update to startframe

    def rewind(self):
        """rewind trajectory file to start"""
        if self.ptf: self.ptf.close() #close the file pointer
        self.ptf=self.opencrd()   #open file pointer and reset counters
        
    def skipframes(self,skip):
        """skip certain number of frames"""
        #self.curr+=self.deltaidx*skip
        #self.currframe+=skip
        #self.ptf.seek(self.curr)
        for iskip in range(skip): self.loadframe() #cumbersome approach:
        
    def loadframe(self,iframe=None):
        """load frame to memory"""
        from utilities.readingWritingFiles import read_to_numpy
        if not self.ptf: self.opencrd()
        if iframe:
            self.currframe=iframe
            self.curr=self.idx[iframe]
        self.ptf.seek(self.curr)
        if self.curr==0 and self.hascomment: self.ptf.readline()
        frame=read_to_numpy(self.ptf,3*self.nat,cformat=self.cformat)
        if not frame.any(): return False
        if self.box: self.box=self.ptf.readline() #advance box line
        self.curr=self.ptf.tell()
        self.currframe+=1
        self.frame=frame.reshape(self.nat,3)
        return True
        
    def loadframes(self,startframe=None,endframe=None,skip=0):
        """load frames to memory"""
        if not startframe: startframe=self.currframe
        if not endframe: endframe=self.nframe
        self.opencrd(startframe=startframe)
        nframes=(endframe-startframe)/(1+skip)
        N=nframes*self.nat*3
        self.frames=numpy.zeros(N).reshape(nframes,self.nat,3)
        jframe=0           #number of loaded frames
        iframe=startframe  #number of frames read over
        for jframe in range(nframes):
            self.loadframe() #update self.frame
            self.frames[jframe]=self.frame
            self.skipframes(skip) #;print 'jframe=',jframe
        return None
        
    def __iter__(self):
        self.rewind()
        return frame

    def next(self):
        if self.loadframe():
            return self.frame
        return numpy.array(None)

    def nextravel(self):
        """same as next, but flatten the array"""
        if self.loadframe():
            return self.frame.ravel()
        return numpy.array(None)
        
    def avframe(self):
        """calculate average frame"""
        avfr=numpy.zeros(self.nat*3).reshape(self.nat,3)
        while self.loadframe(): avfr+=self.frame
        avfr/=self.nframe
        self.rewind()
        return avfr

    def geometricCenterTraj(self):
        """output trajectory coordinates of geometric Center"""
        traj=[]
        while self.loadframe():
            traj.append( self.frame.sum(axis=0)/self.nat )
        return numpy.array(traj)

    def msd2(self):
        """calculate mean square displacements"""
        av=numpy.zeros(self.nat*3).reshape(self.nat,3)
        dev=numpy.zeros(self.nat*3).reshape(self.nat,3)
        while self.loadframe():
            av+=self.frame; dev+=self.frame*self.frame
        av/=self.nframe
        return dev.sum(axis=1)/self.nframe-(av*av).sum(axis=1)
    
    def findCentroid(self):
        """find Centroid of a trajectory"""
        from tertiary.rmsd import rmsd
        #load maxNfr frames memory...DANGER!
        maxNfr=5000; self.rewind(); frames=[]; iframe=0
        delta=float(self.nframe)/maxNfr;
        jframe=0; next=int(round(jframe*delta))
        while self.loadframe():
            if iframe==next: #signal to record frame
                jframe+=1; next=int(round(jframe*delta))
                frames.append(self.frame);
                #print 'iframe=',iframe,'jframe=',jframe
            iframe+=1;
        maxNfr=jframe
        #proceed to calculate cumulative RMSD for the maxNfr frames
        rmsdT=numpy.zeros(self.nframe)
        for i in range(maxNfr-1):
            iframe=frames[i]
            for j in range(i+1,maxNfr):
                r=rmsd(iframe,frames[j]); #print i,j,r
                rmsdT[i]+=r; rmsdT[j]+=r
        #index of trajectory with lowest cumulative RMSD
        I=numpy.argsort(rmsdT)[::-1][0]
        return False
        
    def varcov(self,masses=numpy.array(None)):
        """calculate (mass-weighted) variance-covariance matrix
        varcov_{ij}=<(x_i-av_i)*(x_j-av_j) >
        x_i: 3*self.nat (mass-weighted) coordinates
        av_i=<x_i>
        """
        from stats.statistics import calc_varcov
        if not masses.all(): masses=numpy.ones(self.nat)
        M=self.nat*3; sm=numpy.sqrt(masses);
        tm=numpy.ones(M).reshape(self.nat,3) #tm speeds up calculations
        for i in range(self.nat): tm[i]*=sm[i]
        av=numpy.zeros(M); dev=numpy.zeros(M)
        vcv=numpy.zeros(M*M).reshape(M,M);
        self.rewind(); iframe=0
        while self.loadframe():           
            fr=(self.frame)*tm; iframe+=1;
            #******************* only for debugging
            #print 'ONLY FOR DEBUGGING iframe=',iframe
            #if iframe==100: break
            #                    only for debugging*******************
            fr=fr.ravel(); av+=fr; dev+=fr*fr
            vcv+=numpy.dot(fr.reshape(M,1),fr.reshape(1,M)) #slow
        results=calc_varcov(av,dev,vcv,M,iframe)
        results['av']=results['av'].reshape(self.nat,3)
        results['dev']=results['dev'].reshape(self.nat,3)
        return results

    def quasiharmonic(self,masses=numpy.array(None),evecspt=None):
        """calculate quasiharmonic modes"""
        from numpy.linalg import eigh,inv
        from utilities.readingWritingFiles import write_from_numpy as wfn
        #variance-covariance mass-weighted matrix
        results=self.varcov(masses=masses)
        #diagonalize and sort from bigger to smaller eigenvalue
        [evals,evecs]=eigh(results['varcov']);
        perm=numpy.argsort(evals)[::-1]
        evals=evals[perm]; evecs=evecs[:,perm]
        results.update({'evals':evals,'evecs':evecs})
        #save to evecs file        
        if evecspt:
            M=3*self.nat
            buf='Eigenvector file: MWCOVAR\n'
            buf+=repr(M)+' '+repr(M)
            evecspt.write(buf)
            wfn(evecspt,evals,format=' %8.5f',ncols=7)
            for ivec in range(M):
                evecspt.write('****\n%4d %11.5f\n'%(ivec+1,evals[ivec]))
                wfn(evecspt,evecs[ivec],format=' %8.5f',ncols=7)
        #calculate force constant matrix by inverting variance-covariance
        #matrix. Remove the last 6 eigenvalues, which should be zero
        #because of translational  and rotational invariance
        I=numpy.diag(numpy.append(1/evals[0:-6],[0,]*6))
        Kij=numpy.dot(evecs,numpy.dot(I,evecs.transpose()))
        results.update({'Kij':Kij})
        return results
        
    def eflux(self,velobj,masses=numpy.array(None)):
        """calculate energy flux between particles in the harmonic
        approximation. We assume that all snapshots in the coordinate
        trajectory have been adjusted to the simulation centroid
        velobj: velocity trajectory, instance of crd class
        eflux_ij=0.5*sum_
        eflux_ij=0.5*kij(v_i+v_j)( (r_i-c_i)-(r_j-c_j) )
        eflux_ij: net flux coming to i from j
        v_i: velocity vector of atom i
        r_i: position vector of atom i
        c_i: centroid position vector of atom i
        """
        from numpy import tensordot as tdt
        nn=3*self.nat;
        #find quasiharmonic modes to obtain the variance-covariance matrix
        print '\n\nWARNING: WE USE SELF-QUASIHARMONIC\n\n'
        results=self.varcov(masses=masses)
        print 'eflux performing inverse of varcov'
        Kij=numpy.linalg.pinv(results['varcov']).reshape(self.nat,3,
                                                         self.nat,3)
        #net and deviation flux
        fnet=numpy.zeros(self.nat*self.nat).reshape(self.nat,self.nat)
        fdev=numpy.zeros(self.nat*self.nat).reshape(self.nat,self.nat)
        #average coordinate position
        frameav=results['av']
        #iterate over the simulation
        self.rewind(); iframe=0
        while self.loadframe():
            iframe+=1;
            print 'eflux iframe=',iframe
            #******************* only for debugging
            #print 'ONLY FOR DEBUGGING iframe=',iframe
            #if iframe==100: break
            #                    only for debugging*******************

            #f:current flux
            f=numpy.zeros(self.nat*self.nat).reshape(self.nat,self.nat)
            #coordinate and velocities
            cfr=(self.frame-frameav).ravel()
            velobj.loadframe(); vfr=velobj.frame.ravel()
            #commuter [x_{ia}v_{jb}-v_{ia}x_{jb}], where "i" and "j"
            #are particle index, "a" and "b" are one of X,Y,Z
            #directions
            comm=numpy.outer(cfr,vfr)-numpy.outer(vfr,cfr)
            comm=comm.reshape(self.nat,3,self.nat,3)
            #compute flux=0.5 sum_{ab} K_{iajb} comm_{iajb}
            
            f=0.5*(Kij*comm).sum(axis=1).sum(axis=-1)
            fnet+=f; fdev+=f*f
        fdev=numpy.sqrt(fdev/iframe-(fnet/iframe))
        abs=f.sum(axis=1)-f.diagonal() #net flux coming to each particle
        return {'fnet':fnet,'absorvance':abs,'fdev':fdev,'nframe':iframe}
        
    def coalesce(self,outpt,groups,masses=numpy.array(None)):
        """create trajectory of macro-particles, for instance, the
        trajectory of the center of mass of the residues
        center of masses coordinates
        groups: list of groups [ [1,4,56], [90,200], [0,6] ]
        masses: weight each atom by its mass"""
        from utilities.readingWritingFiles import write_from_numpy
        ngroups=len(groups)
        cm=numpy.zeros(ngroups*3).reshape(ngroups,3) #center of mass
        cmav=numpy.zeros(ngroups*3).reshape(ngroups,3)
        cmdev=numpy.zeros(ngroups*3).reshape(ngroups,3)
        M=numpy.zeros(ngroups) #masses of the groups
        if not masses.all(): masses=numpy.ones(self.nat)
        while self.loadframe():
            igroup=0
            for group in groups: #go residue by residue
                M[igroup]=0.0
                for index in group: #go atom by atom within current residue
                    cm[igroup]+=masses[index]*self.frame[index]
                    M[igroup]+=masses[index]
                cm[igroup]/=M[igroup]; igroup+=1
            cmav+=cm; cmdev+=cm*cm
            write_from_numpy(outpt,cm,ncols=10,format=' %7.4f',comment='')
        cmav/=self.nframe; cmdev=numpy.sqrt(cmdev-cmav*cmav)
        return {'M':M,'cmav':cmav,'cmdev':cmdev}
    
    def kinetic_energy(self,masses,outptf=None):
        """find kinetic energy of each frame 1/2*mass*x^2"""
        M=masses
        if not outptf: #return K of current frame
            v2=(self.frame*self.frame).sum(axis=1)
            return 0.5*M*v2
        else:          #iterate over all frames
            from utilities.readingWritingFiles import write_from_numpy
            Kav=numpy.zeros(self.nat); Kdev=numpy.zeros(self.nat)
            while self.loadframe():
                v2=(self.frame*self.frame).sum(axis=1)
                K=0.5*M*v2; Kav+=K; Kdev+=K*K
                if outptf:
                    write_from_numpy(outptf,K,format=' %8.4f',ncols=10,
                                     comment='# %05d'%self.nframe)
            Kav/=self.nframe; Kdev=numpy.sqrt(Kdev/self.nframe-Kav*Kav)
            return (Kav,Kdev) #return averages and standard deviation

    def Kcov(self,masses):
        """calculate covariance matrix of kinetic energies"""
        from numpy.linalg import eigh
        from stats.statistics import calc_cov
        #average and standar deviation for kinetic energy
        Kav=numpy.zeros(self.nat); Kdev=numpy.zeros(self.nat)
        #covariance of the kinetic energy
        Kcov=numpy.zeros(self.nat*self.nat).reshape(self.nat,self.nat)
        while self.loadframe():
            v2=(self.frame*self.frame).sum(axis=1);
            K=0.5*masses*v2; Kav+=K; Kdev+=K*K
            Kcov+=numpy.dot(K.reshape(self.nat,1),K.reshape(1,self.nat))
        results=calc_cov(Kav,Kdev,Kcov,self.nat,self.nframe)
        [evals,evecs]=eigh(results['cov']);        #diagonalize
        perm=numpy.argsort(evals)[::-1] #sort from bigger to smaller
        evals=evals[perm]; evecs=evecs[:,perm]
        return {'Kav':results['av'], 'Kdev':results['dev'],
                'Kcov':results['cov'], 'evals':evals, 'evecs':evecs}
    
    def dK_dt(self,masses,outptf=None):
        """find derivatives of the kinetic energy"""
        from utilities.readingWritingFiles import write_from_numpy
        self.rewind(); Kd=None
        dKav=numpy.zeros(self.nat); dKdev=numpy.zeros(self.nat)
        self.loadframe(); curr=self.kinetic_energy(masses); iframe=0;
        while self.loadframe():
            iframe+=1; next=self.kinetic_energy(masses)
            Kd=next-curr; dKav+=Kd; dKdev+=Kd*Kd
            write_from_numpy(outptf,Kd,format=' %8.4f',ncols=10,
                             comment='# %05d'%iframe)
            curr=next;
        #last frame assume same Kd as next to last frame
        iframe+=1; dKav+=Kd; dKdev+=Kd*Kd
        write_from_numpy(Kdpt,Kd,format=' %8.4f',ncols=10,
                         comment='# %05d'%iframe)
        #return averages and standard deviations
        dKav/=iframe; dKdev=numpy.sqrt(dKdev/iframe-dKav*dKav)
        return (dKav,dKdev)

    def kflux(self,masses,VdW,neighbors,distances,xyz,outptf):
        """compute flux of kinetic energy
        masses: array of masses
        VdW: array of Vaan der Walls radii"""
        from utilities.readingWritingFiles import write_from_numpy
        outptf.write('trajectory generated by amber10\n')
        Kd=None; VdW2=VdW*VdW #square of the VdW radii
        flav =numpy.zeros(self.nat*3).reshape(self.nat,3);#average flux
        fldev=numpy.zeros(self.nat*3).reshape(self.nat,3);#stdev flux
        fl   =numpy.zeros(self.nat*3).reshape(self.nat,3);#current flux
        #calculate weighted relative distance matrix and extended VdW
        rel=numpy.zeros(self.nat*self.nat*3).reshape(self.nat,self.nat,3)
        for i in range(self.nat):
            r=xyz[i]; neig=neighbors[i]; dists=distances[i]
            for j in range( len(neig) ):
                k=neig[j]; d=dists[j];
                rel[i][k]=(xyz[i]-xyz[k])/(d*d*d)
        rel=rel.transpose( (0,2,1) ) #rearrange indexes
        VdWx=numpy.ones(self.nat*3).reshape(self.nat,3);
        for i in range(self.nat): VdWx[i]*=VdW[i]
        #open velocity file
        self.opencrd();
        self.loadframe();
        curr=self.kinetic_energy(masses);
        iframe=0;
        while self.loadframe():
            iframe+=1; print iframe
            next=self.kinetic_energy(masses)
            Kd=next-curr;
            fl=VdWx*numpy.dot(rel,Kd)
            flav+=fl;
            fldev+=fl*fl
            write_from_numpy(outptf,fl.ravel(),format=' %7.3f',ncols=10)
        #last frame
        flav+=fl;
        fldev+=fl*fl
        write_from_numpy(outptf,flav.ravel(),format=' %7.3f',ncols=10)
        #return net flux (not average) and standard devations
        fldev=numpy.sqrt(fldev/iframe-flav*flav/(iframe*iframe))
        return (flav,fldev,iframe)
    
    def distances(self):
        """calculate all pairs of distances"""
        dij=numpy.zeros(self.nat*self.nat).reshape(self.nat,self.nat)
        for i in range(self.nat):
            r=frame[i]
            for j in range(self.nat):
                rel=r-frame[j];
                dij[i][j]=numpy.sqrt( (rel*rel).sum() )
        return dij

    def writeFrame(self,fr=numpy.array([])):
        """write a frame and the bounding box line, if present"""
        if not fr.any(): fr=self.frame
        fr=list( fr.ravel() )
        buf=''
        while fr:
            buf+=''.join(['%8.3f'%x for x in fr[0:10] ])+'\n'
            fr=fr[10:]
        if self.box: buf+=self.box
        return buf
            
    def translate(self,outf,t=numpy.array([])):
        """translate trajectory by some conformation
        If empty, we translate by minus the average conformation"""
        pt=open(outf,'a') #open for append
        buf='trajectory generated by amber10.py\n'
        L=50E4  #L; maximum buffer size
        self.rewind()
        if not t.any(): t=-self.avframe()
        read=0
        while self.loadframe():
            self.frame+=t
            buf+=self.writeFrame()
            if len(buf)>L:
                pt.write(buf); buf=''
            read+=1
            #for debugging purposes
            #print 'translate %d'%read
            #if read==10000:break
        if buf: pt.write(buf)
        self.rewind()
        
    def removeMode(self,evecsobj,imode,outf=None):

        """remove projection of mode imode from trajectory and store
        in new trajectory nextcrdf. All neccessary info for the mode
        in evecsobj, an instance of class evecs

        NOTE: we previously must have removed the average conformation from
        the trajectory. This can be done using method self.removeAvg and
        average conformation evecsobj.avgX"""
        pt=open(outf,'a') #open for append
        buf='trajectory generated by amber10.py\n'
        L=50E4 #maximum buffer size
        Q=(evecsobj.qs[imode-1]).m.ravel() #quasiharmonic mode
        sM=evecsobj.sM
        isM=evecsobj.isM
        read=0
        while self.loadframe():
            fr=self.frame.ravel()
            q=(fr*sM*Q).sum()  #projection of the mode
            fr-=sM*Q*q         #remove projection from the frame
            buf+=self.writeFrame(fr=fr)
            if len(buf)>L:
                pt.write(buf); buf=''
            read+=1
            #print 'removeMode %d'%read
        if buf: pt.write(buf)
        self.rewind()
        
class sander:
    """dealing with sander input and output"""

    def __init__(d=None):
        
        self.exe=None  #path to executable
        self.inpf=None #path to input file
        self.topf=None #path to topology file
        self.crdf=None #path to output crd file
        self.crd=None  #class crd object
        self.velf=None #path to output vel file
        self.vel=None  #class vel object
        if d: self.__dict__.update(d)


    def initialize(self,attr):
        if attr=='vel':
            if not self.velf: sys.stderr.write('ERROR: no velocity file')
            self.vel=crd(self.velf,self.nat)
            


class ptraj:
    """dealing with ptraj"""

    def __init__(self,top,crd,vel=None):
        """crd can be one string or a list

        top   (str) topoly file
        cmd   (list)
        crds  (list)
        Ncs   (list, number of frames in every coordinate file)
        Nc    (total number of frames)
        vel
        outcrd output trajectory file
        """
        self.cmd=[]   #list of commands
        self.crds=[]  #list of coordinate files
        self.Ncs=[] 
        self.Nc=0   
        self.top=top
        if type(crd)==str:
            self.crds=[crd,]
        elif type(crd)==list:
            self.crds=crd
        self.crd=self.crds[0]
        if vel: self.vel=vel
        
    def iscommitted(self,command):
        """check if command keyword is already in self.cmd"""
        buf='\n'+'\n'.join(self.cmd)
        if '\n'+command in buf:
            sys.stderr.write(command+' already committed\n')
            return True
        return False
        
    def trajin(self,start=None,end=None,offset=None):
        #check first we did not do trajin first
        if self.iscommitted('trajin'): return {'self':self}
        tmp=[]
        for crd in self.crds:
            buf='trajin '+crd
            if start: buf+=' '+`start`
            if end: buf+=' '+`end`
            if offset: buf+=' '+str(offset)
            tmp.append(buf)
        self.cmd=tmp+self.cmd #trajin goes to the beginning
        return {'self':self}

    def trajout(self,outcrd,args=''):
        if self.iscommitted('trajout'): return {'self':self}
        self.cmd.append('trajout '+outcrd+' '+args)
        return {'self':self}

    def strip(self,entities=[':WAT',]):
        """remove entities list"""
        for entity in entities: self.cmd.append('strip '+entity)
        return {'self':self,}

    def reference(self,reffile):
        """set reference structure"""
        self.cmd.append('reference '+reffile)
        return {'self':self}

    def rms(self,args=None,mode='first',outrms=''):
        if args:
            self.cmd.append('rms '+args)
        else:
            self.cmd.append('rms '+mode+' out '+outrms)
        return {'self':self}

    def atomicfluct(self,args=None):
        """atomic fluctuations"""
        if args:
            self.cmd.append('atomicfluct '+args)
        return {'self':self}
    
    def correl(self):
        """add correlation matrix"""
        self.cmd.append('matrix correl name correl out correl.dat byres mass')
        return {'self':self}
    
    def quasiharmonic(self,nvecs=1000):
        """output quasiharmonic analysis"""
        list=['matrix mwcovar name mwcovar out mwcovar.dat byatom',
              'analyze matrix mwcovar out evecs.dat thermo vecs '+`nvecs`]
        if not self.iscommitted('rms'):list=['rms first out rms.dat',]+list
        self.cmd+=list
        return {'self':self}

    def go(self,script=None,outfile=None,outcrd=None,fake=None):
        """execute the stored commands
        outfile: output from ptraj
        """
        output=None
        rmscript=False
        rmoutfile=False
        if not script:
            xxx,script=mkstemp(dir=tmpdir)
            rmscript=True
        if not outfile:
            xxx,outfile=mkstemp(dir=tmpdir)
            rmoutfile=True
        self.trajin()
        if outcrd and not self.iscommitted('trajout'):
            self.outcrd=outcrd
            self.trajout(outcrd)
        open(script,'w').write('\n'.join(self.cmd))
        cmd=ptrajexe+' '+self.top+' < '+script+' &> '+outfile
        if fake:
            print cmd
        else:
            os.system(cmd)
            output=''.join(open(outfile,'r').readlines())
            if rmscript: os.system('/bin/rm '+script)
            if rmoutfile: os.system('/bin/rm '+outfile)
        return {'self':self,'output':output}

    def saveSlice(self,start,end,out):
        """save a slice of the first coordinate file"""
        traj=ptraj(self.top,self.crds[0]).trajin(start,end).trajout(out)
        
    def count_frames(self,crd=None):
        if not crd: crd=self.crds
        if type(crd).__name__!='list': crd=[crd,] #make it a list
        self.Ncs=[]
        for x in crd:
            output=''.join(ptraj(self.top,x).go()['output'])
            pattern=re.compile('Successfully read in (\d+) sets')
            match=pattern.search(output)
            self.Ncs.append( int(match.group(1)) )
        self.Nc=sum(self.Ncs)
        return {'self':self,'N':self.Nc}            
        
    def create_rst(self,rst,offset=0,vel=None):
        """create restart file from crd and vel files"""
        if not vel: vel=self.vel
        Nx=self.count_frames()['N']         #;print 'Nx='+`Nx`
        Nv=self.count_frames(crd=vel)['N']  #;print 'Nv='+`Nv`
        N=min(Nx,Nv)-offset                 #;print 'N ='+`N`
        if N<=0: return {'self':self,'N':0} #empty files or some error
        #do a restart of coordinates
        xxx,restart_crd=mkstemp(dir=tmpdir)
        cmd='trajin '+self.crd+' '+`N`+' '+`N`+'\n'+\
             'trajout '+restart_crd+' restart'
        #print cmd
        xxx,script=mkstemp(dir=tmpdir)
        open(script,'w').write(cmd)
        cmd=ptrajexe+' '+self.top+' < '+script+' &>/dev/null '+\
             ' && /bin/rm '+script
        os.system(cmd)
        restart_crd+='.'+`N`   #annoying twist by ptraj
        #do a restart of velocities
        xxx,restart_vel=mkstemp(dir=tmpdir)
        cmd='trajin '+vel+' '+`N`+' '+`N`+'\n'+\
             'trajout '+restart_vel+' restart\n'
        #print cmd
        xxx,script=mkstemp(dir=tmpdir)
        open(script,'w').write(cmd)
        cmd=ptrajexe+' '+self.top+' < '+script+' &>/dev/null '+\
             ' && /bin/rm '+script
        os.system(cmd)
        restart_vel+='.'+`N`
        #merge both restart files
        box=os.popen('tail -1 '+restart_crd).readline()  #;print 'box=',box
        buf=''
        pt=open(restart_crd,'r') ; line=pt.readline()
        while '   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000' not in line:
            buf+=line ; line=pt.readline()
        pt=open(restart_vel,'r')
        pt.readline() ; pt.readline() ; line=pt.readline() ;
        while '   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000' not in line:
            buf+=line ; line=pt.readline()
        buf+=box
        open(rst,'w').write(buf)
        os.system('/bin/rm '+restart_crd+' '+restart_vel)
        return {'self':self,'Nc':Nx,'Nv':Nv,'N':N}

    def cluster(self,args='out junk representative pdb '+\
                'average pdb means clusters 5 rms @CA'):
        """cluster trajectory""" 
        self.cmd.append('cluster '+args)
        return {'self':self}

class evecs:
    """modes as given by amber format, output of analyze matrix
    Given the Linear equation:
                              C E = M L^2 E
    with mass-weighted covariance matrix C = <X^t M X> and diagonal
    matrix L.
    The eigenvalue matrix E diagonalizes both the mass matrix and the
    mass weighted covariance matrix:
                              E^t M E = I
                              E^t C E = L^2
    Thus, the units of the eigenvector components are M^(-1/2)
    The reported eigenvalues in evecs file (in units of 1/cm)
    correspond to the diagonal (M L^2) matrix, not L^2. 
    """

    def __init__(self,evecsf):
        self.evecsf=evecsf   #evecs file
        self.q=None      #quasiharmonic_mode class instance
        self.pt=None     #pointer to evecs file
        self.iq=0        #current mode
        self.nq=None     #number of modes
        self.qs=numpy.array([])  #quasiharmonic modes
        self.qspf=None   #file of projections of quasiharmonic modes
        self.find_nq()   #find number of modes in evecsf
        self.freqs=numpy.array([])#numpy array of mode frequencies
        self.M=numpy.array([])    #mass of each degree of freedom
        self.sM=numpy.array([])   #sqrt(M)
        self.isM=numpy.array([])  #1/sqrt(M)
        self.natoms=0             #number of atoms
        self.avgX=numpy.array([]) #average coordinates shape=(natoms,3)
        self.load_avgX()

    def open(self):
        if not isthere(self.evecsf):
            sys.stderr.write(self.evecsf+': No such file. Nothing to do')
        return open(self.evecsf,'r')

    def close(self):
        """allows pickling of self"""
        if self.pt: self.pt.close()
        self.pt=None
        return None

    def write2string(self):
        """write object to string"""
        from lib.utilities.readingWritingFiles import write_from_numpy
        self.GetFreqs()
        buf='Eigenvector file: MWCOVAR\n'
        buf+='%5d %5d\n'%(self.nq, self.nq)
        buf+=write_from_numpy(None,self.avgX.ravel(),
                              format='%11.5f',ncols=7)
        for i in range(self.nq):
            buf+=' ****\n'
            buf+='%5d %12.5f\n'%(1+i,self.freqs[i])
            buf+=write_from_numpy(None,self.qs[i].m.ravel(),
                                  format='%11.5f',ncols=7)
        return buf

    def massProcess(self,mlist):
        """obtain several mass vectors from the list of atomic
        masses. It is suggested you pass the mass list from a 'top'
        object."""
        self.M=[]
        for m in mlist: self.M+=[m,m,m]
        self.M=numpy.array(self.M)
        self.sM=numpy.sqrt(self.M)
        self.isM=1.0/self.sM
        return self.M

    def find_nq(self):
        """find number of modes"""
        cmd='grep -A 1 -e "\*\*\*" %s|tail -1'%(self.evecsf,)
        self.nq=int(os.popen(cmd).readline().split()[0])
        return self.nq

    def load_avgX(self):
        pt=self.open()
        pt.readline()
        z=[int(x) for x in pt.readline().split()][0]#degrees of freed
        self.natoms=z/3
        self.avgX=[]; read=0
        while len(self.avgX)<z:
            self.avgX+=[float(x) for x in pt.readline().split()]
        self.avgX=numpy.array(self.avgX).reshape(self.natoms,3)
        self.close()
            
    def loadModes(self,nq=None):
        """load all modes into memory
        nq: load only the first nq modes"""
        if not nq: nq=self.nq
        pt=self.open()
        self.qs=[] ; read=0
        while self.readMode():
            self.qs.append(self.q) #append current mode
            read+=1
            if read==nq: break
        self.nq=len(self.qs)
        self.qs=numpy.array(self.qs)

    def normalize_modes(self,masses=[]):
        """normalize all modes"""
        if not masses:
            if self.M.any(): masses=self.M
        if not self.qs.any(): self.loadModes()
        for i in range(self.nq):            
            self.qs[i].m=(self.qs[i].m.ravel()*sM).reshape(self.qs[i].N,3)

    def return_normal_modes(self,masses=None):
        if not masses:
            if self.M.any(): masses=self.M
        sM=numpy.sqrt(masses)
        if not self.qs.any(): self.loadModes()
        Q=[]
        for i in range(self.nq):
            q=deepcopy(self.qs[i])
            q.m = ( q.m.ravel() * sM ).reshape(q.N,3)
            Q.append( q )
        return Q
    
    def readMode(self):
        """read one mode"""
        if not self.pt: self.pt=self.open()
        self.q=quasiharmonic_mode(self.pt)
        if self.q.error: return None
        self.iq+=1; #print 'iq=%d'%self.iq
        return self.q

    def GetFreqs(self):
        """fill self.freqs"""
        fr=[]
        if not self.freqs.any():
            if self.qs.any():
                for i in range(self.nq): fr.append(self.qs[i].e)
            self.freqs=numpy.array(fr)
        return self.freqs

    def HessSimilarity( self, evB, MA=[], MB=[], mA=-1, mB=-1,
                        RMS=False, lA=[], lB=[]):
        """compare self to another set of modes with the Hess metric:
        Hess = 1/sqrt(nA*nB) Sum_i^nA Sum_j^nB qA_i * qB_j, where:
         - qA_i is ith normal mode from self (qB_j normal mode from evB)
         - nA modes from self and nB modes from evB. We apply Hess
           varying nA from 1 to mA and nB from 1 to mB. mA=-1 means use
           all modes from self (mB=-1 use all evB modes)
         - MA, MB: list of masses for normalization of modes. Default:
           MA=self.M, MB=evB.M
         - RMS: by default, we do not do an RMS superimposition of the
           two average structures. If RMS is set to True, then:
           - find RMS superposition of atoms of evB onto all
             atoms of self, using the coordinates from the corresponding
             average structures. Then, rotate modes of evB using the
             rotation resulting of the RMS superposition.
           - By default, all atoms are selected. If lA and lB are
             supplied, use only the atoms specified in lA and lB for
             the RMS superposition and Hess Metric( len(lA) should be
             same as len(lB) ). lA is a list of atom indexes, with
             first atom having index=0 (analogously for lB).
         - outf: if passed, print results in a buffer string.
        """
        #find out masses
        if not MA: MA=self.M
        if not MB: MB=evB.M
        #maximum number of modes to use
        if mA<0: mA=self.nq
        if mB<0: mB=evB.nq
        #normalize modes
        qsA=self.return_normal_modes(masses=MA)
        qsB=evB.return_normal_modes(masses=MB)

        buf='#Hess similarity metric\n'
        #Do we need to do RMS superposition of qsB to qsA ?           
        if RMS:           
            from tertiary.rmsd import rmsd
            xA = self.avgX
            xB = evB.avgX
            if lA:
                xA = []
                for i in lA: xA.append( self.avgX[i] )
                xA = numpy.array( xA )
                xB = []
                for i in lB: xB.append( evB.agvX[i] )
                xB = numpy.array( xB)
            r,R,t = rmsd( xA, xB )
            buf += '#rot matrix:\n#%s\n' % str(R).replace( '\n', '\n#' )
            #rotate each mode of qsB according to R
            for i in range(mB):
                qsB[i].rotate( R.transpose() )
        #Absolute value of dot product
        D=numpy.zeros(mA*mB).reshape(mA,mB)        
        for I in range(mA):
            for J in range(mB):
                D[I][J]=abs(numpy.dot(qsA[I].m.ravel(),qsB[J].m.ravel()))
        #cummulative summ on the rows
        Q=numpy.zeros(mA*mB).reshape(mA,mB)   
        for I in range(mA):
            h = numpy.dot( qsA[I].m.ravel(),qsB[0].m.ravel() )
            Q[I][0] = h*h
            for J in range(1,mB):
                h = numpy.dot( qsA[I].m.ravel(),qsB[J].m.ravel() )
                Q[I][J] = h*h + Q[I][J-1]
        #cummulative summ on the columns
        Hs=numpy.zeros(mA*mB).reshape(mA,mB)
        for J in range(mB): Hs[0][J] = Q[0][J]
        for I in range(1,mA):
            for J in range(mB): Hs[I][J] = Q[I][J] + Hs[I-1][J]
        #normalize
        for I in range(mA):
            print 1+I,1+I,Hs[I][I]/float(1+I)
            for J in range(mB):
                Hs[I][J] /= numpy.sqrt((1+I)*(1+J))
                buf +=  '%4d %4d %6.4f\n' % (1+I,1+J,Hs[I][J])
        return D,Hs,buf
                
    def rotate_modes(self,R):
        """rotate modes according to rotation R"""
        for i in range(self.nq): self.qs[i].rotate( R )
        return None

    def align_modes_to(self, evB,):
        """rotate modes according according to rotation R resulting
        from superposition of average structure to that of evB
        """
        from tertiary.rmsd import rmsd,euler2rot
        r,R,t = rmsd( self.avgX, evB.avgX )
        W = euler2rot( (20*numpy.pi/180)*numpy.random.rand(3))
        R = numpy.dot(W,R)
        self.rotate_modes( R )
        return r,R,t
        
    def next(self):
        val=self.readMode()
        if val: return self.q
        raise StopIteration
        
    def __iter__(self):
        return self

    def create_PDBs(self, templatePDB, nmodes=None,
                    outdir=None, byres=False):

        """create PBD files for each file in evecs.dat. Each atom
        (residue) will have a b-factor equal to the absolute value of
        its mode coordinate, averaged over X,Y,Z directions, and
        scaled by the corresponding eigenvalue"""
        
        from mypdb.smallThings import insertBfact
        if not nmodes: nmodes = self.nq
        if not outdir: outdir = '.'
        modes = self.return_normal_modes()[ :nmodes ]
        for i in range( nmodes ):
            #displacenment moduli in 3D
            v = numpy.sqrt( (modes[ i ].m* modes[ i ].m).sum(axis=1) )
            #mutated residues, in the ref15-shel sequence numbering
            mutated='1,3,4,6,7,13,14,15,16,23,24,28,30,'+\
                     '31,32,40,43,44,46,47,50,51'
            expected = float(len(mutated)) / len(v)
            #display some enrichment statistics:
            mean= v.mean()
            var = sqrt( v.var() )
            n1 = [] # number of residues avobe mean+var/2
            n2 = [] # number of residues avobe mean+var
            for j in range( len(v) ):
                if v[ j ] > ( mean + var/2.0 ) : n1.append(str(1+j))
                if v[ j ] > ( mean + var ) : n2.append(str(1+j))
            en1 = 0 #number of mutated residues in n1
            for x in n1:
                if x in mutated: en1 += 1
            observed = float( en1 ) / len(n1)
            en1 = 100 * (observed / expected )
            en2 = 0 #number of mutated residues in n2
            for x in n2:
                if x in mutated: en2 += 1
            en1
            observed = float( en2 ) / len(n2)
            en2 = 100*( observed / expected )#enrichment
            #mode avg     var  >avg+var/2  >avg+var  en1  en2
            print '%3d %4.5f %4.5f    %2d         %2d   %5.1f %5.1f'\
                  %(1+i, mean, var, len(n1), len(n2), en1,en2)
            #output to file
            buf = insertBfact( templatePDB, v, byres=byres )
            outf = '%s/%04d.pdb' % ( outdir, 1+i )
            open( outf, 'w' ).write( buf )
        return False

    def gen_movies(self, templatePDB, nmodes=None, outdir=None,
                   nframe=20, byres=False, maxdispl=0.9, mindispl=0.1):
        """Create PDB movies for each mode in evecs.dat
        nmodes=self.nq create movies for all modes (self.nq)
        outdir=? directory to put the PDB movies
        templatePDB must contain as many atoms as evecs

        maxdispl: cut-off so that pymol will show up the secondary
        structure. 0.4 seem to be the highest cut-off."""
        
        if not nmodes: nmodes=self.nq
        #determine output directory
        currdir = os.getcwd()
        if not outdir: outdir = currdir
        os.system( '/bin/mkdir -p ' + outdir )
        #load coordinates
        from Bio.PDB.PDBParser import PDBParser
        from Bio.PDB.Structure import Structure
        from Bio.PDB import PDBIO
        p = PDBParser( PERMISSIVE=1 )
        model0 = p.get_structure( "tmpl", templatePDB )[ 0 ]
        #iterate over the normalized modes
        if not self.qs.any(): self.loadModes()
        for imode in range( 0,nmodes,10 ):
            print 'imode=',imode
            q = self.qs[ imode ]
            z = sqrt( abs( q.e ) )
            #check maximum atomic displacement
            md = numpy.max( z * numpy.sqrt((q.m*q.m).sum(axis=1)))
            if md > maxdispl: z *= maxdispl/md
            if md < mindispl: z *= mindispl/md
            buf = ''
            st = Structure( id=str( imode ) )
            prevmodel = model0
            for iframe in range( nframe ):
                model = deepcopy( prevmodel )
                model.id = iframe
                x = z * cos( pi/2 + ( 2 * pi * iframe ) / nframe )
                iat = 0  #residue counter
                ires = 0 #atom counter
                for residue in model.get_residues():
                    for atom in residue:
                        if byres: v = q.m[ ires ]
                        else: v = q.m[ iat ]
                        xyz = numpy.array(atom.get_coord())
                        xyz += x*v
                        atom.coord = xyz
                        iat += 1
                    ires += 1
                st.add( model ) #insert model into structure
                prevmodel = model
            io = PDBIO()
            io.set_structure( st )
            io.save( '%s/%04d.pdb' % ( outdir, 1+imode ) )
        return False
        
    def create_movies(self,topfile,nmodes=None,outdir=None):
        """NOTE: legacy_code. Left for compatibility
        Create PDB movies for each mode in evecs.dat
        nmodes=self.nq create movies for all modes (self.nq)
        outdir=? directory to put the PDB movies"""
        if not nmodes: nmodes=self.nq
        currdir=os.getcwd()
        if not outdir: outdir=currdir
        os.system('/bin/mkdir -p '+outdir)
        xxx,workdir=mkdtemp(dir=currdir)
        arvindfiles=joink(codedir,'software/amber9/arvind_movies/*')
        os.system('mkdir -p '+workdir)
        cmd= 'cd '+workdir+' && /bin/cp '+arvindfiles+' . && '+\
             'ln -s '+topfile+' topology ; '+\
             './vecs_to_binpos1.0 < '+self.evecsf+' > binpos 2> junkX ; '+\
             'tail -1 junkX|tr -s " "|awk \'{print $6}\' > bytes'
        os.system(cmd)
        bytes=open(joink(workdir,'bytes'),'r').readline().strip()
        cmd= 'cd '+workdir+' && split -a 3 -b '+bytes+' binpos && '+\
             './doit '+`nmodes`+' '+outdir+' && '+\
             'cd '+currdir+' && '+\
             '/bin/rm -rf '+workdir
        os.system(cmd)
        #create log files for every movie
        template='cmd.load("_MOD_.pdb")\ncmd.dss("_MOD_")\ncmd.show_as("cartoon"   ,"_MOD_")\ncmd.spectrum("count",selection="(_MOD_)&e. c")\ncmd.mplay()'
        for i in range(nmodes):
            mod='mod%04d'%(i+1); logf=outdir+'/'+mod+'.pml'
            open(logf,'w').write(template.replace('_MOD_',mod))
        
    def projectCRD(self,crdobj,outf,nq=None):

        """Find the projection of each trajectory snapshot onto the
        quasi-harmonic modes from an evec file. Format is same as
        'projection' command in ptraj.
        
        outf: output file to store projections
        nq: find projections only on the first nmodes
        masses: list of atomic masses from the topology file
        """
        
        if not nq: nq=self.nq
        pt=None
        if isthere(outf): #check outf is new file
            sys.stderr.write('file %s exists: overwrite?:(y/n)'%(outf,))
            if sys.stdin.readline()[0].lower()!='y':
                sys.stderr.write('OK. Nothing to do!\n')
                return False
        pt=open(outf,'a') #open for append
        iframe=0; buf=''; bufmaxsize=int(50E4)#50E6);#50MB
        #Write file header
        buf='Projection of snapshots onto modes\n  Snapshot'
        for i in range(1,nq+1):
            word='Mode%d'%i; buf+='%10s'%word
        buf+='\n'
        #Write projections for each snapshot
        while crdobj.loadframe():
            iframe+=1; buf+='#frame %05d\n'%(iframe,) #store number
            fr=crdobj.frame.ravel()-self.avgX.ravel()
            for q in self.qs[0:nq]:
                buf+='%10.3f'%(fr*self.sM*q.m.ravel()).sum()
            buf+='\n'
            if len(buf)>bufmaxsize:  #flush buffer
                pt.write(buf); buf=''
        if buf: pt.write(buf)
        return True
    
class mode:
    import numpy
    """mode, normal or quasi-harmonic"""
    def __init__(self,l):
        
        self.error='' #error messages
        self.i=None #mode number
        self.e=None #mode eigenvalue
        self.m=None #store mode as numpy array (N,3)
        self.n=None #store normalized mode as numpy array (N,3)

    def __mul__(self,other):
        """overloading '*' multiplication"""
        return (self.m.ravel()*other.m.ravel()).sum()    

    def norm(self):
        from math import sqrt
        return sqrt(self*self)

    def normalize(self,masses=numpy.array(None)):
        """normalize the vector"""
        norm=self.norm()
        self.n=self.m/norm
        if masses.all(): self.n=self.n*numpy.sqrt(masses)
        return norm

class normal_mode(mode):
    """normal mode, for harmonic systems"""
    pass
    
class quasiharmonic_mode(mode):

    """quasiharmonic mode from ptraj, the result of the following linear
    transformation:
                       X = M^{1/2} * E * Q
    M: mass matrix (3*nparticles x 3*nparticles)
    X: coordinates vector (3*nparticles)
    Q: quasiharmonic modes coordinate vector (3*nparticles)
    E: E diagonalizes the mass weighted covariance matrix
                   varcov = E * evals * E^t
                       E^t * M * E = I
       its columns are the quasi-harmonic eigenvectors
    """
    def __init__(self,pt):
        """pt is a file pointer to and evecs file
        self.u #Mass normalized mode numpy array (N,3)"""
        import numpy
        self.error=''
        l=pt.readline()
        while l and '***' not in l: l=pt.readline() #beginning of new mode
        if not l:
            self.error='reached end of file\n'
            return
        self.i,self.w=[ float(x) for x in pt.readline().split() ]
        self.i=int(self.i)
        #"0.6" is conversion of kT for 300K into kcal/mol
        #from 1/cm to sqrt(Kgr)*Angstrom
        self.  e = 0.6 * ( 108.587 / self.w )**2
        if self.w < 0: self.e = -self.e
        c=[]; curr=pt.tell(); l=pt.readline()
        while l and '***' not in l:
            c+=[float(x) for x in l.split()]
            curr=pt.tell(); l=pt.readline()
        if '***' in l: pt.seek(curr) #rewind one line
        self.N=len(c)/3
        self.m=numpy.array(c).reshape(self.N,3)

    def Massnormalize(self,masses):
        """normalize the vector
        masses is a self.N vector"""
        from numpy import sqrt
        norm=self.norm()
        self.u=self.m/norm
        self.u=self.u*sqrt(masses)
        return norm
        
    def massmul(self,other,masses):
        a=(self.m*other.m).sum(axis=1); N=len(a)
        return (masses[0:N]*a).sum()

    def rotate(self,R):
        """rotate the modes by the provided rotation matrix R"""
        if R.shape==(3,3):
            self.m=numpy.dot( R,self.m.transpose() ).transpose()
        elif R.shape==(3*self.N,3*self.N):
            self.m=numpy.dot( R, self.m.ravel() ).reshape(self.N,3)
            
class evecs_projection:

    """Object for the ouput of 'projection' tool of 'ptraj'
    NOTE: INDEXES IN ITEMS AND SLICES BEGIN WITH ONE, NOT ZERO
                    Q = E^t * M^{1/2} * X
    For mode 'i', projection Q_i = e_i^t * M^{1/2} * X, where e_i is
    the ith quasi-harmonic eigenmode stored in the evecs file, here
    represented as the transpose of a vector of length 3*nparticles.
    """
    def __init__(self,projf):
        """
        projf       #projections file
        pin         #pointer to projf
        imodes      #list of mode numbers
        nsnapshots  #number of frames, or snapshots
        index_projf #index of the snapshots
        stats       #several statistics from the projections
         |_mean     #list of means of q-values
         |_std      #list of standar deviation of the q-values
         |_histogram #list of histograms of the q-values
         |_xvals    #list of q-values for each of the previous histograms
         |_chi      #list of fit of each histogram to a Normal distribution
        """
        self.projf=projf   #projections file
        if '.gz' in projf: self.pin=gzip.open(projf,'r')
        else:              self.pin=open(projf,'r');
        self.pin.readline()#this is a comment
        self.imodes=[int(x) for x in \
                     re.compile('Mode(\d+)').findall(self.pin.readline())]
        self.nmodes=len(self.imodes)
        self.nsnapshots=0
        self.index_projf() #index the snapshots
        self.stats={}

    def close(self):
        """allows pickle saving"""
        if self.pin: self.pin.close()
        self.pin=None
        
    def index_projf(self):
        """Find indexes for every snapshot"""
        self.idx=[0,self.pin.tell()]
        n=0; l=self.pin.readline()
        while l:
            self.idx.append(self.pin.tell()); l=self.pin.readline(); n+=1
        self.pin.seek(0)
        self.nsnapshots=n
        
    def next(self):
        l=self.pin.readline()
        if l:
            return [float(x) for x in l.split()[1:]]
        raise StopIteration

    def __iter__(self):
        self.pin.seek(0)
        self.pin.readline()#this is a comment
        self.pin.readline()#this is list of mode numbers
        return self

    def load_snapshot(self,i):
        """go to appropriate snapshot and load it"""
        self.pin.seek( self.idx[i] )
        l=self.pin.readline()
        return [float(x) for x in l.split()[1:]]
    
    def __getitem__(self,i):
        """return snapshot given by index"""
        return self.load_snapshot(i)

    def __getslice__(self,i,j):
        """return a slice from indexes i to j-1"""
        if i==0: i=1 #snaphots begin at index 1, not zero
        if j>self.nsnapshots: j=self.nsnapshots
        if j<0: j=self.nsnapshots-j
        self.pin.seek( self.idx[i] )
        z=[]        
        for k in range(i,j):
            l=self.pin.readline()
            z.append( [float(x) for x in l.split()[1:]] )
        return z

    def fit2Normal(self):
        """calculate fitting of each projection to a normal
        distribution"""
        from stats.histograms import histogram,chiSquare
        r={'mean':[],'std':[],'chi':[],'histogram':[],'xvals':[]}
        for imode in self.imodes:
            sys.stdout.write('..'+str(imode))
            Qs=[]
            for snapshot in self: Qs.append( snapshot[imode-1] )
            Qs=numpy.array(Qs)
            mean=Qs.mean()
            std=Qs.std()
            h,xvals=histogram(Qs,bins=100)
            chi=chiSquare(Qs,xvals,h,lowcut=0.05)
            r['mean'].append(mean)
            r['std'].append(std)
            r['chi'].append(chi)
            r['histogram'].append(h)
            r['xvals'].append(xvals)
        self.stats.update(r)
        self.close()
        return r
            
class qhm_proj:
    """For every snapshot, find the projection of each trajectory
        snapshot onto the quasi-harmonic modes"""
    def __init__(self,projf):
        self.projf=projf   #projections file
        self.zip=False     #check if file is zipped
        self.tmpf=None     #uncompressed projf, if compressed projf
        self.ptf=None      #file pointer        
        self.crdf=None     #trajectory file
        self.nat=None      #number of atoms
        self.evecsf=None   #evecs file
        self.nvec=None     #number of modes
        self.nfr=None      #number of frames
        
        self.currfr=0      #current frame number
        self.fr=None       #current array of projections
        self.frs=None      #array of all projections

    def close(self):
        if self.ptf: self.ptf.close()
        if self.zip: os.system('/bin/rm '+self.tmpf)

    def __del__(self):
        self.close()

    def open(self):
        """open and prepare projections file"""
        if self.ptf: self.ptf.seek(0)
        else:
            for suffix,tool in {'.gz':'gunzip', '.bz2':'bunzip2'}.items():
                if suffix in self.projf:
                    self.zip=suffix
                    x,self.tmpf=mkstemp(dir=tmpdir)
                    os.system('%s -c %s > %s'%(tool,self.projf,self.tmpf))
                    self.ptf=open(self.tmpf)
                    break
            self.read_nvec()
            self.read_nfr()
        self.currfr=0
        
    def read_nvec(self):
        """read number of modes"""                 
        self.open(); l=self.ptf.readline()
        while '#' in l: l=self.ptf.readline()
        self.nvec=len(l.split())
        
    def read_nfr(self):
        """read number of frames"""
        self.open();
        cmd='tail -2 %s|grep -e "#"'%(self.tmpf)
        z=os.popen(cmd).readline()
        self.nfr=int(z.split()[-1])
        
    def populate(self,crdf,nat,evecsf):
        """calculate the projections"""
        self.crdf=crdf     #trajectory file
        self.nat=nat       #number of atoms
        crdobj=crd(self.crdf,nat)
        self.evecsf=evecsf #evecs file
        evecobj=evecs(evecsf)
        evecobj.projectCRD(crdobj,self.projf,nq=1000)

    def readFrame(self,skip=0):
        """read next frame.
        NOTE: not implemented skipping frames"""
        deltaf=1+skip
        l=self.ptf.readline()
        while l and '#' in l: l=self.ptf.readline()
        if not l: return 0
        self.fr=numpy.array( [float(x) for x in l.split()] )
        self.currfr+=deltaf
        return deltaf

    def loadProjections(self):
        """read all modes into memory"""
        self.open()
        self.prjs=numpy.zeros(self.nvec*self.nfr)
        self.prjs=self.prjs.reshape(self.nvec,self.nfr)
        ifr=0
        while self.readFrame():
            for ivec in range(self.nvec):
                self.prjs[ivec][ifr]=self.fr[ivec]
            ifr+=1#; print ifr
        return self.prjs

class mwcovar:

    """mwcovar is the output mass weighted covariant matrix from ptraj
    command:
                          C = <X^t M X>
    UNITS: Angstroms^2 * gram/mole """

    def __init__( self, inp ):
        
        init_dict={
            'file' : '',              #file storing the matrix
            'n' : 0,                  #dimension
            'mwc' : numpy.zeros( 1 ), #the covarian matrix
            'M': numpy.zeros( 1 ),    #vector of masses
            'sM': numpy.zeros( 1 ),   #sqrt(mass)
            'isM': numpy.zeros( 1 ),   #1/sqrt(mass)
            }

        self.__dict__.update( init_dict )

        #elucidate type of input
        if inp==None:
            pass
        if type( inp ) == type( '' ):
            self.file = inp
            self.load_mwc( inp )
        elif type( inp ) in ( type( [] ),type( () ) ):
            self.n = int( sqrt( len( inp) ) )
            self.mwc = numpy.array( inp ).reshape( self.n, self.n )
        elif type( inp ) == type( numpy.zeros( 1 ) ):
            self.n   = int( sqrt( inp.size ) )
            self.mwc = deepcopy( inp )
            self.mwc.reshape( self.n, self.n )

    def load_mwc( self, inp ):
        """load matrix onto memory"""
        x = []
        for line in open( inp ).readlines():
            x += [ float(z) for z in line.split() ]
        n = int( sqrt( len(x) ) )
        self.mwc = numpy.array( x ).reshape( n, n )
        self.n = n
        return self.mwc
    
    def normalize( self ):
        """divide each entry by the uncorrelated correlation"""
        mwc = numpy.zeros( self.n * self.n ).reshape( self.n, self.n )
        for i in range(self.n):
            a=sqrt( self.mwc[i][i] )
            for j in range(self.n):
                b = sqrt( self.mwc[j][j] )
                mwc[i][j] = self.mwc[i][j]/(a*b)
        self.mwc = mwc
        return self.mwc

    def findModes(self):
        """NOTE: this modes are normalized to 1 (q^t * q = 1) as
        opposed to the class evecs normalization (q^t * M * q = 1)
        
        Find normal modes given a variance-covariance matrix
        varcovar=E( (x-mu)^t (x-mu) )    (NOTE: ^t means transpose)
        evals,evecs are nd, (nd x nd) matrices, respectively"""
        from numpy.linalg import eigh
        #varcov = evecs * evals * evecs^+
        [self.evals,self.evecs]=eigh( self.mwc ) #diagonalize varcovar
        """
        order from higher to smaller correlation. This corresponds to
        the introduction of a 'permutation matrix' P, ie:
        C =  E * L * E^+ = (E * P^+) * ( P * L * P^+) * (P * E^+),
        where E^+ represents the hermitian conjugate of E
         P * E^+ = (E * P^+)^+ permutes the rows of E^+
         E * P^+ permutes the columns of E, ie, permutes the eigenvectors
         P * L * P^+ permutes columns and rows of L to reorder the eigenvalues

        This reordering leaves unaltered the 'frame of reference' of
        the coordinates from which varcovar was generated, that is,
        if we think of the quadratic potential 
         V  =  x^+ * C^+ * x  =  x^+ * (E * L * E^+)^+ * x  =  x^+ * (E * L^+ * E^+) * x, we know have
         V  =  x^+ * ( (E * P^+) * ( P * L * P^+) * (E * P^+)^+ )^+ * x = 
            =  x^+ * ( (E * P^+) * ( P * L^+ * P^+) * (E * P^+)^+ ) * x, and we see that coordinates 'x'
        do not neet be multiplied by any matrix.
        """
        perm=numpy.argsort(self.evals)[::-1]
        self.evals=self.evals[perm]
        self.evecs=self.evecs[:,perm]
        return { 'evals':self.evals, 'evecs':self.evecs }

    def massProcess(self,mlist):
        """obtain several mass vectors from the list of atomic
        masses. It is suggested you pass the mass list from a 'top'
        object."""
        self.M=[]
        for m in mlist: self.M+=[m,m,m]
        self.M=numpy.array(self.M)
        self.sM=numpy.sqrt(self.M)
        self.isM=1.0/self.sM
        return self.M

    def pseudoinverse(self,co=None,tZero=True):
        """find the pseudoinverse
        co: set eigenvalue exactly to zero if absolute value smaller than co
        Tzero: set the last 6 modes, assumed to be the translation ones, to zero
        """
        self.findModes()
        K = self.evals
        if tZero: K[-6:] = 0.0  #set translation modes exactly to zero
        if co: K[ numpy.abs(K) < co ] = 0.0 #zero out elements smaller than co
        nonzeroelements = numpy.where(numpy.abs(K)>0)
        K[ nonzeroelements ] = 1 / K[ nonzeroelements ] #invert non-zero elements
        return K, self.evecs
        
    def force_constant_matrix( self, T=300, mlist=[],
                               co_map=numpy.zeros(1), model=None ):
        """find the force constant matrix, according to:
        K = beta * Mass * mwcovar^(-1)  [Kcal/(mol*Angstrom^2)]
        co_map : cut-off map between degrees of freedom"""
        dot = numpy.dot
        if not self.M.any(): self.massProcess(mlist)
        K, D = self.pseudoinverse( )
        K = numpy.diag( K )
        K = numpy.dot( numpy.diag( self.M ), K )
        K = dot( D.transpose(), dot( K, D ) )
        if co_map.any():
            K = K * co_map
        K = K * (T/300.) * 1.99E-03
        return K
    
    def saveModes(self, outf, N=None, avg=None, atomic_masses=[],
                  eigFormat = 'frequency', format = 'ptraj evecs',
                  T = 300 ):
        """save eigenvectors of the matrix"""
        from utilities.readingWritingFiles import write_from_numpy
        if 'evecs' not in self.__dict__.keys( ): self.findModes( )
        if not N: N = self.n
        buf = ''
        #determine the output format for the eigenvalues
        eig = list( self.evals )[:N]
        if eigFormat == 'frequency':
            """output eigenvalues as frequencies"""
            for i in range(N):
                """conversion from mass-weighted fluctuation, in
                Angstrom^2 * (gram/mole), to wavenumbers, in 1/cm:
                equipartition theorem:
                      kT/2 = 1/2 m w^2 dx^2 = 1/2 w^2 * eval
                      w = sqrt( kT / eval )
                      kT in Kcal/mol
                      eval in Anstrom^2 * (gram/mole)
                from frequency to wavenumber in :
                      w*L = c (light speed)
                      wavenumber = 1 / (2*pi*L)
                """
                eig[i] = 108.587 * sqrt( (0.6*T/300) / abs(self.evals[i]) )
                if self.evals[ i ] < 0: eig[ i ] = -eig[ i ]
        #determine the output format
        if format == 'ptraj evecs':
            #ptraj evecs have the normalization q^t * M * q = 1
            if not atomic_masses:
                print 'ERROR: need atomic masses for this format'
                return buf
            isM = []
            for i in range( len( atomic_masses ) ):
                for j in range( 3 ): isM.append(1/sqrt(atomic_masses[ i ]))
            isM = numpy.array( isM )
            #write header
            buf += ' Eigenvector file: MWCOVAR\n'
            buf += '%5d %5d\n'%( self.n, self.n )
            #write average snapshot
            if avg == None:
                avg = numpy.zeros( self.n )
            else:
                avg = numpy.array( avg ).flatten( )
            buf += write_from_numpy( None, avg, format='%11.5f', ncols=7 )
            buf += '****\n'
            #write eigenvectors
            for i in range( N ):
                buf += '%5d %11.5f\n'%( 1+i, eig[ i ] )
                v = isM * self.evecs[ i ]
                buf += write_from_numpy(None, v, format='%11.5f', ncols=7)
                if i < N-1: buf += '****\n'
        if outf: open(outf,'w').write(buf)
        return buf
                
    def sav2fil( self, out, format=None ):
        """save matrix to file in a variety of formats
         format=
         None: as a square matrix
         as 3 columns: suitable for sigmaplot visualization"""
        from utilities.readingWritingFiles import write_from_numpy
        pt = open( out, 'w' )
        if format==None:
            write_from_numpy (pt, self.mwc, format=' %6.3f', ncols=self.n)
        elif format=='full precision':
            write_from_numpy (pt, self.mwc, format=' %14.8f', ncols=self.n)
        elif format=='as 3 columns':
            buf=''
            for i in range(self.n):
                for j in range(self.n):
                    buf+='%3d %3d %6.3f\n'%( i+1, j+1, self.mwc[ i ][ j ] )
            pt.write( buf )
        elif format=='isotropic average':
            """average over X,Y,Z directions"""
            m = self.n / 3
            buf = numpy.zeros( m * m ).reshape( m, m )
            for i in range(0, self.n, 3 ):
                for j in range(0, self.n, 3 ):
                    x = 0.0
                    for k in range( 3 ):
                        x += self.mwc[ i+k ][ j+k ]
                    if i==j: x=0.00 #ignore self-correlation
                    buf[i/3][j/3] = x/3
            write_from_numpy (pt, buf, format=' %6.3f', ncols=m)
        elif format=='as 3 columns, isotropic average':
            """average over X,Y,Z directions"""
            buf = ''
            for i in range(0, self.n, 3 ):
                for j in range(0, self.n, 3 ):
                    x = 0.0
                    for k in range( 3 ):
                        x += self.mwc[ i+k ][ j+k ]
                    if i==j: x=0.00 #ignore self-correlation
                    buf += '%3d %3d %6.3f\n'%( 1+(i/3), 1+(j/3), x/3.0 )
            pt.write( buf )
        else:
            print 'ERROR: format %s not found!'%format

    def total_fluct( self ):
        """output the trace of the matrix, which is nothing but the
        sum of the eigenvectors, that is, the sum of the particle
        fluctuations"""
        return self.mwc.trace()        

class secstruct:

    """class to handle ouput of sescstruct output from the secstruct
    command of ptraj.
    Attributes of class secstruct:
    N      : number of residues
    resNum : numpy.array[N] residue numbers
    profile: numpy.array[N][7]. The row index corresponds to the
             percentages of the different secondary structure types,
             which are (1)3-10-Helix, (2)alpha-helix, (3)PI-Helix,
             (4)parallel-Sheet, (5)antip.-Sheet, (6)Turn, and (7)Coil
    """

    def __init__( self ):

        pass

    def load( self, datFile=None, sumFile=None ):
        if sumFile:
            """Structure of the sumFile
            #ResNum 3-10-Helix alpha-helix PI-Helix parallel-Sheet antip.-Sheet Turn
            """
            Y = []; R = []
            for line in open(sumFile).readlines():
                if line[0] == '#': continue
                items = line.split()
                R.append( int( items[0] ) ) #residue number
                y = [ float(x)/100 for x in items[1:] ]

                y.append( 1.0 - sum(y) ) #coil
                Y.append( y )
            self.resNum  = numpy.array( R )                
            self.profile = numpy.array( Y )
            self.N = len( Y )
            return self.profile

    def varS( self, sel=numpy.zeros(1) ):
        """compute the following variance matrix from the profile
        sel: numpy list of indexes, if we want a subset"""

        if not sel.any(): sel = numpy.arange( self.N ) #all residues
        N = len(sel)
        P = self.profile[sel] #profile of the selection
        
        V = numpy.zeros( N * N ).reshape(N,N)
        X2 = P * P
        XR = numpy.sqrt( X2.sum( axis=1 ) )
        for i in range(N -1):
            for j in range(i+1,N):
                d = X[i]-X[j]
                V[i][j] = (d*d).sum()/(XR[i]*XR[j])
                V[j][i] = V[i][j]
        return V
    
    def varX( self, s2, sel1=numpy.zeros(1), sel2=numpy.zeros(1) ):
        """same as varS but now between two different profiles
        sel1: numpy list of indexes in self, if we want a subset
        self2: same as sel1 but for s2
        """

        if not sel1.any(): sel1 = numpy.arange( self.N ) #all residues
        N = len(sel1)
        P1 = self.profile[sel1] #profile of the selection
        
        if not sel2.any(): sel2 = numpy.arange( s2.N ) #all residues
        N2 = len(sel2)
        P2 = s2.profile[sel2] #profile of the selection

        V = numpy.zeros( N * N2 ).reshape(N,N2)
        X = P1 * P1
        XR = numpy.sqrt( X.sum( axis=1 ) )
        X2 = P2 * P2
        XR2 = numpy.sqrt( X2.sum( axis=1 ) )

        for i in range(N):
            for j in range(N2):
                d = X[i]-X2[j]
                V[i][j] = (d*d).sum()/(XR[i]*XR2[j])
        return V
