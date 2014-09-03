
# 2007-12-12 
# written by Krzysztof Murzyn, krzymu@gmail.com
from MMTK import *
from Scientific.IO import NetCDF
import re, time, Numeric, string

# ac/protA:1_248;ac/POPC:OLE,GLC,PCH,PAL;m/water:SOL
# ac/protA:1_248;ac/protA;ac/POPC:OLE,GLC,PCH,PAL;m/water:SOL

TRUE = 1
FALSE = 0

def parseMolSpec(pdb,specstr):
    mapping = {}
    chunks = [x.split(':') for x in specstr.split(';')]
    order = []
    for name, spec in chunks:
        if mapping.has_key(name) and not spec: # ambiguity in names
            raise ' %s name is not unique' % name
        if not spec: # name already known
            if not mapping.has_key(name):
                raise ' %s name lacks definition' % name
            order.append(name)
            # check this: molA, molA, molB, molA, molC, molC
        elif '_' in spec: # numeric range of residues, e.g. in peptide chains
            first, last = map(int,spec.split('_'))
            mapping[name] = [pdb[i-1].name for i in range(first,last+1)]
            order.append(name)
        else: # list of residue names or single residue name
            mapping[name] = [x.strip() for x in spec.split(',')]
            order.append(name)
    uord = [order[0]] # I hope this works :)
    for i in range(1,len(order)):
        if uord[-1] != order[i]: uord.append(order[i])
    return mapping, uord

def getMolecularComposition(pdb,mapping,order,objects=None,verbose=0):
    """ mapping['mmtk_molname'] = ['P21','P22','P23',...,'PHG',...]
    """
    
    molecules = {}
    for m in mapping.keys(): molecules[m] = []

    i = 0
    added = 0
    iname = 0
    while iname < len(order):
        name = order[iname]
        molsize = len(mapping[name])
        pdbnames = [x.name for x in pdb[i:i+molsize]]
        mismatch = pdbnames != mapping[name]
        if not added:
            if mismatch:
                raise ' Residue names mismatch:\n [pdb]: %s\n [spec] %s' % \
                      ('-'.join(pdbnames), '-'.join(mapping[name]))
        else:
            if mismatch:
                added = 0
                iname += 1
                continue
        molecules[name].append(pdb[i:i+molsize])
        i += molsize
        added = 1
    if objects is None:
        for m in molecules.keys(): molecules[m] = len(molecules[m])
    # else return objects and not number of objects
    return molecules

def __guessElementName(pdbatom):

    atomname = pdbatom.name
    while atomname[0] in string.digits:
        atomname = atomname[1:] + atomname[0]
    try:
        element = pdbatom['element']
        a = ChemicalObjects.Atom(element, name = atomname)
    except KeyError:
        try:
            element = atomname[:2]
            a = ChemicalObjects.Atom(element, name = atomname)
            
        except IOError:
            element = atomname[:1]
            a = ChemicalObjects.Atom(element, name = atomname)
    return element

def pdb2mmtk(pdb,mapping,order=None):

    collection = Collection()
    mc = getMolecularComposition(pdb,mapping,order,objects=TRUE)
    if order is None: order = mc.keys()
    if not [x for x in order if x in mc.keys()]:
        raise ' Error: Invalid order (%s vs %s)' % (','.join(order),
                                                    ','.join(mc.keys()))
    for name in order:
        for pdb_cluster in mc[name]:
            otype, oname = name.split('/')
            if otype.upper() == 'AC': # AtomCluster
                pdbmap, i, atoms = [], 0, []
                for residue in pdb_cluster:
                    pdbdict = {}
                    for atom in residue:
                        aname = atom.name
                        element = __guessElementName(atom)
                        a = ChemicalObjects.Atom(element,name=aname)
                        a.setPosition(atom.position*Units.Ang)
                        anum = atom.properties['serial_number']
                        a.setIndex(anum - 1)
                        atoms.append(a)
                        pdbdict[atom.name] = Database.AtomReference(i)
                        i = i + 1
                    pdbmap.append((residue.name, pdbdict))
                m = ChemicalObjects.AtomCluster(atoms, name = oname)
                #m.description dc=[(0,1,.124),...]
                # bewere that in MMTK2.4.2 neither dc nor pdbmap is stored
                # in universe description
                m.pdbmap = pdbmap
            elif otype.upper() == 'M': # Molecule
                m = Molecule(oname)
                ires = 0
                resnames = [x[0] for x in m.pdbmap] + mapping[name]
                atom_list = m.atomList()
                for residue in pdb_cluster:
                    if residue.name not in resnames:
                        raise ' Error: residue name mismatch (%s vs %s)' % \
                              (residue.name,','.join(resnames))
                    for atom in residue:
                        if hasattr(m,'pdb_alternative') and \
                           m.pdb_alternative.has_key(atom.name):
                            aname = m.pdb_alternative[atom.name]
                        else:
                            aname = atom.name
                        mmtk_atom = atom_list[m.pdbmap[ires][1][aname].number]
                        anum = atom.properties['serial_number']
                        mmtk_atom.setIndex(anum - 1)
                    ires += 1
            else:
                raise ' Error: Unrecognised object specification (%s)' % item
            collection.addObject(m)
    return collection


###

def readAmberTop(topfn):
    f = open(topfn,'r')
    header = f.readline()
    if not re.match('^%VERSION',header):
        raise 'Invalid AmberTopology file'
    dtop = f.read()
    f.close()
    #
    split_pattern = re.compile('^%FLAG',re.M)
    chunk_pattern = re.compile(' (\S+)\s+^%FORMAT\((\S+)\)\s+$(.+)',re.M+re.S)
    chunks = [re.match(chunk_pattern,x).groups() \
              for x in re.split(split_pattern,dtop) if x]
    result = {}
    format_pattern = re.compile('\d*([aeif])(\d+)\.{0,1}\d*',re.I)
    field_types = {'f': float, 'i': int}
    for name,format,data in chunks:
        if result.has_key(name):
            raise 'Duplicated key found:',name
        field_type, field_length = re.match(format_pattern,format).groups()
        field_pattern = re.compile('.{%s}' % field_length)
        if field_type.lower() in field_types.keys():
            try:
                result[name] = map(field_types[field_type.lower()],
                                   re.findall(field_pattern,data))
            except:
                print name, format
                return data
        else: result[name] = re.findall(field_pattern,data)
    return result

def readCRDFrame(nlines,filename):
    """read first frame of a crd file"""
    import gzip
    nc=None
    if filename[-3:]=='.gz': nc=gzip.open(filename)
    else: nc=open(filename)
    xyz=[]
    for iline in range(nlines): xyz+=nc.readline().split()
    return {'COORDINATES':xyz,'BOX':[0.0,0.0,0.0],'TIME':1}
    
def readNetCDFFrame(nc,frame=0):
    if Numeric.nonzero(Numeric.not_equal(\
        Numeric.ravel(nc.variables['cell_angles'][:]),90.)):
        raise 'Only rectangular box is implemented'
    return {'COORDINATES': nc.variables['coordinates'][frame],
            'TIME': nc.variables['time'][frame],
            'BOX': nc.variables['cell_lengths'][frame]}

def dumpPDB(pdbfn,topology,frame):
    dump = ['REMARK    GENERATED BY CNVA9',
            'HEADER    %s' % '%s :: t=%-14.5f' % \
            (''.join(topology['TITLE'])[:40], frame['TIME'],),
            'CRYST1 %8.3f%8.3f%8.3f  90.00  90.00  90.00 P 1' % \
            tuple(frame['BOX'])]
    atom_line = 'ATOM  %5d  %-3s %-3s %5d    %8.3f%8.3f%8.3f'
    natom = topology['POINTERS'][0]
    nres = topology['POINTERS'][11]
    resnumb = 0
    first = 1
    for i in range(1,natom+1):
        index = i-1
        x,y,z = frame['COORDINATES'][index]
        if i == first:
            resname = topology['RESIDUE_LABEL'][resnumb].strip()[:3]
            resnumb += 1
            if resnumb == nres: first = natom+1
            else: first = topology['RESIDUE_POINTER'][resnumb]
        dump.append(atom_line % (i,topology['ATOM_NAME'][index].strip()[:3],
                                 resname,resnumb,x,y,z))
    dump.append('')
    f = open(pdbfn,'w')
    f.write('\n'.join(dump))
    f.close()

if __name__ == '__main__':
    t = time.time()
    topology = readAmberTop('/projects/rubredoxin/out/1vcx_Fe/T_100/Prod/init/1vcx_Fe.uns.top')
    print ' Read in %.3f' % (time.time()-t,)
    print topology.keys()
    #nc = NetCDF.NetCDFFile('mdcrd','r')
    #frame = readNetCDFFrame(nc)
    nlines=int(natoms*3/10)
    if natoms%10: nlines+=1
    frame=readCRDFrame(nlines,'/projects/rubredoxin/out/1vcx_Fe/T_100/Prod/prev/1vcx_Fe.00001_25000.crd.gz')
    nc.close()
    dumpPDB('test.pdb',topology,frame)
    print 'Bye!' ; raise SystemExit,1
    spec = 'ac/protA:1_248;ac/PC:LIP;m/water:WAT'
    pdb = PDB.Structure('test.pdb')
    mapping,order = parseMolSpec(pdb,spec)
    mc = getMolecularComposition(pdb,mapping,order)
    collection = cnva9.pdb2mmtk(pdb,mapping,order=order)
    universe = Universe.OrthorhombicPeriodicUniverse((0., 0., 0.))
    universe.addObject(collection)
    universe.initializeVelocitiesToTemperature(0.)

