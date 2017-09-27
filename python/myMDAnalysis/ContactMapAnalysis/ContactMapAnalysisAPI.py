from __future__ import print_function

from tqdm import tqdm
import numpy as np
import sys
from tqdm import tqdm
from scipy.sparse import csr_matrix
import MDAnalysis as MDA
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch

import myMDAnalysis.ContactMapAnalysis.ContactMapAnalysis as CMA

from pdb import set_trace as trace

"""
Desired properties:
(1) parser defining how to evaluate a contact
(2) iterate over the trajectory frames
(3) compute the contact map with the given parser and frame, returning a crs_ContactMap object
"""

def AtomsListPerResidue(atoms, shift=True):
    """Returns a list whose elements are also lists. Each element contains the indices                                  
    of the atoms contained in 'atoms' that belong to a particular residue.                                              
    Note: not all atoms belonging to a particular residue may be contained in 'atoms'                                   
    shift=True, by default, numbering of atoms and residues is lost.                                                    
          Example: if 'atoms' contains atoms with numbers 3,5,8,10,12, the numbers become 0,1,2,3,4.                    
                   Similary, if 'atoms' contains atoms in residues 5,9,12, the number become 0,1,2                      
    """
    ResDict={} #as many entries as residues in atoms                                                               
    for resNum in atoms.resnums():
        ResDict[resNum]=[]
    iat=0 #we loose the atom number convention in the PSF file                                                          
    for atom in atoms:
        ResDict[ atom.resid ].append(iat)
        iat+=1
    ResList=[]
    for key in sorted(ResDict.iterkeys()): ResList.append(ResDict[key])
    return ResList    


class ContactMapProtocol:
    """
    Rules to generate a contact map between two sets of atoms
    """

    def __init__(self, row=None, col=None, **kwargs):
        """
        Instantiate the various attributes
        :param row: atom selection string for the row index of the Contact map
        :param col: atom selection string for the column index of the Contact map
        :param kwargs: additional arguments to update self.__dict__
        """
        self.row_selection = row
        self.col_selection = col
        self.cutoff = None  # cut off defining the contact, in Angstroms
        self.row_group = None  # atom group object, from applying selectAtoms to row_selection
        self.col_group = None  # object, from applying selectAtoms to col_selection
        self.shift = None  # set to 1 if the first atomic index in Universe is zero
        self.isInitialized=False  # initialized() hast not been called
        self.byres=False  #if set to True, contact map between residues instead of between atoms
        # reducer is an algorithm to define the contact between two residues based
        # on the atomic contact map between atoms of the two residues
        # Example: if residue is made up of atoms 0, 1, 2 and residue j made up of atoms 4, 5
        # then C[[0,1,2]][:,[4,5]] is the desired C_ij submatrix.
        self.reducer=np.any  #  np.any( C_ij ) is evaluated by default
        # filter is a function to filter non-wanted contacts (e.g self-contacts).
        # Takes a csr_matrix as argument
        self.filter = None
        self.n_residues = None
        # np array where index is atom number and item is resid attribute
        # of MDAnalysis atom object.The first populated index is self.shift
        self.atom2res = list()
        """
        #One example of filtering contacts
        def filter(csr):
            # filter out contacts that are too close along the sequence
            # csr: contact map, a csr_matrix object
            cutoff=1
            i=0
            for irow in range( len(csr.indptr) - 1 ):
                for icol in csr.indices[ csr.indptr[irow]:csr.indptr[irow+1]]:
                    if abs(irow-icol)<=cutoff: csr.data[i] = 0
                    i += 1
            csr.eliminate_zeros()
            return csr
        """
        self.__dict__.update(kwargs)

    def initialize(self, a_universe):
        """
        Populate various attributes starting from the row and column selection strings
        :param a_universe: MDAnalysis.Universe object
        """
        if self.isInitialized: return None
        try:
            self.row_group = a_universe.selectAtoms(self.row_selection)
            self.col_group = a_universe.selectAtoms(self.col_selection)
            self.shift=0
            if a_universe.atoms[0].number==0:
                self.shift=1
        except AttributeError:
            print('row_selection and/or col_selection undefined in protocol\n')
        if self.byres:
            """generate atom number to residue number conversion array"""
            self.n_residues = a_universe.residues.numberOfResidues()
            tmp={}
            for atom in a_universe.atoms:
                tmp[atom.number] = atom.resid
            N = max(tmp.keys()) + self.shift
            self.atom2res = np.zeros(N, dtype=int)
            for (atom_number, residue_id) in tmp.items():
                self.atom2res[atom_number] = residue_id


def GenerateContactMap(time_step, protocol):
    """
    Read a frame and compute the contactMap between two groups of atoms.
    :param time_step: MDAnalysis frame object
    :param protocol: ContactMapProtocol object. if protocol.byres=True, then
    report the number of atomic contacts for a given residue pair in contact
    Numbering of atoms and residues is the one given by the PSF file
    :return csr_ContactMap instance
    """
    row_group_tree = AtomNeighborSearch(protocol.row_group) #KDE-tree object for row_group
    col_group_tree = AtomNeighborSearch(protocol.col_group) #KDE-tree object for col_group
    # atoms of row_group in contact with col_group
    row_close = row_group_tree.search_list(protocol.col_group, protocol.cutoff)
    if not row_close: #no contacts
        shape=(time_step.numatoms, time_step.numatoms) #all atoms in the system
        if protocol.byres:
            shape=(protocol.n_residues,protocol.n_residues)
        csr=csr_matrix( ([],([],[])), shape, dtype='int32')  #empty map
        cma=CMA.csr_ContactMap(csr)
        cma.cut_off = protocol.cutoff
        return cma
    # atoms of col_group in contact with row_close/row_group
    col_close = col_group_tree.search_list(row_close,protocol.cutoff)
    dd = distance_array(row_close.coordinates(), col_close.coordinates(), time_step.dimensions[0:3])
    (row_indices, col_indices) = np.where(dd < protocol.cutoff)
    row_indices = row_close.indices()[row_indices]
    col_indices = col_close.indices()[col_indices]
    # switch from atomic indices to residue numbers if needed
    if protocol.byres:
        row_indices = protocol.atom2res[row_indices]
        col_indices = protocol.atom2res[col_indices]
        shape = (1 + row_indices.max(), 1 + col_indices.max())
    # Take into account if the first atomic index was zero
    else:
        row_indices += protocol.shift
        col_indices += protocol.shift
        shape = (time_step.numatoms, time_step.numatoms)  # all atoms in the system
    #create sparse matrix for the contact map. 
    data=np.ones(len(row_indices), dtype='int32')  # just signaling a contact map
    csr = csr_matrix((data, (row_indices, col_indices)), shape, dtype='int32')
    if protocol.byres:
        csr.data[:] = 1  #overwrite number of atomic contacts per residue pair with 1 only
    if protocol.filter:
        csr = protocol.filter(csr)  # Filtering
    cma = CMA.csr_ContactMap(csr)
    cma.cut_off = protocol.cutoff
    return cma
    
def GenerateContactMapList(PSF_file, traj_file, protocol, start_frame=0, end_frame=None, step=1):
    """Read a trajectory and compute the contactMap between two groups of atoms.
    :param PSF_file: topology file
    :param traj_file: trajectory file
    :param protocol: ContactMapProtocol instance
    :param start_frame: first frame to find the contact map
    :param end_frame: last frame to find the contact map. End of trajectory if default value
    :param step: calculate map every step number of frames
    :return a ContactMapList instance 
    """
    u=MDA.Universe(PSF_file, traj_file)
    if end_frame is None:
        end_frame = len(u)
    protocol.initialize(u)
    cml=CMA.ContactMapList()
    for ts in tqdm(u.trajectory[start_frame:end_frame:step]): # iterate through all frames
        cm=GenerateContactMap(ts,protocol)
        cml.append(cm)
    return cml



def solvent_bound(a_universe, solute, solvent, cut_off, box=None,
                  start_frame=0, end_frame=None, step=1):
    """
    Generator to read a trajectory, and for certain frames find
    whether each atom of solvent is bound any atom of the solute.
    :param a_universe: MDAnalysis.Universe instance
    :param solute: string representing solute atoms 
    :param solvent: string representing solvent atoms
    :param cut_off: maximum distance for considering a solvent atom
      bound to any solute atom.
    :param start_frame: first frame of the trajectory to process
    :param end_frame: last frame of the trajectory to process
    :param step: process trajectory every step frames
    :return: generator yielding for each processed frame the group of
    solvent atoms bound to any solute atom.
    """
    solute_group = a_universe.selectAtoms(solute)
    solvent_group = a_universe.selectAtoms(solvent)
    for _ in a_universe.trajectory[start_frame:end_frame:step]:
        solute_tree = AtomNeighborSearch(solute_group, box=box) # KDE-tree
        yield solute_tree.search(solvent_group, cut_off, level='A')


def solvent_bound_flag(a_universe, solvent, solute, cut_off,
                       start_frame=0, end_frame=None, step=1):
    """
    Read a trajectory and for selected frames find whether each atom of
    solvent is bound any atom of the solute.
    :param a_universe: MDAnalysis.Universe instance
    :param solvent: string representing solvent atoms
    :param solute: string representing solute atoms 
    :param cut_off: maximum distance for considering a solvent atom
      bound to any solute atom.
    :param start_frame: first frame of the trajectory to process
    :param end_frame: frame after last frame of the trajectory to process
    :param step: process trajectory every step frames
    :return: Lists of lists. Each list for a particular solute atom.
      Contents are a sequence of either True or False values, indicating
      if atom was bound to solvent at a particular frame.
    """
    if end_frame is None:
        end_frame = len(a_universe.trajectory)
    solvent_group = a_universe.select_atoms(solvent)
    solute_group = a_universe.select_atoms(solute)
    flags = [[] for _ in range(len(solvent_group))]
    nframe = int((end_frame - start_frame)/step)
    for _ in tqdm(a_universe.trajectory[start_frame:end_frame:step],
                  total=nframe):
        # creating the tree is costly. Done once for each frame
#        tree = AtomNeighborSearch(solvent_group,
#                                  box=a_universe.dimensions).kdtree
        tree = AtomNeighborSearch(solute_group).kdtree
        # parallelizable
        for i, atom in enumerate(solvent_group):
            tree.search(atom.position, cut_off)
            indices = np.array(tree.get_indices())
            is_bound = bool(indices.any())
            flags[i].append(is_bound)
    return flags
