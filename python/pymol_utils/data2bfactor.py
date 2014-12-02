"""
Copyright (c) 2003 Robert L. Campbell
Copyright (c) 2010 Thomas Holder

Original by RLC, with extensive modifications and improvements by TH.

Please read below for instructions

contains the functions 
   data2b_atom(mol='',data_file='')
   data2b_res(mol='',data_file='')
   data2q_atom(mol='',data_file='')
   data2q_res(mol='',data_file='')
"""

import sys, re

comment = re.compile('^\s*$|^\s*#')

def atom_data_extract(data_file):
    """
    Read the specified 'by-atom' data file and extract the data from it
    and store it in parallel dictionaries specifying the data
    and residue names (both with keys of chain and residue number and atom name).

    The data file can contain comment lines starting with \"#\" (on
    lines by themselves).
    
    These comment lines are ignored.
    """
    bdat = {}
    chain = ''

    data_lines = file(data_file, 'rU')

    for line in data_lines:
        # ignore comment lines (beginning with a '#') or blank lines
        if not comment.match(line):
            words = line.split()

            # check number of columns of data
            if len(words) == 5:
                chain = words[0]
                resi = words[1]
                resn = words[2]
                name = words[3]
                if chain == '-':
                    chain = ''
                data = float(words[4])
            elif len(words) == 4:
                resi = words[0]
                resn = words[1]
                name = words[2]
                data = float(words[3])
            else:
                sys.stderr.write("Error in reading data files -- check number of columns")
                sys.exit(1)

            bdat.setdefault(chain, {}).setdefault(resi, {})[name] = (data, resn)

    return bdat

def residue_data_extract(data_file):
    """
    Read the specified 'by-residue' data file and extract the data from it
    and store it in parallel dictionaries specifying the data
    and residue names (both with keys of chain and residue number).
    The data file can contain comment lines starting with \"#\" (on
    lines by themselves).
    These comment lines are ignored.
    """
    bdat = {}
    chain = ''

    data_lines = file(data_file, 'rU')

    for line in data_lines:
        # ignore comment lines (beginning with a '#') or blank lines
        if not comment.match(line):
            words = line.split()

            # check number of columns of data
            if len(words) == 4:
                chain = words[0]
                resi = words[1]
                resn = words[2]
                if chain == '-':
                    chain = ''
                data = float(words[3])
            elif len(words) == 3:
                resi = words[0]
                resn = words[1]
                data = float(words[2])
            elif len(words) == 2:
                resi = words[0]
                data = float(words[1])
                resn = ''
            else:
                sys.stderr.write("Error in reading data files -- check number of columns\n")
                sys.exit(1)

            bdat.setdefault(chain, {})[resi] = (data, resn)
    
    return bdat

###########################################################################################
# for testing purposes:
# if calling this as a program on its own, read the pdb_file name from
# the command line and run residue_data_extract on it. (does not require
# importing cmd from pymol

if __name__ == '__main__':
    pdb_file = sys.argv[1]
    b_dict = residue_data_extract(pdb_file)
    for chain in sorted(b_dict):
        for resi in sorted(b_dict[chain]):
            b, resn = b_dict[chain][resi]
            print "b-factors %s %s %s %s  new B='%s'" % (pdb_file, chain, resn, resi, b)
    sys.exit()


###########################################################################################
# PyMOL stuff

from pymol import cmd

def data2b_atom(mol='', data_file='', prop='b', quiet=0):
    """
DESCRIPTION

    Alters the B-factor data by atom.

USAGE

    data2b_atom <mol>, <data_file>

    where <mol> is the molecular object whose B-factor data you wish
    to modify and <data_file> is a file contain the data (one value
    for each atom) The format of <data_file> should be:

         chain resi resn name data
    or
         resi resn name data

    (i.e. "chain" is optional if all atoms are in one chain). 
    Lines beginning with '#' are ignored as comments.

SEE ALSO

    data2b_res, data2q_atom, data2q_res
    """

    b_dict = atom_data_extract(data_file)
    quiet = int(quiet) == 1
    
    def b_lookup(chain, resi, name, b):
        try:
            if chain in b_dict:
                b = b_dict[chain][resi][name][0]
            else:
                b = b_dict[''][resi][name][0]
            if not quiet: print '///%s/%s/%s new: %f' % (chain, resi, name, b)
        except KeyError:
            if not quiet: print '///%s/%s/%s keeping: %f' % (chain, resi, name, b)
        return b
    stored.b = b_lookup

    cmd.alter(mol, '%s=stored.b(chain, resi, name, %s)' % (prop, prop))
    cmd.rebuild()

def data2b_res(mol='', data_file='', prop='b', quiet=0):
    """
DESCRIPTION

    Alters the B-factor data by residue.

USAGE

    data2b_res <mol>, <data_file>

    where <mol> is the molecular object whose B-factor data you wish to modify
    and <data_file> is a file contain the data (one value for each residue)
    The format of <data_file> should be:

         chain resi resn data
    or
         resi resn data

    (i.e. "chain" is optional). Lines beginning with '#' are ignored as comments.

SEE ALSO

    data2b_atom, data2q_atom, data2q_res
    """

    b_dict = residue_data_extract(data_file)
    quiet = int(quiet) == 1

    def b_lookup(chain, resi, name, b):
        try:
            if chain in b_dict:
                b = b_dict[chain][resi][0]
            else:
                b = b_dict[''][resi][0]
            if not quiet: print '///%s/%s/%s new: %f' % (chain, resi, name, b)
        except KeyError:
            if not quiet: print '///%s/%s/%s keeping: %f' % (chain, resi, name, b)
        return b
    stored.b = b_lookup

    cmd.alter(mol, '%s=stored.b(chain, resi, name, %s)' % (prop, prop))
    cmd.rebuild()

def data2q_atom(mol='',data_file=''):
    """
DESCRIPTION

    Alters the occupancy data by atom.

USAGE

    See data2b_atom

SEE ALSO

    data2q_res, data2b_atom, data2b_res
    """
    data2b_atom(mol, data_file, prop='q')

def data2q_res(mol='',data_file=''):
    """
DESCRIPTION

    Alters the occupancy data by residue.

USAGE

    See data2b_res

SEE ALSO

    data2q_atom, data2b_atom, data2b_res
    """
    data2b_res(mol, data_file, prop='q')

cmd.extend('data2b_res',data2b_res)
cmd.extend('data2b_atom',data2b_atom)
cmd.extend('data2q_res',data2q_res)
cmd.extend('data2q_atom',data2q_atom)