from pymol import cmd
 
def find_buried_waters(sele='all', cutoff=-1, state=1, quiet=1, _self=cmd):
    '''
DESCRIPTION
 
    Finds and turns a selection of waters determined to be buried (no solvent
    accessibility) in the original selection.
 
ARGUMENTS
 
    sele = string: atom selection {default: all}
 
    cutoff = float: threshold on what one considers an "exposed"
    atom (in A**2) {default: surface_residue_cutoff}
 
    state = int: object state {default: 1}
    '''
    cutoff, state, quiet = float(cutoff), int(state), int(quiet)
 
    tmpObj=_self.get_unused_name("__tmp")
    _self.create(tmpObj, sele, state, 1, zoom=0)
 
    _self.set("dot_solvent", 1, tmpObj);
    _self.get_area(tmpObj, state=1, load_b=1)
 
    if cutoff < 0:
        cutoff = _self.get("surface_residue_cutoff")
    _self.remove(tmpObj + " and not solvent")
    _self.remove(tmpObj + " and b > %s" % cutoff)
 
    exposed = set()
    _self.iterate(tmpObj, "exposed.add((chain,resv))", space=locals())
 
    selName = _self.get_unused_name("buried")
    _self.select(selName, "(%s) in %s" % (sele, tmpObj))
 
    # clean up
    _self.delete(tmpObj)
 
    if not quiet:
        print ' Found %d buried water atoms' % (len(exposed))
 
    return sorted(exposed)
 
cmd.extend('find_buried_waters', find_buried_waters)
