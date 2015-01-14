#############################
#  Run in text mode
#  vmd -dispdev text -eofexit -e script.tcl
#############################

set pdbfile "8styrene32.pdb"
set indcd   "equil.in.dcd"
set outdcd  "equil.out.dcd"

# The PDB file adds one unwanted conformation.
set begframe 1
# -1 signal last frame
set endframe -1

mol new $pdbfile
animate read dcd $indcd beg 0 end -1 waitfor all
animate write dcd $outdcdfile beg $begframe end $endframe waitfor all
quit

