#############################
#  Run in text mode
#  vmd -dispdev text -eofexit -e script.tcl
#############################

set pdbfile "8styrene32.pdb"
set indcdfile1 "equil.1.dcd"
set indcdfile2 "equil.2.dcd"
set indcdfile3 "equil.3.dcd"
set outdcdfile "equil.dcd"

# The PDB file adds one unwanted conformation
set begframe 1
# Filter out the first N conformations
set filtered  0
set begframe [expr {$begframe+$filtered}]
# -1 signal last frame
set endframe -1

mol new $pdbfile
animate read dcd $indcdfile1 beg 0 end -1 waitfor all
animate read dcd $indcdfile2 beg 0 end -1 waitfor all
animate read dcd $indcdfile3 beg 0 end -1 waitfor all

animate write dcd $outdcdfile beg $begframe end $endframe waitfor all

quit

