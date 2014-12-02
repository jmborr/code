#############################
#  Run in text mode
#  vmd -dispdev text -eofexit -e script.tcl
#############################

set pdbfile "8styrene32.pdb"
set indcdfile1 "equil.1.dcd"
set indcdfile2 "equil.2.dcd"
set indcdfile3 "equil.3.dcd"
set outdcdfile "equil.dcd"
# It could be each loaded file will add one unwanted conformation at the beginning of the trajectory
# Plus, the PDB file will also add and unwanted conformation. Hence the 5004 below, instead of 5000
set begframe 5004
# -1 signal last frame
set endframe -1

mol new $pdbfile
animate read dcd $indcdfile1 beg 0 end -1 waitfor all
animate read dcd $indcdfile2 beg 0 end -1 waitfor all
animate read dcd $indcdfile3 beg 0 end -1 waitfor all

animate write dcd $outdcdfile beg $begframe end $endframe waitfor all

quit

