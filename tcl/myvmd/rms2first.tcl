# Fittrajectory to selected frame by RMSD transformation
# Use as: vmd -dispdev win -eofexit -e rms2first.tcl

# Call the RMSD trajectory tool
::rmsdtt::rmsdtt

######### Here set all variables

set firstframe 5002

#  frame to align the other frames to
set refframe $firstframe

#  write alined trajectory starting at this frame 
set begframe $firstframe

#  input and output files
set pdbfile "8styrene32.pdb"
set indcdfile "equil.dcd"
set outdcdfile "junk.dcd"

#  select all atoms
$::rmsdtt::w.top.left.sel delete 1.0 end
$::rmsdtt::w.top.left.sel insert end "all"

# Set modifiers
set ::rmsdtt::bb_only 0
set ::rmsdtt::trace_only 0
set ::rmsdtt::noh 0

# Set reference
set ::rmsdtt::rmsd_base top
set ::rmsdtt::traj_ref $firstframe
set ::rmsdtt::traj_sw 1

######### Here perform all transformations

# Load the trajectory
mol new $pdbfile
animate read dcd $indcdfile beg 0 end -1 waitfor all

# Align
::rmsdtt::doAlign

# Write trajectory
animate write dcd $outdcdfile beg $begframe end -1 waitfor all

#exit script
quit
