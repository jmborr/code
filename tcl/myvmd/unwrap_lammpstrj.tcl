mol new _PSF_FILE_
animate read lammpstrj _LAMMPSTRJ_FILE_ beg 0 end -1 waitfor all
pbc wrap -all
pbc unwrap -all
#This commands only for visualizating purposes
#pbc box_draw
#pbc join residue -all
animate write dcd _DCDFILE_ beg 1 end -1 waitfor all
quit
