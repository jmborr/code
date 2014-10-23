mol new _pdbfile_
animate read crd _mdcrd_file_ beg 0 end -1 waitfor all
animate write dcd _dcd_file_ beg 1 end -1 waitfor all
quit
