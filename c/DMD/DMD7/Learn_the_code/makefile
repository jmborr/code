OFLAG = -g

dmd.out: bcp.o make_system.o search.o controls.o cluster.o corr_func.o movie.o bonds.o rms.o fs_output.o schedule.o 
	gcc -o dmd.linux ${OFLAG} -lm bcp.o make_system.o search.o controls.o cluster.o movie.o bonds.o rms.o corr_func.o fs_output.o schedule.o

bcp.o: bcp.c bcp.h make_system.h controls.h cluster.h movie.h search.h bonds.h rms.h corr_func.h fs_output.h schedule.h
	gcc -c ${OFLAG} bcp.c

make_system.o: make_system.c bcp.h make_system.h controls.h cluster.h search.h movie.h bonds.h rms.h corr_func.h  
	gcc -c ${OFLAG} make_system.c

search.o: bcp.h controls.h search.h search.c 
	gcc -c ${OFLAG} search.c

controls.o: controls.c bcp.h controls.h cluster.h fs_output.h
	gcc -c ${OFLAG} controls.c

cluster.o: cluster.c bcp.h bonds.h cluster.h controls.h rms.h
	gcc -c ${OFLAG} cluster.c

movie.o: movie.c bcp.h controls.h rms.h bonds.h make_system.h
	gcc -c ${OFLAG} movie.c

bonds.o: bonds.c bonds.h 
	gcc -c ${OFLAG} bonds.c

rms.o: rms.c rms.h bcp.h controls.h search.h
	gcc -c ${OFLAG} rms.c

corr_func.o: corr_func.c bcp.h corr_func.h controls.h bonds.h 
	gcc -c ${OFLAG} corr_func.c

fs_output.o: fs_output.c fs_output.h bcp.h controls.h make_system.h
	gcc -c ${OFLAG} fs_output.c

schedule.o: schedule.c schedule.h bcp.h controls.h
	gcc -c ${OFLAG} schedule.c

clean:
	rm -rf *.o
