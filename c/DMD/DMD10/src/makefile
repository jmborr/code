all: ../bin/preRelax.linux ../bin/txt2txt.linux ../bin/genePDB.linux ../bin/geneNative.linux ../bin/pdb2fyspec.linux ../bin/fy2prop.linux ../bin/movie2pdbs.linux ../bin/txt2pdb.linux ../bin/readMovHead.linux ../bin/movie_analysis.linux ../bin/movie2pdbs_ewin.linux

GCC = g++
LIBF = ../include/pseudoAA.h ../include/pseudoPep.h ../include/pseudoTXT.h
IFLAG = -I../include
OFLAG = -O


../bin/preRelax.linux: preRelax.o
	${GCC} ${OFLAG} -o ../bin/preRelax.linux preRelax.o
preRelax.o: preRelax.cpp ${LIBF}
	${GCC} -c ${OFLAG} ${IFLAG} preRelax.cpp


../bin/txt2txt.linux: txt2txt.o
	${GCC} ${OFLAG} -o ../bin/txt2txt.linux txt2txt.o
txt2txt.o: txt2txt.cpp ${LIBF}
	${GCC} -c ${OFLAG} ${IFLAG} txt2txt.cpp


../bin/genePDB.linux: genePDB.o
	${GCC} ${OFLAG} -o ../bin/genePDB.linux genePDB.o
genePDB.o: genePDB.cpp ${LIBF}
	${GCC} -c ${OFLAG} ${IFLAG} genePDB.cpp


../bin/geneNative.linux: geneNative.o
	${GCC} ${OFLAG} -o ../bin/geneNative.linux geneNative.o
geneNative.o: geneNative.cpp ${LIBF}
	${GCC} -c ${OFLAG} ${IFLAG} geneNative.cpp


../bin/pdb2fyspec.linux: pdb2fyspec.o
	${GCC} ${OFLAG} -o ../bin/pdb2fyspec.linux pdb2fyspec.o
pdb2fyspec.o: pdb2fyspec.cpp ${LIBF}
	${GCC} -c ${OFLAG} ${IFLAG} pdb2fyspec.cpp


../bin/fy2prop.linux: fy2prop.o
	${GCC} ${OFLAG} -o ../bin/fy2prop.linux fy2prop.o
fy2prop.o: fy2prop.cpp ../include/second.h
	${GCC} -c ${OFLAG} ${IFLAG} fy2prop.cpp


../bin/movie2pdbs.linux: movie2pdbs.o movie.o
	${GCC} ${OFLAG} -o ../bin/movie2pdbs.linux movie2pdbs.o movie.o
movie2pdbs.o: movie2pdbs.cpp ../include/movie.h ${LIBF}
	${GCC} -c ${OFLAG} ${IFLAG} movie2pdbs.cpp
movie.o: movie.cpp ../include/movie.h
	${GCC} -c ${OFLAG} ${IFLAG} movie.cpp

../bin/movie2pdbs_ewin.linux: movie2pdbs_ewin.o movie.o
	${GCC} ${OFLAG} -o ../bin/movie2pdbs_ewin.linux movie2pdbs_ewin.o movie.o
movie2pdbs_ewin.o: movie2pdbs_ewin.cpp ../include/movie.h ${LIBF}
	${GCC} -c ${OFLAG} ${IFLAG} movie2pdbs_ewin.cpp


../bin/txt2pdb.linux: txt2pdb.o
	${GCC} ${OFLAG} -o ../bin/txt2pdb.linux txt2pdb.o
txt2pdb.o: txt2pdb.cpp ${LIBF}
	${GCC} -c ${OFLAG} ${IFLAG} txt2pdb.cpp


../bin/readMovHead.linux: readMovHead.o
	${GCC} ${OFLAG} -o ../bin/readMovHead.linux readMovHead.o movie.o
readMovHead.o: readMovHead.cpp ../include/movie.h
	${GCC} -c ${OFLAG} ${IFLAG} readMovHead.cpp


../bin/movie_analysis.linux: movie_analysis.o
	${GCC} ${OFLAG} -o ../bin/movie_analysis.linux movie_analysis.o movie.o
movie_analysis.o: movie_analysis.cpp ../include/movie.h
	${GCC} -c ${OFLAG} ${IFLAG} movie_analysis.cpp

clean:
	rm -rf *.o

