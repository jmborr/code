CC = g++
CFLAGS = -Wno-deprecated -c
CODE = /users/seldon/jmborr/Code
VPATH = $(CODE)/PDBClasses
INCLUDES = -I$(CODE)/PDBClasses
OBJS = amino_acid_param.o atom_param.o salt.o miscellanea.o \
       bibliography.o salt.o

#default target
all: pdbClasses $(OBJS) libheavy_atom.a message

#auxiliary targets
message:
	echo BUILDING HEAVY ATOM MODEL

pdbClasses: 
	cd $(CODE)/PDBClasses ; make all ;

amino_acid_param.o: amino_acid_param.cpp amino_acid_param.h miscellanea.h
	$(CC) $(CFLAGS) $<

atom_param.o: atom_param.cpp atom_param.h miscellanea.h bibliography.h
	$(CC) $(CFLAGS) $<

salt.o: salt.cpp salt.h
	$(CC) $(CFLAGS) $(INCLUDES) $<

miscellanea.o: miscellanea.cpp miscellanea.h bibliography.h pdbClasses2.h
	$(CC) $(CFLAGS) $(INCLUDES) $<

bibliography.o: bibliography.cpp bibliography.h
	$(CC) $(CFLAGS) $<

libheavy_atom.a:
	'rm' -f libheavy_atom.a;
	ar -r $@ $(OBJS)

clean:
	rm -f *.o core
	cd $(CODE)/PDBClasses ; make clean ;
