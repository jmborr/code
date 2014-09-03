CC = g++
CFLAGS = -Wno-deprecated -c

#Default target
all: pdbClasses2.o #message

message:
	echo BUILDING PDBCLASSES

pdbClasses2.o: pdbClasses2.cpp pdbClasses2.h
	$(CC) $(CFLAGS) $<
	'rm' libpdbClasses.a; ar -r libpdbClasses.a pdbClasses2.o

#auxiliary targets
clean:
	rm -f *.o core
