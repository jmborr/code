#include "pseudoPep.h"
#include "PDBLib.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;


int main(int argc, char* argv[]){
  if(argc<3){
    cout << "usage: command PDB nat_f" << endl;
    exit(1);
  }
  
  pseudoPep prot(argv[1], _CRYSTAL_);
  prot.makeIndex();
  int naa = prot.getLength();

  protein p(argv[1]);
  double size = 10.0*p.getRg();

  ofstream out(argv[2]);
  char buf[1024];
  sprintf(buf, "%ld %lf", naa, size);
  out << buf << endl;
  
  atom* theAtom;
  for(int i=0; i<naa; i++){
    theAtom = prot.getResidue(i)->getCA();
    sprintf(buf, "%ld %lf %lf %lf", 
	    theAtom->getIndex()-1,
	    theAtom->getR()->getX(),
	    theAtom->getR()->getY(),
	    theAtom->getR()->getZ());
    out << buf << endl;
  }
  out.close();
}
