#include "pdbClasses2.h"

int createCAchain(char *inpdb, char *outpdb){
  PDBchain chain(inpdb);
  listPDBatom ca;
  chain.createCAchain(ca);
  ca.renumberFully(1);
  ca.renumberFullyAtomSerialNumber(1);
  ofstream out(outpdb);
  out<<ca;
}
/*=====================================================*/
double **getCAcoords(char *inpdb){
  PDBchain chain(inpdb);
  listPDBatom ca;
  chain.createCAchain(ca);
  return ca.dump_coordinates_to_array();
}
/*=====================================================*/
int getCAcoords2(double **x, char *inpdb){
  PDBchain chain(inpdb);
  listPDBatom ca;
  chain.createCAchain(ca);
  return ca.dump_coordinates_to_array2(x);
}
/*=====================================================*/
