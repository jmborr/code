/*return CA chain from a PDB*/
/*g++ -Wno-deprecated -o createCAcontactMap.x createCAcontactMap.cpp pdbClasses2.cpp*/
#include "pdbClasses2.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 3 ){
    system("clear");
    cout << "Usage: ./createCAhcin.x PDBfile mapFile co\n" ;
    cout << "PDBfile: absolute file name of the pdb file\n";
    cout << "CAfile:  output file containing CA contact map\n";
    cout << "co: cut-off, in angstroms (suggested 8.0A)\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){
  if( test_input(argc, argv) == false ) { return 1 ; }
  PDBchain chain(argv[1]);
  listPDBatom ca;
  chain.createCAchain(ca);
  ca.renumberFully(1);
  ca.renumberFullyAtomSerialNumber(1);
  ofstream out(argv[2]);
  double co=atof(argv[3]);
  ca.createContactMap( out, co ) ;
  return 0;
}
