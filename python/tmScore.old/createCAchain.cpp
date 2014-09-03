/*return CA chain from a PDB*/
#include "createCAchain.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./createCAhcin.x PDBfile CAfile\n" ;
    cout << "PDBfile: absolute file name of the pdb file\n";
    cout << "CAfile:  output file containing CA atoms only\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){
  if( test_input(argc, argv) == false ) { return 1 ; }
  createCAchain(argv[1],argv[2]);
  return 0;
}
