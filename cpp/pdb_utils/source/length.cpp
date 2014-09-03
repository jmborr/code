#include "pdbClasses2.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./length.x pdb ID\n" ;
    cout << "pdb: absolute file name of the pdb file\n";
    cout << "ID:one letter that specifies the chain (\"_\" also valid) \n";
    cout << "output to STDOUT the lenght\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){

  if( test_input(argc, argv) == false ) { return 1 ; }

  PDBchains prot(argv[1]); /*cout<<prot;exit(0);*/
  string id(argv[2]);
  PDBchain *chain=prot.pointToPDBchainFromChainID(id);
  cout<<chain->length()<<endl;

  return 0;
}
