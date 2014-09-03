#include "pdbClasses2.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 3 ){
    system("clear");
    cout << "Usage: ./output_single_chain.x pdb ID file\n" ;
    cout << "pdb: absolute file name of the pdb file\n";
    cout << "ID:one letter that specifies the chain (\"_\" also valid) \n";
    cout << "file: name of output file containing only chain with ID\n";
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
  ofstream out(argv[3]);
  if(chain){
    chain->renumberFully(1);/*first amino acid has index 1*/
    chain->renumberFullyAtomSerialNumber(1);/*first atom has serial number 1*/
    out<<*chain;
    out<<"END\n";
  }
  else cout<<id<<" of "<<argv[1]<<"not found!\n";
  out.close();
  
  return 0;
}
