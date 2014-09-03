/*assume input PDB file contains the chain of interest as first chain*/
#include "pdbClasses2.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./add_HN.x pdb file\n" ;
    cout << "pdb: absolute file name of the pdb file\n";
    cout << "file: output file including hydrogens in the backbone\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){

  if( test_input(argc, argv) == false ) { return 1 ; }

  PDBchain chain(argv[1]); /*cout<<prot;exit(0);*/
  ofstream out(argv[2]);
  if(!chain.isEmpty()){
    chain.insert_H( );
    out<<chain;
    out<<"END\n";
  }
  else cout<<"ERROR: no chain found!\n";
  out.close(); 
  return 0;
}
