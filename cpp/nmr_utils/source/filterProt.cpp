/*return 0 iff heavy atom contact map for side chains has more than N/4
  long-range native contacts (long-range iff (i,j>i+4)*/
/*input is first absolute path to pdb file, second is the chain index*/
#include "pdbClasses2.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./filterProt.x filename ID\n" ;
    cout << "filename: absolute file name of the pdb file\n";
    cout << "ID:one letter that specifies the chain (\"_\" also valid) \n";
    cout << "return 0 iff heavy atom contact map for side chains has more\nthan N/4 long-range native contacts (long-range iff (i,j>i+4)\ninput is first absolute path to pdb file, second is the chain\nindex\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){

  if( test_input(argc, argv) == false ) { return 1 ; }

  /*cout<<"filename="<<argv[1]<<" id="<<id<<endl;exit(0);*/
  double co=4.5;/*two heavy atoms do contact if relat. dist.<4.5Angstroms*/
  string id(argv[2]);
  PDBchains prot(argv[1]);
  PDBchain *chain=prot.pointToPDBchainFromChainID(id);
  int N=chain->length();
  int skip=4;
  int nc=chain->number_of_contacts(skip,_HV_HV_,co);
  int x=4;
  if(nc/x>=N) return 0;
  return 1;
}
