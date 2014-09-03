/*return long-range contacts between hydrogen bonds*/
/*input is first absolute path to pdb file, second is the chain index*/
#include "pdbClasses2.h"
/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./outputHNHN.x file HNHNmap\n" ;
    cout << "filename: absolute file name of the pdb file\n";
    cout << "HNHNmap:output file containing HN-HN contacts\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){
  if( test_input(argc, argv) == false ) { return 1 ; }

  /*cout<<"filename="<<argv[1]<<" id="<<id<<endl;exit(0);*/
  double co=5.0;/*two heavy atoms do contact if relat. dist.<4.5Angstroms*/
  PDBchain chain(argv[1]);
  int N=chain.length();
  int skip=3;
  string type(" HH ");
  cont_sel cs=str2cs(type);
  double **map=alloc_array(N,N);
  int nc=chain.createContactMap(skip,cs,map,co);
  ofstream out(argv[2]);
  char buff[32];
  sprintf(buff,"%d\n",nc); out<<buff;
  print_symmetric_map(out,map,N);
  out.close();

  return 0;
}
