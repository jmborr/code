/*return long-range contacts between carbons of two methyl groups*/
/*input is first absolute path to pdb file*/
#include "pdbClasses2.h"
/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./outputMETHYLMETHYL.x file HNHNmap\n" ;
    cout << "filename: absolute file name of the pdb file\n";
    cout << "METHYLMETHYLmap:output file containing HN-HN contacts\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){
  if( test_input(argc, argv) == false ) { return 1 ; }

  /*cout<<"filename="<<argv[1]<<" id="<<id<<endl;exit(0);*/
  double co=5.0;/*two heavy atoms do contact if relat. dist.<5Angstroms*/
  PDBchain chain(argv[1]);
  int N=chain.length();
  int skip=3;
  double **map=alloc_array(N,N);
  int nc=chain.createContactMap(skip,_METHYLC_METHYLC_,map,co);
  ofstream out(argv[2]);
  char buff[32];
  sprintf(buff,"%d\n",nc); out<<buff;
  print_symmetric_map(out,map,N);
  out.close();

  return 0;
}
