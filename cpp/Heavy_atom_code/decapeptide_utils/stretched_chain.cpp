#include<iostream>
using namespace std;
#include<fstream>
#include<string>
#include <cstdlib> 
#include"pdbClasses2.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( (number_of_arguments<2) || (number_of_arguments>4) ){
    system("clear");
    cout << "Usage: ./stretched_chain.x sequence chain.pdb [DPhi DPsi]\n" ;
    cout << "sequence: file with a single column of 3-letter code amino acids\n";
    cout << "chain.pdb: file in PDB format with the strectched chain\n";
    cout << "DPhi, Dpsi: the default chain is and stretched chain, ie, Phi=Psi=180. To change this enter additional rotations DPhi Dpsi.\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(  int argc, char ** argv ){

  /*check number of arguments*/
  if( test_input(argc, argv) == false ) { exit(0); }
  
  /*retrieve input*/
  int n_arg = argc -1 ;
  ifstream SEQ( argv[ 1 ] ) ; 
  ofstream OUT( argv[ 2 ] ) ;
  double DPhi = 0 ;
  double DPsi = 0 ;
  if( n_arg==4 ){
    DPhi = atof( argv[ 3 ] ) ;
    DPsi = atof( argv[ 4 ] ) ;
  }
  string aa3letter ;
  PDBchain chain ;
  while( getline( SEQ, aa3letter) ){
    /*cout << aa3letter <<endl;*/
    PDBamino_acid aa( aa3letter ) ;
    if( aa3letter == "PRO" ){
      cout << "It is PRO!\n" ;
      chain.insertAtBack( aa, 120, -60 ) ; 
    }
    else{ 
      chain.insertAtBack( aa, DPhi, DPsi ) ;
    }
  }
  chain.remove_all_hydrogens( ) ;
  chain.renumberFully( 1 ) ;
  chain.renumberFullyAtomSerialNumber( 1 ) ;
  SEQ.close( ) ;
  OUT << chain ;
  OUT.close( ) ;
  return 1 ;
}
