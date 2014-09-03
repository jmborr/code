#include<iostream>
using namespace std;
#include<fstream>
#include"pdbClasses2.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 3 ){
    system("clear");
    cout << "Usage: ./txt2pdb_Heavy_Atom.x txt pdb out\n" ;
    cout << "txt: initial configuration file\n";
    cout << "pdb: any PDB snapshot of the protein\n";
    cout << "out: output configuration in PDB format\n" ;
    cout << "pre-relaxed.txt: configuration suitable for DMD, with steric clashes\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/

int main(  int argc, char ** argv ){

  /*check number of arguments*/
  if( test_input(argc, argv) == false ) { exit(0); }

  /*retrieve input*/
  ifstream TXT( argv[1] ) ;
  ifstream PDB( argv[ 2 ] ) ;
  ofstream OUT( argv[ 3 ] ) ;
  /*store PDB in memory*/
  PDBchain chain ;   PDB >> chain ;   PDB.close( ) ;
  chain.remove_all_hydrogens( ) ;
  chain.renumberFully( 1 ) ;
  chain.renumberFullyAtomSerialNumber( 1 ) ; 
  /*input coordinates from TXT*/
  int l=chain.length();
  /*scroll TXT until system size*/
  string line;
  double box;
  do{getline(TXT,line);}while(line.find("A.SYSTEM SIZE")==string::npos);
  int p=TXT.tellg(); getline(TXT,line);
  if(line.find("//")==string::npos){ TXT.seekg(p);}
  TXT>>box>>box>>box;
  /*scroll TXT until coordinates*/
  do{getline(TXT,line);}while(line.find("H.LIST OF ATOMS")==string::npos);
  p=TXT.tellg(); getline(TXT,line);
  if(line.find("//")==string::npos){ TXT.seekg(p);}
  /*change coordinates of "chain"*/
  int a,b;  double **r=NULL, v[3];  r=alloc_array(l,3);
  for(int i=0;i<l;i++){TXT>>a>>b>>r[i][0]>>r[i][1]>>r[i][2]>>v[0]>>v[1]>>v[2];}
  TXT.close();
  remove_bound_cond(r,l,box);
  for(int i=0;i<l;i++){chain.changeCoordOfAtomWithIndex(i+1,r[i]);}
  /*output new conformation*/
  OUT<<chain<<"END\n";

  return 0;
}
