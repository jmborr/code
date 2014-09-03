#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./angleLYSII.x pdbmov atf\n";
    cout << "pdbmov : movie en PDB format\n" ;
    cout << "atf    : atom file in PDB format containing six ATOM entries\n";
    cout << "         ordered as a1, a2, a3, a4, a5, a6\n";
    cout << "         angle = sign*acos( (a1-a2)*(a3-a2) )\n" ;
    cout << "         sign determined by the sign of the scalar product\n";
    cout << "          between vector (a3-a2) and (a6-a4)^(a6-a5)\n";
    cout << "\n";
    cout << "\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char* argv[]){

  /*check arguments*/
  if( test_input(argc, argv) == false ) { return 1 ; }

  /*retrieve input*/
  ifstream PDBMOV(argv[1]);
  listPDBatom list_ats(argv[2]);
  /*retrieve the five atom numbers*/
  if(list_ats.length()!=6){ 
    cout<<"ERROR: we need five atoms in "<<argv[2]<<" !\n";
    return 1 ;
  }
  int ats_n[6];
  PDBatom at ;
  for(int i=1;i<=6;i++){
    at = list_ats.getAtomAtII(i);
    ats_n[i-1]=at.getSerialNumber();
  }
  PDBchain prot ;
  PDBvector v[6] ;
  double angle ;
  /*string line ; getline(PDBMOV,line);cout<<line<<endl;*/
  while( PDBMOV>>prot ){
    for(int i=0;i<6;i++){
      at = prot.getAtomWithIndex(ats_n[i]) ;
      v[i] = at.getCoord( ) ; 
    }
    v[0]=v[0]-v[1]; v[0].normalize() ;
    v[2]=v[2]-v[1]; v[2].normalize();
    v[4]=(v[5]-v[3])^(v[5]-v[4]); v[4].normalize();
    angle = acos( v[0]*v[2] ) ;
    if(v[2]*v[4]<0){angle *= -1;}
    cout<<angle*180/PI<<endl;
  }

  return 0 ;
}
