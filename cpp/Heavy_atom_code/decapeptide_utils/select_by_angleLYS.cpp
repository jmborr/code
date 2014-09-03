#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 5 ){
    system("clear");
    cout << "Usage: ./angleLYS.x pdbmov atf Am AM out\n";
    cout << "pdbmov : movie en PDB format\n" ;
    cout << "atf    : atom file in PDB format containing five ATOM entries\n";
    cout << "         The first three entries to create the plane. The\n";
    cout << "a2     : remaining two entries to determine beginning and end\n";
    cout << "a3     : points of a vector. We then determine angle between\n";
    cout << "a4     : the plane and the vector\n";
    cout << "       : If the five entries are a1,a2,a3,a4,a5, then\n"; 
    cout << "         angle = pi/2-acos((a5-a4) * [ (a1-a3) x (a2-a3) ])\n";
    cout << "Am     : minimum angle\n";
    cout << "AM     : maximum angle\n";
    cout << "out    : snapshots for which angle is in the Am-AM range\n";
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
  ifstream PDBMOV(argv[1]); /*cout<<"PDBMOV="<<argv[1]<<endl;*/
  listPDBatom list_ats(argv[2]); /*cout<<"five atoms file="<<argv[2]<<endl;*/
  double Am=atof(argv[3]), AM=atof(argv[4]);
  ofstream OUT(argv[5]);
  /*retrieve the five atom numbers*/

  if(list_ats.length()!=5){ 
    cout<<"list_ats.length()="<<list_ats.length()<<endl;
    cout<<"ERROR: we need five atoms in "<<argv[2]<<" !\n";
    return 1 ;
  }
  int ats_n[5];
  PDBatom at ;
  for(int i=1;i<=5;i++){
    at = list_ats.getAtomAtII(i); ats_n[i-1]=at.getSerialNumber();
  }
  PDBchain prot ;
  PDBvector v[5] ;
  double angle ;  /*string line ; getline(PDBMOV,line);cout<<line<<endl;*/
  while( PDBMOV>>prot ){
    for(int i=0;i<5;i++){
      at = prot.getAtomWithIndex(ats_n[i]) ;
      v[i] = at.getCoord( ) ;
    }
    v[1]=(v[0]-v[2])^(v[1]-v[2]); v[1].normalize() ;
    v[2]=v[4]-v[3]; v[2].normalize();
    angle = (PI/2-acos( v[1]*v[2] ))*180/PI ;/*angle with the plane*/
    if(angle>Am && angle <AM){ OUT<<prot<<"END \n" ; }
  }
  OUT.close();
  return 0 ;
}
