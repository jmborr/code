#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

#define NO_m  3.12
#define NO_M  3.20
#define NC_m  3.80
#define NC_M  4.35
#define CO_m  3.60
#define CO_M  4.12
#define CA0_m 3.60  /*minimal CA-O distance*/
#define CAO_M 4.17

double no_m=NO_m*NO_n;
double no_M=NO_M*NO_M;
double nc_m=NC_m*NC_m;
double nc_M=NC_M*NC_M;
double co_m=CO_m*CO_m;
double co_M=CO_M*CO_M;
double cao_m=CAO_m*CAO_m;
double cao_M=CAO_M*CAO_M;

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./n_HB.x pdbmov out\n";
    cout << "pdbmov : movie en PDB format\n" ;
    cout << "out    : number of hydrogen bonds per frame printed to \"out\" file\n";
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
  ofstream OUT(argv[2]);
  /*retrieve the five atom numbers*/

  PDBchain prot;
  PDBamino_acid *aa1,*aa2, *aaX;
  PDBatom *N1,*O1,*C1,*CA1,*N2,*O2,*C2,*CA2,*CX;
  int n_HB;
  while( PDBMOV>>prot ){
    n_HB=0;
    prot.renumberFully(1);
    for( int i=1;i<prot.length()-1;i++){
      aa1=prot.getAminoFromIndex(i);
      N1 =aa1->getAtomFromName(" N  ") ;
      O1 =aa1->getAtomFromName(" O  ") ;
      C1 =aa1->getAtomFromName(" C  ") ;
      CA1=aa1->getAtomFromName(" CA ");
      for( int j=i+1;j<prot.length();j++){
	aa2=prot.getAminoFromIndex(j);
	N2 =aa2->getAtomFromName(" N  ") ;
	O2 =aa2->getAtomFromName(" O  ") ;
	C2 =aa2->getAtomFromName(" C  ") ;
	CA2=aa2->getAtomFromName(" CA ");
	if(i>1 && N1->d2(*O2)>no_m && N1->d2(*O2)<no_M && 
	   N1->d2(*C2)>nc_m &&  N1->d2(*C2)<nc_M &&
	   CA1->d2(*O2)>cao_m && CA1->d2(*O2)<cao_M){
	  aaX=prot.getAminoFromIndex(i-1);
	  CX=aaX->getAtomFromName(" C  ");
	  if(CX->d2(O2)>co_m && CX->d2(O2)<co_M){ n_HB++; }
	  delete aaX; delete CX;
	}
	if(O1->d2(N2)>no_m && O1->d2(N2)<no_M &&
	   O1->d2(CA2)>cao_m && O1->d2(CA2)<cao_M &&
	   C1->d2(N2)>nc_m  && C1->d2(N2)<nc_M ){
	  aaX=prot.getAminoFromIndex(j-1);
	  CX=aaX->getAtomFromName(" C  ");
	  if(CX->d2(O1)>co_m && CX->d2(O1)<co_M){ n_HB++; }
	  delete aaX; delete CX;
	}
	delete aa2; delete N2; delete O2; delete C2; delete CA2;
      }
      delete aa1; delete N1; delete O1; delete C1; delete CA1;
    }
  }
  OUT.close();
  return 0 ;
}
