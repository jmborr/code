#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>
#include<cstring>
/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./turn_propensities.x pdbmov out\n";
    cout << "pdbmov : movie en PDB format\n" ;
    cout << "out    : \"angle\" printed to \"out\" file\n";
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

  PDBatom prev_2,prev_1,curr,next_1,next_2 ;
  PDBchain prot;
  listPDBatom listCA;
  int l;
  double r;
  char buff[5120],buff2[128];
  double angle ;  /*string line ; getline(PDBMOV,line);cout<<line<<endl;*/
  while( PDBMOV>>prot ){
    sprintf(buff,"");
    prot.createCAchain(listCA);
    listCA.renumberFully(1);
    l=listCA.length();
    prev_2=listCA.getAtomAtII(1);
    prev_1=listCA.getAtomAtII(2);
    curr  =listCA.getAtomAtII(3);
    next_1=listCA.getAtomAtII(4);
    for(int i=5;i<=l;i++){
      next_2=listCA.getAtomAtII(i);
      /*cout<<prev<<endl<<curr<<endl<<next<<endl;*/
      r=curr.radius(prev_2,next_2);
      /*sprintf(buff2,"%6.2lf",r);*/
      sprintf(buff2,"%4d %6.2lf\n",i-2,r);
      strncat(buff,buff2,strlen(buff2));
      prev_2=prev_1;  prev_1=curr;  curr=next_1;  next_1=next_2;
    }
    OUT<<buff<<endl;
  }
  OUT.close();
  return 0;
}
