#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 5 ){
    system("clear");
    cout << "Usage: ./ses.x pdbmov atf Sm SM out\n";
    cout << "pdbmov : movie en PDB format\n" ;
    cout << "atf    : atom file in PDB format containing the atoms for\n";
    cout << "         which we compute the Solvent Excluded Surface\n";
    cout << "Sm     : minimal SES\n";
    cout << "SM     : maximal SES\n";
    cout << "out    : configurations with Sm<SES<SM\n";
    cout << "\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char* argv[]){

  /*check arguments*/
  if( test_input(argc, argv) == false ) { return 1 ; }
  
  double rC=1.87,rO=1.40,rN=1.65;
  string name;
  /*retrieve input*/
  ifstream PDBMOV(argv[1]);
  listPDBatom list_ats(argv[2]);
  double Sm=atof(argv[3]),SM=atof(argv[4]);
  ofstream OUT(argv[5]);
  int n_ats=list_ats.length();
  double r[n_ats];
  int ats_n[n_ats];
  PDBatom at ;
  for(int i=1;i<=n_ats;i++){
    at = list_ats.getAtomAtII(i);
    at.getAtomName(name);
    if(name.find("N")!=string::npos){r[i-1]=rN;}
    else if(name.find("O")!=string::npos){r[i-1]=rO;}
    else {r[i-1]=rC;}
    ats_n[i-1]=at.getSerialNumber();
  }
  PDBchain prot ;
  double ses,coord[3] ;
  ofstream JUNK;
  ifstream JUNK2;
  char command[200];
  system(command);
  /*string line ; getline(PDBMOV,line);cout<<line<<endl;*/
  while( PDBMOV>>prot ){
    JUNK.open("ses_junk");
    for(int i=0;i<n_ats;i++){
      prot.getCoordOfAtomWithAtomIndex(ats_n[i],coord) ;
      JUNK<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r[i]<<endl;
    }
    JUNK.close();
    sprintf(command,"msms -if ses_junk > ses_junk2;wait");
    system(command);
    sprintf(command,"perl -e '@_=<>;split(\" \",$_[$#{@_}-3]);print \"$_[$#{@_}]\\n\"' ses_junk2 > ses_junk;wait;");
    system(command);
    JUNK2.open("ses_junk"); 
    JUNK2>>ses; 
    JUNK2.close();
    if(ses>Sm && ses<SM){ OUT<<prot<<"END \n";}
  }
  
  system("rm -f ses_junk;rm -f ses_junk2");
  return 0 ;
}
