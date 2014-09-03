#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 3 ){
    system("clear");
    cout << "Usage: ./sas.x pdbmov atf out\n";
    cout << "pdbmov : movie en PDB format\n" ;
    cout << "atf    : atom file in PDB format containing the atoms for\n";
    cout << "         which we compute the Solvent Accessible Surface\n";
    cout << "out    : output file\n";
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
  double sas,coord[3] ;
  ofstream JUNK;
  char command[200];
  sprintf(command,"'rm' %s >& /dev/null",argv[3]);
  system(command);
  /*string line ; getline(PDBMOV,line);cout<<line<<endl;*/
  while( PDBMOV>>prot ){
    JUNK.open("sas_junk");
    for(int i=0;i<n_ats;i++){
      prot.getCoordOfAtomWithAtomIndex(ats_n[i],coord) ;
      JUNK<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r[i]<<endl;
    }
    JUNK.close();
    exit(0);
    sprintf(command,"msms -if sas_junk > sas_junk2;wait");
    system(command);
    sprintf(command,"perl -e 'while(<>){if($_=~\"ANALYTICAL SURFACE AREA\"){$_=<>;$_=<>;chomp;split(\" \",$_);print \"$_[$#{@_}]\\n\"}}' sas_junk2 >> %s;wait;",argv[3]);
    system(command);
  }
  
  system("rm -f sas_junk;rm -f sas_junk2");
  return 0 ;
}
