#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 5 ){
    system("clear");
    cout << "Usage: ./sas.x pdbmov first_aa second_aa both_aa's out\n";
    cout << "pdbmov : movie en PDB format\n" ;
    cout << "first_aa: first side-chain atoms in PDB format \n";
    cout << "second_aa: scond side-chain atoms in PDB format \n";
    cout << "both_aa's: concatenation of the previous two files\n";
    cout << "out    : output file with SAS values for each frame of the previous movie\n";
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
  ofstream OUT(argv[5]);
  listPDBatom list_ats_first(argv[2]);
  listPDBatom list_ats_second(argv[3]);
  listPDBatom list_ats_both(argv[4]);
  int n_ats_first=list_ats_first.length();
  int n_ats_second=list_ats_second.length();
  int n_ats_both=list_ats_both.length();
  double r_first[n_ats_first], r_second[n_ats_second], r_both[n_ats_both];
  int ats_n_first[n_ats_first], ats_n_second[n_ats_second], ats_n_both[n_ats_both];
  PDBatom at ;
  for(int i=1;i<=n_ats_first;i++){
    at = list_ats_first.getAtomAtII(i);
    at.getAtomName(name);
    if(name.find("N")!=string::npos){r_first[i-1]=rN;}
    else if(name.find("O")!=string::npos){r_first[i-1]=rO;}
    else {r_first[i-1]=rC;}
    ats_n_first[i-1]=at.getSerialNumber();
  }
  for(int i=1;i<=n_ats_second;i++){
    at = list_ats_second.getAtomAtII(i);
    at.getAtomName(name);
    if(name.find("N")!=string::npos){r_second[i-1]=rN;}
    else if(name.find("O")!=string::npos){r_second[i-1]=rO;}
    else {r_second[i-1]=rC;}
    ats_n_second[i-1]=at.getSerialNumber();
  }
  for(int i=1;i<=n_ats_both;i++){
    at = list_ats_both.getAtomAtII(i);
    at.getAtomName(name);
    if(name.find("N")!=string::npos){r_both[i-1]=rN;}
    else if(name.find("O")!=string::npos){r_both[i-1]=rO;}
    else {r_both[i-1]=rC;}
    ats_n_both[i-1]=at.getSerialNumber();
  }
  PDBchain prot ;
  double sas, sas_2,coord[3] ;
  ofstream JUNK_OUT;
  ifstream JUNK_IN;
  char command[200];

  /*string line ; getline(PDBMOV,line);cout<<line<<endl;*/
  while( PDBMOV>>prot ){
    JUNK_OUT.open("sas_junk");
    for(int i=0;i<n_ats_both;i++){
      prot.getCoordOfAtomWithAtomIndex(ats_n_both[i],coord) ;
      JUNK_OUT<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r_both[i]<<endl;
    }
    JUNK_OUT.close();
    sprintf(command,"msms -if sas_junk > sas_junk2;wait");
    system(command);
    sprintf(command,"perl -e 'while(<>){if($_=~\"ANALYTICAL SURFACE AREA\"){$_=<>;$_=<>;chomp;split(\" \",$_);print \"$_[$#{@_}]\\n\"}}' sas_junk2 > sas_junk3");
    system(command);
    JUNK_IN.open("sas_junk3"); JUNK_IN>>sas;  JUNK_IN.close();
    if(sas<325){/*calculate FIRST and SECOND separately, then add*/
      JUNK_OUT.open("sas_junk");
      for(int i=0;i<n_ats_first;i++){
      prot.getCoordOfAtomWithAtomIndex(ats_n_first[i],coord) ;
      JUNK_OUT<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r_first[i]<<endl;
      }
      JUNK_OUT.close();
      sprintf(command,"msms -if sas_junk > sas_junk2;wait");
      system(command);
      sprintf(command,"perl -e 'while(<>){if($_=~\"ANALYTICAL SURFACE AREA\"){$_=<>;$_=<>;chomp;split(\" \",$_);print \"$_[$#{@_}]\\n\"}}' sas_junk2 > sas_junk3");
      system(command);
      JUNK_IN.open("sas_junk3"); JUNK_IN>>sas;  JUNK_IN.close();

      JUNK_OUT.open("sas_junk");
      for(int i=0;i<n_ats_second;i++){
      prot.getCoordOfAtomWithAtomIndex(ats_n_second[i],coord) ;
      JUNK_OUT<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r_second[i]<<endl;
      }
      JUNK_OUT.close();
      sprintf(command,"msms -if sas_junk > sas_junk2;wait");
      system(command);
      sprintf(command,"perl -e 'while(<>){if($_=~\"ANALYTICAL SURFACE AREA\"){$_=<>;$_=<>;chomp;split(\" \",$_);print \"$_[$#{@_}]\\n\"}}' sas_junk2 > sas_junk3");
      system(command);
      JUNK_IN.open("sas_junk3"); JUNK_IN>>sas_2;  JUNK_IN.close();
      sas += sas_2;
    }
    OUT<<sas<<endl;
  }
  system("rm -f sas_junk;rm -f sas_junk2;rm -f sas_junk3");
  return 0 ;
}
