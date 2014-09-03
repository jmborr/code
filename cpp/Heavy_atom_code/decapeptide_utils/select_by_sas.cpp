#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 7 ){
    system("clear");
    cout << "Usage: ./sas.x pdbmov val lys both Smin Smax out\n";
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
  ofstream OUT(argv[7]);
  listPDBatom list_ats_val(argv[2]);
  listPDBatom list_ats_lys(argv[3]);
  listPDBatom list_ats_both(argv[4]);
  int n_ats_val=list_ats_val.length();
  int n_ats_lys=list_ats_lys.length();
  int n_ats_both=list_ats_both.length();
  double r_val[n_ats_val], r_lys[n_ats_lys], r_both[n_ats_both];
  int ats_n_val[n_ats_val], ats_n_lys[n_ats_lys], ats_n_both[n_ats_both];
  PDBatom at ;
  for(int i=1;i<=n_ats_val;i++){
    at = list_ats_val.getAtomAtII(i);
    at.getAtomName(name);
    if(name.find("N")!=string::npos){r_val[i-1]=rN;}
    else if(name.find("O")!=string::npos){r_val[i-1]=rO;}
    else {r_val[i-1]=rC;}
    ats_n_val[i-1]=at.getSerialNumber();
  }
  for(int i=1;i<=n_ats_lys;i++){
    at = list_ats_lys.getAtomAtII(i);
    at.getAtomName(name);
    if(name.find("N")!=string::npos){r_lys[i-1]=rN;}
    else if(name.find("O")!=string::npos){r_lys[i-1]=rO;}
    else {r_lys[i-1]=rC;}
    ats_n_lys[i-1]=at.getSerialNumber();
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
  double minsas=atof(argv[5]);
  double maxsas=atof(argv[6]); 
  double sas, sas_2,coord[3] ;
  ofstream JUNK_OUT;
  ifstream JUNK_IN;
  char command[200];

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
    if(sas<325){/*calculate VAL and LYS separately, then add*/
      JUNK_OUT.open("sas_junk");
      for(int i=0;i<n_ats_val;i++){
      prot.getCoordOfAtomWithAtomIndex(ats_n_val[i],coord) ;
      JUNK_OUT<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r_val[i]<<endl;
      }
      JUNK_OUT.close();
      sprintf(command,"msms -if sas_junk > sas_junk2;wait");
      system(command);
      sprintf(command,"perl -e 'while(<>){if($_=~\"ANALYTICAL SURFACE AREA\"){$_=<>;$_=<>;chomp;split(\" \",$_);print \"$_[$#{@_}]\\n\"}}' sas_junk2 > sas_junk3");
      system(command);
      JUNK_IN.open("sas_junk3"); JUNK_IN>>sas;  JUNK_IN.close();

      JUNK_OUT.open("sas_junk");
      for(int i=0;i<n_ats_lys;i++){
      prot.getCoordOfAtomWithAtomIndex(ats_n_lys[i],coord) ;
      JUNK_OUT<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r_lys[i]<<endl;
      }
      JUNK_OUT.close();
      sprintf(command,"msms -if sas_junk > sas_junk2;wait");
      system(command);
      sprintf(command,"perl -e 'while(<>){if($_=~\"ANALYTICAL SURFACE AREA\"){$_=<>;$_=<>;chomp;split(\" \",$_);print \"$_[$#{@_}]\\n\"}}' sas_junk2 > sas_junk3");
      system(command);
      JUNK_IN.open("sas_junk3"); JUNK_IN>>sas_2;  JUNK_IN.close();
      sas += sas_2;
    }
    if(minsas<sas && sas<maxsas){OUT<<prot<<"END\n";}
  }
  system("rm -f sas_junk;rm -f sas_junk2;rm -f sas_junk3");

  
  return 0 ;
}
