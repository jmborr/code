#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

static double rC=1.87,rO=1.40,rN=1.65;

/*=================================================================*/
void initialize(listPDBatom &list,double *r,int *ats_n,int n_ats, 
		double *r_both,int *ats_n_both,int shift         ){
  string name;
  PDBatom at ;
  for(int i=1;i<=n_ats;i++){
    at = list.getAtomAtII(i);
    at.getAtomName(name);
    if(name.find("N")!=string::npos){r_both[shift+i-1]=r[i-1]=rN;}
    else if(name.find("O")!=string::npos){r_both[shift+i-1]=r[i-1]=rO;}
    else {r_both[shift+i-1]=r[i-1]=rC;}
    ats_n_both[shift+i-1]=ats_n[i-1]=at.getSerialNumber();
  }
}
/*=================================================================*/
double get_sas(PDBchain &prot, double *r, int *ats_n, int n_ats){
  ifstream JUNK_IN;
  ofstream JUNK_OUT;
  char command[200];
  double coord[3],sas;
  
  JUNK_OUT.open("sas_junk");
  for(int i=0;i<n_ats;i++){
    prot.getCoordOfAtomWithAtomIndex(ats_n[i],coord) ;
    JUNK_OUT<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r[i]<<endl;
    cout<<coord[0]<<"  "<<coord[1]<<"  "<<coord[2]<<"  "<<r[i]<<endl;
  }
  JUNK_OUT.close();
  sprintf(command,"msms -if sas_junk > sas_junk2;wait");
  system(command);
  sprintf(command,"perl -e 'while(<>){if($_=~\"ANALYTICAL SURFACE AREA\"){$_=<>;$_=<>;chomp;split(\" \",$_);print \"$_[$#{@_}]\\n\"}}' sas_junk2 > sas_junk3");
  system(command);
  JUNK_IN.open("sas_junk3"); JUNK_IN>>sas;  JUNK_IN.close();
  return sas;
}
/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 4 ){
    system("clear");
    cout << "Usage: ./sas.x pdbmov first_aa second_aa out\n";
    cout << "pdbmov : movie en PDB format\n" ;
    cout << "first_aa: first side-chain atoms in PDB format \n";
    cout << "second_aa: scond side-chain atoms in PDB format \n";
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
  
  /*retrieve input*/
  ifstream PDBMOV(argv[1]);
  listPDBatom list_ats_first(argv[2]);
  listPDBatom list_ats_second(argv[3]);
  int n_ats_first=list_ats_first.length();
  int n_ats_second=list_ats_second.length();
  int n_ats_both= n_ats_first+n_ats_second;
  double r_first[n_ats_first], r_second[n_ats_second], r_both[n_ats_both];
  int ats_n_first[n_ats_first],ats_n_second[n_ats_second],ats_n_both[n_ats_both];
  ofstream OUT(argv[4]);

  initialize(list_ats_first,r_first,ats_n_first,n_ats_first,
	     r_both,ats_n_both,0);
  initialize(list_ats_second,r_second,ats_n_second,n_ats_second,
	     r_both,ats_n_both,n_ats_first);

  PDBchain prot ;
  double sas_all, sas_first, sas_second;

  /*string line ; getline(PDBMOV,line);cout<<line<<endl;*/
  string line;
  while( PDBMOV>>prot ){
    sas_all=get_sas(prot,r_both,ats_n_both,n_ats_both);
    cout<<"sas_all="<<sas_all<<endl;
    sas_first=get_sas(prot,r_first,ats_n_first,n_ats_first);
    cout<<"sas_first="<<sas_first<<endl;
    if(sas_all*0.99<sas_first && sas_all*1.01>sas_first){
      /*first and second are separated. Then msms can't just calculate
	SAS of first, then SAS of second, then add them*/
      sas_second=get_sas(prot,r_second,ats_n_second,n_ats_second);
      sas_all+=sas_second;
      cout<<"sas_second="<<sas_second<<endl;
    }
    OUT<<sas_all<<endl;
    getline(PDBMOV,line);
  }
  
  system("rm -f sas_junk;rm -f sas_junk2;rm -f sas_junk3");
  return 0 ;
}
