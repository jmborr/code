/*return amino acid indexes lacking the particular atom*/
#include "pdbClasses2.h"
/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 4 ){
    system("clear");
    cout << "Usage: ./do_chain_lack_AtomName.x filename header type ATOM\n" ;
    cout << "filename: absolute file name of the pdb file\n";
    cout << "type: type of NOE ( HVSC_HVSC,HVBK_HVBK,F_HVSC_S_HVBK)\n"; 
    cout << "ATOM: 4-letter code for atom";
    cout << "\nreturn amino acid indexes lacking the particular atom\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){
  if( test_input(argc, argv) == false ) { return 1 ; }
  PDBchains prot(argv[1]); /*cout<<prot;exit(0);*/
  string header(argv[2]), id(argv[3]), atom(argv[4]);
  /*cout<<"header="<<header<<" id="<<id<<" |"<<atom<<"|"<<endl;exit(0);*/
  PDBchain *chain=prot.pointToPDBchainFromChainID(id);
  cout<<header<<id<<":";
  if(chain){
    chain->renumberFully(1);/*first amino acid has index 1*/
    chain->renumberFullyAtomSerialNumber(1);/*first atom has serial number 1*/
    int N=chain->length();
    for(int i=1;i<=N;i++)
      if(!chain->isAtomNameAtAminoAcidIndex(atom,i)) printf("%4d",i);
    cout<<endl;
  }
  else cout<<" no chain found\n";
  return 0;
}
