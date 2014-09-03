/*return amino acid indexes lacking the particular atom*/
#include "pdbClasses2.h"
/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./do_chain_lack_HN.x filename chainID\n" ;
    cout << "filename: absolute file name of the pdb file\n";
    cout << "chainID: ID of the chain to select\n";
    cout << "Output amino acid indexes lacking backbone hydrogen, except for the first amino acid and for PRO amino acids\n\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){
  if( test_input(argc, argv) == false ) { return 1 ; }
  PDBchains prot(argv[1]); /*cout<<prot;exit(0);*/
  string id(argv[2]);
  string name1(" H  "),name2("1H  ");
  PDBchain *chain=prot.pointToPDBchainFromChainID(id); 
  PDBamino_acid *aa;
  if(chain){
    int N=chain->length();
    for(int i=2;i<=N;i++){ /*i==1 has no HN*/
      aa=chain->pointToAminoFromIndex(i);
      if(!aa->is("PRO") && !aa->hasAtomWithName(name1) &&
	 !aa->hasAtomWithName(name2) ) /*PRO has no HN*/
	printf("%4d",i);
    }
  }
  return 0;
}
