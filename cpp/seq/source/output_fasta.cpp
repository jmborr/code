#include "pdbClasses2.h"

/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 4 ){
    system("clear");
    cout << "Usage: ./output_seq.x pdb header ID seq\n" ;
    cout << "pdb: absolute file name of the pdb file\n";
    cout << "header: header string that heads the seq file\n";
    cout << "ID:one letter that specifies the chain (\"_\" also valid) \n";
    cout << "seq: name of output sequence in fasta format\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/

int main(int argc, char ** argv){

  if( test_input(argc, argv) == false ) { return 1 ; }

  /*cout<<"pdb="<<argv[1]<<endl;exit(0);*/
  PDBchains prot(argv[1]); /*cout<<prot;exit(0);*/
  string header(argv[2]),id(argv[3]);
  /*cout<<"id="<<id<<" header="<<header<<endl; exit(0);*/
  PDBchain *chain;
  chain=prot.pointToPDBchainFromChainID(id); 
  ofstream out(argv[4]);
  if(chain)chain->output_fasta_seq(out,header);
  else cout<<header<<"not found!\n";
  out.close();
  return 0;
}
