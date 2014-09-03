#include "pdbClasses2.h"
/*
  For a given pdb chain, outputs a string of integers, as long as
  the chain length. '1' indicates that the particular residue is
  making a contact with some other residue that belongs to another
  chain in the same pdb file.

  For example, pdb1gh6.ent is composed of two chains. The C-terminal
  of 1gh6A makes contacts with 1gh6B. The command line is:

  interChainConts.x 1gh6A.pdb A 1gh6

  The output is:

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001302424032310000000

 last '4' means this particular residue in 1gh6A makes 4 contacts with some residues in 1gh6B
 NOTE: maximum reported number of contacts is 9, because we want one digit per residue.
*/
/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./interchainConts PDBfile chainID header\n" ;
    cout << "header: 4LETTER pdb code for the pdb containing the chain\n";
    cout << "chainID: one letter for chain identifier\n";
    cout << "\n\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char ** argv){
  if( test_input(argc, argv) == false ) { return 1 ; }
  int minl=40; /*neglect chains smaller than this*/
  int maxl=1000; /*neglecth chains longer than this*/
  double co=8.5; /*CA-CA distance cutoff*/
  string lib("/library/pdb");
  string pdbID(argv[1]);
  string chainID(argv[2]);
  string pdbf; 
  pdbf=lib+"/pdb"+pdbID+".ent";  /*cout<<pdbf<<endl;exit(1);*/
  PDBchains chains(pdbf);        /*cout<<chains<<endl;*/
  int l,nchains=chains.length(); /*cout<<nchains<<endl;*/
  PDBchain *chain;
  chain=chains.pointToPDBchainFromChainID(chainID);
  int nres=chain->length();       /*cout<<nres<<endl;*/
  double *doContact=new double[nres];
  for(int i=0;i<nres;i++) doContact[i]=0.0;
  double **cmap=new double*[nres]; /*interchain contact map*/
  cmap[0]=new double[nres*maxl];
  for(int i=1;i<nres;i++) cmap[i]=cmap[i-1]+maxl;
  PDBchain *pchain;
  for(int i=1;i<=nchains;i++){
    pchain=chains.pointToPDBchainFromIndex(i);   /*cout<<*pchain<<endl;*/
    l=pchain->length();
    if(pchain->printChainID() != chainID) /*don't compare to itself or NRM isomers*/
       if(l>minl){ /*neglect comparison to small chains*/
	 for(int i=0;i<nres*maxl;i++){ cmap[0][i]=0.0;} /*erase contact map*/
	 chain->createCAcontactMap( *pchain,cmap,co);
	 for(int i=0;i<nres;i++)
	   for(int j=0;j<l;j++){
	     doContact[i]+=cmap[i][j];
	     if(doContact[i]>9) doContact[i]=9; /*maximum reported number of contacts*/
	   }
       }
  }

  cout<<"> "<<pdbID<<chainID<<endl;
  for(int i=0;i<nres;i++) cout<<int(doContact[i]);
  cout<<endl;
  return 0;
}
