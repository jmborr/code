#include<iostream>
using namespace std;
#include<fstream>
#include<string>
#include <cstdlib> 
#include<fstream>
#include"pdbClasses2.h"
#include"atom_param.h"
#include"amino_acid_param.h"
#include"miscellanea.h"
#include"random.h"
/*=====================================================*/
void printSYS_SIZE( ostream& out, double &size, int dimension=3 ){
  cout << "SYS_SIZE ...\n" ;
  out << txtKeyWords[SYS_SIZE] << endl;
  for(int i=0; i<dimension; i++){ out << size << " "; }
  out << endl;
}
/*=====================================================*/
void printNUM_ATOMS(ostream& out, PDBchain& chain ){
  cout << "NUM_ATOMS ...\n" ;
  out << txtKeyWords[NUM_ATOMS] << endl ;
  out << chain.numberOfAtoms( ) << endl ;
}
/*=====================================================*/
void printATOM_TYPE( ifstream &DAT, ostream& out, double *hcr_list ){
  cout << "ATOM_TYPE ...\n" ;
  out << txtKeyWords[TYPE_ATOMS] << endl;
  dmd_atom_t type ;
  string amino_and_atom ;
  char cbuf[300] ;
  for( int i=1; i<=n_dmd_atom_t; i++){
    sprintf( cbuf, "%3d ", i ) ; 
    out << cbuf ;
    type = static_cast<dmd_atom_t>( i ) ; 
    amino_and_atom = dmd_atom_t2string( type ) ;
    print_atom_type_line( cbuf, DAT, amino_and_atom ) ;
    hcr_list[ i ] = out_hard_core_radius( DAT, amino_and_atom ) ;
    out << cbuf ;
  }
}
/*=====================================================*/
void  printNONEL_COL( ifstream &DAT, ofstream &out, double *hcr_list ){
  cout << "NONEL_COL ...\n" ;
  out << txtKeyWords[NONEL_COL] << endl;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  double f = get_inner_VW_factor( DAT ) ;
  char cbuf[256] ;
  double hcd ;
  for( int i=1; i<n_dmd_atom_t; i++){
    for( int j=i; j<=n_dmd_atom_t; j++){
      sprintf( cbuf, "%3d %3d ", i, j ) ;
      out << cbuf ;
      hcd = (hcr_list[ i ] + hcr_list[ j ]) * f ;
      sprintf( cbuf, " %9.6f %9.6f 1.0 %9.6f 1.0 %9.6f 1.0 %9.6f 1.0 %9.6f 1.0 %9.6f 1.0 %9.6f 1.0\n",
	       hcd*0.50, hcd*0.60, hcd*0.65, hcd*0.70,
	       hcd*0.75, hcd*0.80, hcd*0.85, hcd) ;
      out << cbuf ;
    }
  }
  DAT.seekg( pos ) ;
}
/*=====================================================*/
void  printLINKED_PAIRS( ifstream &DAT, ofstream &DMD ){
  cout << "LINKED_PAIRS ...\n" ;
  DMD << txtKeyWords[LINK_PAIRS] << endl;
  int pos = DAT.tellg( ), i, j ;
  char cbuf[256], cbuf2[256] ; ;
  dmd_atom_t type1, type2 ;
  string start, end, link, aa_at1, aa_at2 ;
  link = start = "#" + txtKeyWords[LINK_PAIRS] ; 
  end = "#END " + txtKeyWords[LINK_PAIRS] ; /*cout<<"end="<<end<<endl;*/
  go_to( DAT, start ) ;/* cout<<"start="<<start<<endl;*/
  while( link != end ){
    if(out_prerelaxed_linked_pair( link, aa_at1, aa_at2, cbuf ) ){
      type1 = string2dmd_atom_t( aa_at1 ) ;
      type2 = string2dmd_atom_t( aa_at2 ) ;
      i = static_cast<int>( type1 ) ;
      j = static_cast<int>( type2 ) ;
      if( i==0 ){ cout <<"|"<<aa_at1<<"| has no type!\n"; exit(1); }
      if( j==0 ){ cout <<"|"<<aa_at2<<"| has no type!\n"; exit(1); }
      sprintf( cbuf2, "%3d %3d  %s\n", i, j, cbuf ) ;
      DMD << cbuf2 ;
    }
    getline( DAT, link ) ;
  }
  DAT.seekg( pos ) ;
}
/*=====================================================*/
void printLIST_ATOMS( PDBchain &chain, double &box, ofstream &DMD ){
  cout << "LIST ATOMS ...\n" ;
  DMD << txtKeyWords[LIST_ATOMS] << endl;
  double boxH = box/2.0 ;
  int nAts ;
  string amino_and_atom ;
  dmd_atom_t type ;
  PDBatom at ;
  PDBvector CMfinal( boxH, boxH, boxH ) ;
  randomGenerator v(_RAN2_, -1);
  double *r = new double[3] ;
  char cbuf[256] ;
  chain.translate( CMfinal - chain.get_CA_CM( ) ) ;
  nAts = chain.numberOfAtoms( ) ;
  for( int i=1; i<=nAts; i++ ){
    at = chain.getAtomWithIndex( i ) ;
    at.getCoordToDoubleArrayII( r ) ;
    amino_and_atom = at.get_resName_nameII( ) ;
    type = string2dmd_atom_t( amino_and_atom ) ;
    sprintf( cbuf, "%4d %3d %9.5f %9.5f %9.5f %9.4f %9.4f %9.4f\n",
	     i, static_cast<int>( type ), r[0], r[1], r[2],
	     v.nextGauss( ), v.nextGauss( ), v.nextGauss( ) ) ;
    DMD << cbuf ;
  }
}
/*=====================================================*/
void printLIST_OF_BONDS( PDBchain &chain, ifstream &DAT, ofstream &DMD ){
  cout << "LIST OF BONDS ...\n" ;
  DMD << txtKeyWords[LIST_BONDS] << endl;
  int nAts = chain.numberOfAtoms( ), l=chain.length( ) ; 
  int nBonds_same, nBonds_next, resSeq, resSeq2, index, index2 ;
  string aa_at, resName, resName2, *friends_same=NULL, *friends_next=NULL  ;
  string companion, same_or_next, at2_name ;
  char cbuf[32] ;
  PDBatom at ;
  for( index=1; index<=nAts; index++ ){
    /*cout<<"l="<<l<<", index="<<index<<endl;*/
    at=chain.getAtomWithIndex(index);  /*cout<<at<<endl;*/
    resSeq = at.getResSeq( )   ;  resSeq2 = resSeq+1 ;
    resName = at.getResName( ) ;  
    if( resSeq < l ){ resName2 = chain.getNameFromIndex( resSeq2 ) ; }
    aa_at = at.get_resName_nameII( ) ;
    /*cout<<"resSeq="<<resSeq<<" resName="<<resName<<endl;*/
    friends_same=out_all_same_permanent_bonds(aa_at,resName,DAT,nBonds_same);
    if( friends_same ){
      for( int i=0; i<nBonds_same; i++ ){
	at2_name = friends_same[ i ].substr( 4, 4 ) ;
	index2 = chain.getSerialNumberOfAtomNameAtIndex( at2_name, resSeq );
	sprintf( cbuf, "%3d %3d\n", index, index2 ) ;
	DMD << cbuf ;
      }
      delete [ ] friends_same ;  friends_same = NULL ; 
    }
    if( resSeq < l ){/*there is no next amino acid when resSeq==l*/
      /*cout<<"resSeq="<<resSeq<<endl;*/
      friends_next = out_all_next_permanent_bonds( aa_at, resName, resName2, DAT, nBonds_next ) ;
      if( friends_next ){
	for( int i=0; i<nBonds_next; i++ ){
	  at2_name = friends_next[ i ].substr( 4, 4 ) ;
	  index2 = chain.getSerialNumberOfAtomNameAtIndex( at2_name, resSeq2 );
	  /*cout<<friends_next[ i ]<<" index2="<<index2<<endl;*/
	  sprintf( cbuf, "%3d %3d\n", index, index2 ) ;
	  DMD << cbuf ; 
	}
	delete [ ] friends_next ; friends_next = NULL ; 
      }
    }
  }/*matches for( index=1; index<=nAts; index++ )*/
}
/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 4 ){
    system("clear");
    cout << "Usage: ./pdb2txt_relax.x pre-relaxed.pdb box-size force-field pre-relaxed.txt\n" ;
    cout << "pre-relaxed.pdb: initial PDB file of the peptide\n";
    cout << "box-size: dimensions ( in Angstroms) of the simulation cube\n";
    cout << "force-field: field with info on atoms (mass, radii, bonded and non-bonded interactions\n" ;
    cout << "pre-relaxed.txt: configuration suitable for DMD, with steric clashes\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(  int argc, char ** argv ){

  /*check number of arguments*/
  if( test_input(argc, argv) == false ) { exit(0); }
  
  /*retrieve input*/
  ifstream PDB( argv[ 1 ] ) ;
  PDBchain chain ;   PDB >> chain ;   PDB.close( ) ;
  double cubeSide = atof( argv[2] ) ;
  ifstream DAT( argv[ 3 ] ) ;
  ofstream DMD( argv[ 4 ] ) ;

  chain.remove_all_hydrogens( ) ;
  chain.renumberFully( 1 ) ;
  chain.renumberFullyAtomSerialNumber( 1 ) ;

  double *hcr_list = new double[ n_dmd_atom_t + 1] ;
  printSYS_SIZE( DMD, cubeSide ) ;
  printNUM_ATOMS( DMD, chain ) ;
  printATOM_TYPE( DAT, DMD, hcr_list ) ;
  printNONEL_COL( DAT, DMD, hcr_list ) ;
  printLINKED_PAIRS( DAT, DMD ) ;
  printLIST_ATOMS( chain, cubeSide, DMD ) ;  /*cout<<chain ;*/
  printLIST_OF_BONDS( chain, DAT, DMD ) ;

  return 0 ;
}
