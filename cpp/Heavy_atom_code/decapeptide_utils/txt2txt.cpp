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
void printSYS_SIZE( ifstream &RLX, ofstream &DMD ){
  DMD << txtKeyWords[SYS_SIZE] << endl;
  int pos = RLX.tellg( ) ;
  RLX.seekg( 0 ) ;
  string line ;
  double box[3] ;
  char cbuf[128] ;
  while( getline( RLX, line ) ){
    if( line.find(txtKeyWords[SYS_SIZE]) != string::npos ){ break ; }
  }
  RLX >> box[ 0 ] >> box[ 1 ] >> box[ 2 ] ; 
  sprintf( cbuf, "%9.4f %9.4f %9.4f\n", box[0], box[1], box[2]) ;
  DMD << cbuf ;
  RLX.seekg( pos ) ;
  cout << "SYS_SIZE ...\n" ;
}
/*=====================================================*/
void  printNUM_ATOMS( PDBchain& chain, ofstream &DMD ){
  cout << "NUM_ATOMS ...\n" ;
  DMD << txtKeyWords[NUM_ATOMS] << endl ;
  DMD << chain.numberOfAtoms( ) << endl ;
}
/*=====================================================*/
void printATOM_TYPE( ifstream &DAT, ostream& DMD, double *hcr_list ){
  cout << "ATOM_TYPE ...\n" ;
  DMD << txtKeyWords[TYPE_ATOMS] << endl;
  dmd_atom_t type ;
  string amino_and_atom ;
  char cbuf[300] ;
  for( int i=1; i<=n_dmd_atom_t; i++){
    sprintf( cbuf, "%3d ", i ) ; 
    DMD << cbuf ;
    type = static_cast<dmd_atom_t>( i ) ;/*amino_param.h*/
    amino_and_atom = dmd_atom_t2string( type ) ;/*amino_param.cpp*/
    print_atom_type_line( cbuf, DAT, amino_and_atom ) ;/*atom_param.cpp*/
    hcr_list[ i ] = out_hard_core_radius( DAT, amino_and_atom ) ;
    DMD << cbuf ;
  }
}
/*=====================================================*/
void printNONEL_COL( ifstream &DAT, ostream& DMD, double *hcr_list ){
  cout << "NONEL_COL ...\n" ;
  DMD << txtKeyWords[NONEL_COL] << endl;
  string nonel_col[n_dmd_atom_t+1][n_dmd_atom_t+1] ;
  dmd_atom_t type1, type2 ;
  string stype1, stype2 ;
  char cbuf[256], cbuf2[256] ;
  for( int i=1; i<=n_dmd_atom_t; i++ ){
    for( int j=i; j<=n_dmd_atom_t; j++ ){
      type1 = static_cast<dmd_atom_t>( i ); stype1=dmd_atom_t2string( type1 ) ;
      type2 = static_cast<dmd_atom_t>( j ); stype2=dmd_atom_t2string( type2 ) ;
      /*cout<<"stype1=|"<<stype1<<"| stype2=|"<<stype2<<"|\n";*/
      get_default_Van_der_Waals( stype1, stype2, nonel_col[i][j], DAT ) ;
      get_default_hydrophobicity( stype1, stype2, nonel_col[i][j], DAT ) ;
      get_default_electrostatics( stype1, stype2, nonel_col[i][j], DAT ) ;
      get_signal_Hydrogen_Bond( stype1, stype2, nonel_col[i][j], DAT ) ;
      string2char( nonel_col[i][j], cbuf2 ) ;
      sprintf( cbuf, "%3d %3d %s\n", i, j, cbuf2 ) ;
      DMD << cbuf ;
    }
  }
}
/*=====================================================*/
void  printLINKED_PAIRS( ifstream &DAT, ofstream &DMD ){
  cout << "LINKED_PAIRS ...\n" ;
  DMD << txtKeyWords[LINK_PAIRS] << endl;
  int i, j, pos = DAT.tellg( ) ; DAT.seekg( 0 ) ;
  double HBstrength = get_HB_strength( DAT ) ;
  string link, aa_at1, aa_at2 ;
  char cbuf[256], cbuf2[256] ;
  dmd_atom_t type1, type2 ;
  /*output hydrogen bond links*/
  goto_dat_key( DAT, LINK_HB_DAT ) ; getline( DAT, link ) ;
  while( !is_end( link, LINK_HB_DAT ) ){
    /*cout<<"link=|"<<link<<"|\n";*/
    if( !is_comment(link) ){
      rescale_barriers_of_link_by_factor( link, HBstrength ) ;
    }
    if( out_linked_pair( link, aa_at1, aa_at2, cbuf ) ){
      type1 = string2dmd_atom_t( aa_at1 ) ;  i = static_cast<int>( type1 ) ;
      type2 = string2dmd_atom_t( aa_at2 ) ;  j = static_cast<int>( type2 ) ;
      if( i==0 ){ cout <<"|"<<aa_at1<<"| has no type!\n"; exit(1); }
      if( j==0 ){ cout <<"|"<<aa_at2<<"| has no type!\n"; exit(1); }
      sprintf( cbuf2, "%3d %3d  %s\n", i, j, cbuf ) ;
      DMD << cbuf2 ;
    }
    getline( DAT, link ) ;
  }
  DAT.seekg( 0 ) ;
  /*output permanent links*/
  goto_dat_key( DAT, LINK_PAIRS_DAT ) ; getline( DAT, link ) ;
  while( !is_end( link, LINK_PAIRS_DAT ) ){
    if( out_linked_pair( link, aa_at1, aa_at2, cbuf ) ){
      type1 = string2dmd_atom_t( aa_at1 ) ;  i = static_cast<int>( type1 ) ;
      type2 = string2dmd_atom_t( aa_at2 ) ;  j = static_cast<int>( type2 ) ;
      if( i==0 ){ cout <<"|"<<aa_at1<<"| has no type!\n"; exit(1); }
      if( j==0 ){ cout <<"|"<<aa_at2<<"| has no type!\n"; exit(1); }
      if( have_rotamers( DAT, aa_at1, aa_at2 ) ){
	/*the default_rotamer_barrier function works only if "link" is only
	  a torsion bond, ie, no mixture of angle and torsion bonds.*/
	default_rotamer_barrier( DAT, link, cbuf ) ;
      }
      sprintf( cbuf2, "%3d %3d  %s\n", i, j, cbuf ) ;
      DMD << cbuf2 ;
    }
    getline( DAT, link ) ; /*cout<<"link=|"<<link<<"|\n" ;*/
  }
  DAT.seekg( pos ) ;
}
/*=====================================================*/
void printREACT( ifstream &DAT, ofstream &DMD ){
  cout << "REACT ...\n" ;
  string line, sold1, sold2, snew1, snew2 ;
  dmd_atom_t told1, told2, tnew1, tnew2 ;
  char cbuf[ 128 ] ;
  DMD << txtKeyWords[REACT] << endl;
  goto_dat_key( DAT, REACT_DAT ) ; getline( DAT, line ) ;
  while( !is_end( line, REACT_DAT ) ){
    if( get_reactants_products( line, sold1, sold2, snew1, snew2 ) ){
      told1 = string2dmd_atom_t( sold1 ) ;
      told2 = string2dmd_atom_t( sold2 ) ;
      tnew1 = string2dmd_atom_t( snew1 ) ;
      tnew2 = string2dmd_atom_t( snew2 ) ;
      sprintf( cbuf, "%3d %3d %3d %3d    1\n", 
	       static_cast<int>( told1 ), static_cast<int>( told2 ), 
	       static_cast<int>( tnew1 ), static_cast<int>( tnew2 )  ) ;
      DMD << cbuf ;
    }
    getline( DAT, line ) ;
  }
}
/*=====================================================
In the relaxed conformations, there are no hydrogen bonds formed, thus all bonds 
are permanent bonds                                     */
void printLIST_ATOMS_LIST_OF_BONDS( ifstream &RLX, ofstream &DMD ){
  cout << "LIST_ATOMS ...\n" ;
  cout << "LIST_BONDS ...\n" ;
  int pos = RLX.tellg( ) ; RLX.seekg( 0 ) ;
  string line ;
  goto_txt_key( RLX, LIST_ATOMS ) ;
  DMD << txtKeyWords[LIST_ATOMS] << endl ;
  while( getline( RLX, line ) ){ DMD << line << endl ; }
  RLX.clear( ios::goodbit ) ;  RLX.seekg( pos ) ;
}
/*=====================================================*/
void  printLIST_PERM_BONDS( ifstream &RLX, PDBchain& chain, ofstream &DMD ){
  cout << "LIST_PERM_BONDS ...\n" ;
  int pos = RLX.tellg( ) ; RLX.seekg( 0 ) ;
  string line, aa_at1, aa_at2 ;
  dmd_atom_t type1, type2 ;
  int *ids = NULL, n=0 ;
  char cbuf[32] ;
  PDBatom at1, at2 ;
  /*no list of permanent bonds in the  prerelaxed conformation, thus go to
   list of bonds*/
  goto_txt_key( RLX, LIST_BONDS ) ;
  DMD << txtKeyWords[LIST_PERM_BONDS] << endl ;
  while( getline( RLX, line ) ){
    /*DMD <<"line="<<line<<"\n" ;*/
    if(line.find("//") == string::npos){/*line is not a comment*/
      ids = split_line2ints( line, n ); /*cout<<"line=|"<<line<<"|, n="<<n<<endl;*/
      if( n==2 && ids ){
	at1=chain.getAtomWithIndex(ids[0]) ; at2=chain.getAtomWithIndex(ids[1]);
	aa_at1 = at1.get_resName_nameII( );  aa_at2 = at2.get_resName_nameII( );
	type1 = string2dmd_atom_t( aa_at1 );  type2 = string2dmd_atom_t( aa_at2 );
	sprintf( cbuf, "%4d %4d %3d %3d\n",ids[0], ids[1], 
		 static_cast<int>(type1), static_cast<int>(type2) ) ;
	DMD << cbuf ;
	delete [] ids ;  ids = NULL ; n = 0 ;
      }
    }
  }
  RLX.clear( ios::goodbit ) ;  RLX.seekg( pos ) ;  
}
/*=====================================================*/
void printHBA_LIST( ifstream &DAT, PDBchain& chain, ofstream &DMD ){
  cout << "HB_LIST ...\n" ;
  string line, *associates=NULL, *atoms=NULL, sbuf, sbuf2 ;
  int *rel_index=NULL, n=0, l = chain.length( ), ix ;
  bool flag ;
  char cbuf[64] ;
  DMD << txtKeyWords[HB_LIST] << endl ;
  goto_dat_key( DAT, HB_LIST_DAT ) ;
  getline( DAT, line ) ;
  while( !is_end( line, HB_LIST_DAT ) ){
    if( get_associates( line, associates, rel_index, atoms, n ) ){

      for( int i=1; i<=l; i++ ){/*go through all amino acids in chain*/
    	flag = true ;
	for( int j=0; j<n; j++ ){/*check presence of all associates*/
	  if(!chain.is_there_resName_name_at_resIndexII( associates[j], 
						       i+rel_index[j]) ){
	    flag = false ; break ;
	  }
	}
	if( flag ){/*all associates present. Retrieve their atom indexes*/
	  sbuf.assign("") ;
	  ix=chain.getSerialNumberOfAtomNameAtIndex(atoms[0],i+rel_index[0]);
	  /*below, "-1" indicates amino acid index for DMD simulations is shifted*/
	  sprintf( cbuf, "%3d %4d ", ix, i+rel_index[0]-1 ) ;
	  sbuf2.assign( cbuf ) ; sbuf += sbuf2 ;
	  for( int j=1; j<n; j++ ){/*get atom serial of all associates*/
	    ix=chain.getSerialNumberOfAtomNameAtIndex(atoms[j],i+rel_index[j]);
	    sprintf( cbuf, "%3d ", ix ) ; sbuf2.assign( cbuf ) ; sbuf += sbuf2 ;
	  }
	  DMD << sbuf << endl ;
	}
      }
    }
    getline( DAT, line ) ;
  }
}
/*=====================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 4 ){
    system("clear");
    cout << "Usage: ./txt2txt.x relaxed.txt model.pdb force_field init.txt\n" ;
    cout << "relaxed.txt: DMD conformation with no hydrogen bond formed and with all permanent bonds present (outcome of a relaxation procedure).\n";
    cout << "model.pdb: any conformation of the protein in PDB format\n";
    cout << "force_field: file with info on atoms (mass, radius, bonded and nonbonded interactions)\n";
    cout << "init.txt: output conformation where we include hydrogen bonds in the backbone and physical interactions.\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(  int argc, char ** argv ){

  /*check number of arguments*/
  if( test_input(argc, argv) == false ) { return 1 ; }

  /*retrieve input*/
  int n_arg = argc -1 ;
  ifstream RLX( argv[ 1 ] ) ;
  ifstream PDB( argv[ 2 ] ) ;
  PDBchain chain ;   PDB >> chain ;   PDB.close( ) ;
  ifstream DAT( argv[ 3 ] ) ;
  ofstream DMD( argv[ 4 ] ) ;

  chain.remove_all_hydrogens( ) ;
  chain.renumberFully( 1 ) ;
  chain.renumberFullyAtomSerialNumber( 1 ) ;

  double *hcr_list = new double[ n_dmd_atom_t + 1] ;

  printSYS_SIZE( RLX, DMD ) ;
  printNUM_ATOMS( chain, DMD ) ;
  printATOM_TYPE( DAT, DMD, hcr_list ) ;
  printNONEL_COL( DAT, DMD, hcr_list ) ;
  printLINKED_PAIRS( DAT, DMD ) ;
  printREACT( DAT, DMD ) ;
  printLIST_ATOMS_LIST_OF_BONDS( RLX, DMD ) ;
  printLIST_PERM_BONDS( RLX, chain, DMD ) ;
  printHBA_LIST( DAT, chain, DMD ) ;

  return 0 ;
}
