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

static const int n_bkb_t = 7 ; 

typedef enum {
  _UNIDENTIFIED_CODE=0,
  _FILL_,
  _EPSILON_
}error_code ;

struct go_atom{
  PDBatom atom ;
  string aa_at ;
  int id ;
  dmd_atom_t type ;
  int hv_type ;
  int go_type ;
  bool has_go ;
} ;

/*=====================================================*/
void print_go_atom( const go_atom &the_atom ){
  cout<<the_atom.atom<<endl<<"aa_at=|"<<the_atom.aa_at<<"|, id="<<the_atom.id
      <<", type="<<static_cast<int>(the_atom.type)
      <<", hv_type="<<the_atom.hv_type<<", go_type="<<the_atom.go_type<<endl ;
}
/*=====================================================*/
void print_chain_go( go_atom *chain_go, const int &n_atoms ){
  for( int i=1; i<=n_atoms; i++ ){  print_go_atom( chain_go[i] );  }
}
/*=====================================================*/
bool exit_error( const error_code &the_error){
  switch( the_error ){
  case _FILL_: cerr<<"\nERROR #"<<static_cast<int>( _FILL_ )
		   <<" : number of filled go atoms different than n_atoms!\n" ;
  case _EPSILON_:cerr<<"\nERROR #"<<static_cast<int>( _EPSILON_ )
		     <<" : energies must be negative!\n" ;
  default : cerr<<"\nunidentified error!\n" ;
  }
  exit(1);
}
/*=====================================================*/
int out_n_go_types( PDBchain &chain, const int &n_atoms ){
  int n_go_types = n_bkb_t ;
  PDBatom atom ;
  for( int i=1; i<=n_atoms; i++ ){
    atom = chain.getAtomWithIndex( i ) ;
    if( !atom.is_backbone( ) ){ n_go_types++ ; }
  }
  return n_go_types ;
}
/*=====================================================*/
bool fill_arrays( go_atom *chain_go, const int &n_atoms, int *go2hv,
		  dmd_atom_t *go2dmd_atom_t, int *go2id, int &n_go_types,
		  PDBchain &chain ){
  n_go_types = n_bkb_t ;
  for( int i=1; i<=n_atoms; i++ ){
    chain_go[i].atom = chain.getAtomWithIndex( i ) ;
    chain_go[i].aa_at = chain_go[i].atom.get_resName_nameII( ) ;
    chain_go[i].id = i ;
    chain_go[i].type = string2dmd_atom_t( chain_go[i].aa_at ) ;
    chain_go[i].hv_type = static_cast<int>( chain_go[i].type ) ;
    if( !chain_go[i].atom.is_backbone( ) ){
      chain_go[i].go_type = ++n_go_types ;
      chain_go[i].has_go = true ;
    }
    else{ /*for backbone atoms, heavy atom type same as Go type*/
      chain_go[i].go_type = chain_go[i].hv_type ;
      chain_go[i].has_go = false ; 
    }
    go2hv[ chain_go[i].go_type ] = chain_go[i].hv_type ;
    go2dmd_atom_t[ chain_go[i].go_type ] = chain_go[i].type ;
    go2id[ chain_go[i].go_type ] = chain_go[i].id ;
  }
  for( int i=1; i<=n_bkb_t; i++){/*fill all blanks in go2hv and go2dmd_atom_t*/
    go2dmd_atom_t[ i ] = static_cast<dmd_atom_t>( i ) ;
    go2hv[ i ] = i ;
    go2id[ i ] = 0 ;
  }
  return true ;
}
/*=====================================================*/
void printSYS_SIZE( ifstream &RLX, ofstream &DMD ){
  cout << "SYS_SIZE ...\n" ;
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
}
/*=====================================================*/
void  printNUM_ATOMS( const int &n_atoms, ofstream &DMD ){
  cout << "NUM_ATOMS ...\n" ;
  DMD << txtKeyWords[NUM_ATOMS] << endl ;
  DMD << n_atoms << endl ;
}
/*=====================================================*/
void printATOM_TYPE( ifstream &DAT, ostream& DMD, const int &n_go_types,
		     dmd_atom_t *go2dmd_atom_t, double *hcr_list ){
  cout << "ATOM_TYPE ...\n" ;
  DMD << txtKeyWords[TYPE_ATOMS] << endl;
  dmd_atom_t type ;
  string amino_and_atom ;
  char cbuf[300] ;
  for( int i=1; i<=n_go_types; i++){
    sprintf( cbuf, "%3d ", i ) ; 
    DMD << cbuf ;
    type = go2dmd_atom_t[ i ] ; 
    /*cout<<i<<" "<<static_cast<int>(go2dmd_atom_t[ i ])<<endl;*/
    amino_and_atom = dmd_atom_t2string( type ) ;/*amino_param.cpp*/
    print_atom_type_line( cbuf, DAT, amino_and_atom ) ;/*atom_param.cpp*/
    hcr_list[ i ] = out_hard_core_radius( DAT, amino_and_atom ) ;
    DMD << cbuf ;
  }
}
/*=====================================================*/
void get_Go_interaction( const int &idi, const int &idj, const double &go_co_f,
			 const double &epsilon, PDBchain &chain,
			 string &nonel_col, go_atom *chain_go,
			 double *hcr_list, ifstream &DAT ){
  PDBatom atomi, atomj ;
  double mind, maxd, avd ;
  static const double f = get_inner_VW_factor( DAT ) ;
  static const int min_sep = 2 ;
  char cbuf[256] ;
  int resIi, resIj ;
  if( idi*idj > 0 ){/*idi=0 when atom is backbone. Same for idj*/
    atomi = chain.getAtomWithIndex( idi ) ; resIi = atomi.getResSeq( ) ;
    atomj = chain.getAtomWithIndex( idj ) ; resIj = atomj.getResSeq( ) ; 
    avd = hcr_list[chain_go[idi].go_type] + hcr_list[chain_go[idj].go_type] ;
    mind = f * avd ; /*(1-go_co_f) * avd ;*/
    maxd = (1+go_co_f)* avd ;
    if( (atomi.d2( atomj )<maxd*maxd) ){
      if( (abs(resIj-resIi)<min_sep) ){/*neighbor residues do not interact*/
	sprintf( cbuf, "%8.6f %8.6f 0.000000", mind, maxd) ;
      }
      else{ sprintf( cbuf, "%8.6f %8.6f %8.6f", mind, maxd, epsilon ) ; }
    }
    else{
      sprintf( cbuf, "%8.6f %8.6f %8.6f", mind, maxd, -epsilon ) ;      
    }
    nonel_col.assign( cbuf ) ;
  }
}
/*=====================================================*/
void printNONEL_COL( ifstream &DAT, ostream& DMD, const int &n_go_types,
		     dmd_atom_t *go2dmd_atom_t, int *go2id, PDBchain &chain,
		     const double &go_co_f, const double &epsilon, 
		     go_atom *chain_go, double *hcr_list ){
  cout << "NONEL_COL ...\n" ;
  DMD << txtKeyWords[NONEL_COL] << endl;
  string nonel_col[n_go_types+1][n_go_types+1] ;
  dmd_atom_t type1, type2 ;
  string stype1, stype2 ;
  int idi, idj ;
  char cbuf[256], cbuf2[256] ;
  for( int i=1; i<=n_go_types; i++ ){
    for( int j=i; j<=n_go_types; j++ ){  /*cout<<"i="<<i<<", j="<<j<<endl;*/
      idi = go2id[ i ] ;  idj = go2id[ j ] ;
      type1 = go2dmd_atom_t[ i ]; stype1=dmd_atom_t2string( type1 ) ;
      type2 = go2dmd_atom_t[ j ]; stype2=dmd_atom_t2string( type2 ) ;
      /*cout<<"stype1=|"<<stype1<<"| stype2=|"<<stype2<<"|\n";*/
      get_default_Van_der_Waals( stype1, stype2, nonel_col[i][j], DAT ) ;
      get_Go_interaction( idi, idj, go_co_f, epsilon, chain,
			  nonel_col[i][j], chain_go, hcr_list, DAT ) ;
      get_signal_Hydrogen_Bond( stype1, stype2, nonel_col[i][j], DAT ) ;
      string2char( nonel_col[i][j], cbuf2 ) ;
      sprintf( cbuf, "%3d %3d %s\n", i, j, cbuf2 ) ;
      DMD << cbuf ;
    }
  }
}
/*=====================================================*/
void printLINKED_PAIRS( ifstream &DAT, ofstream &DMD, const int &n_go_types,
			dmd_atom_t *go2dmd_atom_t, const double &factor,
			const double &fr){
  cout << "LINKED_PAIRS ...\n" ;
  DMD << txtKeyWords[LINK_PAIRS] << endl;
  dmd_atom_t type1, type2 ;
  string stype1, stype2, link ;
  char cbuf[256], cbuf2[256] ;
  double factor2=get_default_rotamer_barrier(DAT); factor*=fr;
  
  for( int i=1; i<=n_go_types; i++ ){
    for( int j=i; j<=n_go_types; j++ ){  /*cout<<"i="<<i<<", j="<<j<<endl;*/
      type1 = go2dmd_atom_t[ i ]; stype1=dmd_atom_t2string( type1 ) ;
      type2 = go2dmd_atom_t[ j ]; stype2=dmd_atom_t2string( type2 ) ;
      if( is_there_linked_pair( DAT, stype1, stype2, link ) ){
	/*cout<<"stype1=|"<<stype1<<"|, stype2=|"<<stype2<<"|\n"<<link<<endl;*/
	if( is_HB_link( link, DAT ) ){
	  /*cout <<"they are HB_link\n";*/
	  rescale_barriers_of_link_by_factor( link, factor ) ;	
	}
	else{
	  rescale_barriers_of_link_by_factor( link, factor2 ) ;
	}
	out_linked_pair( link, stype1, stype2, cbuf ) ;
	sprintf( cbuf2, "%3d %3d  %s\n", i, j, cbuf ) ;
	DMD << cbuf2 ;
      }
    }
  }
}
/*=====================================================*/
void printREACT( ifstream &DAT, ofstream &DMD, const int &n_go_types,
		 dmd_atom_t *go2dmd_atom_t ){
  cout << "REACT ...\n" ;
  DMD << txtKeyWords[REACT] << endl;
  dmd_atom_t type1, type2, new_type1, new_type2 ;
  string stype1, stype2, reaction, snew1, snew2 ;
  char cbuf[64] ;
  for( int i=1; i<=n_go_types; i++ ){
    for( int j=i; j<=n_go_types; j++ ){  /*cout<<"i="<<i<<", j="<<j<<endl;*/
      type1 = go2dmd_atom_t[ i ]; stype1=dmd_atom_t2string( type1 ) ;
      type2 = go2dmd_atom_t[ j ]; stype2=dmd_atom_t2string( type2 ) ;
      if( do_they_react( DAT, stype1, stype2, reaction ) ){
	get_reactants_products(  reaction, stype1, stype2, snew1, snew2 ) ;
	/*we assume only backbone atoms react, thus heavy atom type 
	  same as go_type*/
	type1 = string2dmd_atom_t( stype1 ) ;
	type2 = string2dmd_atom_t( stype2 ) ;
	new_type1 = string2dmd_atom_t( snew1 ) ;
	new_type2 = string2dmd_atom_t( snew2 ) ;
	sprintf( cbuf, "%3d %3d %3d %3d\n",
		 static_cast<int>(type1), static_cast<int>(type2),
		 static_cast<int>(new_type1), static_cast<int>(new_type2) ) ;
	DMD<<cbuf;
      }
    }
  }
}
/*=====================================================*/
void  printLIST_ATOMS_LIST_OF_BONDS( ifstream &RLX, ofstream &DMD,
				     go_atom *chain_go ){
  cout << "LIST_ATOMS ...\n" ;
  cout << "LIST_BONDS ...\n" ;
  int a, b, pos = RLX.tellg( ), n_w ; RLX.seekg( 0 ) ;
  double *w = NULL ;
  string line ;
  char cbuf[256] ;
  goto_txt_key( RLX, LIST_ATOMS ) ;
  DMD << txtKeyWords[LIST_ATOMS] << endl ;
  while( getline( RLX, line ) ){
    if( line.find("LIST OF BONDS") != string::npos ){ break ; }
    if( line.find("//") != string::npos ){ DMD << line << endl ; continue ; }
    w = split_line2doubles( line, n_w ) ; /*cout<<"n_w="<<n_w<<endl;*/
    a = static_cast<int>( w[0] ) ;
    b = chain_go[ a ].go_type ;
    /*cout<<"a="<<a<<", b="<<b<<", w[7]="<<w[7]<<endl;*/
    sprintf( cbuf, "%4d %4d %18.13f %18.13f %18.13f %18.13f %18.13f %18.13f\n",
	     a, b, w[2], w[3], w[4], w[5], w[6], w[7] ) ;
    delete [ ] w ;
    DMD << cbuf ; 
  }
  DMD << line << endl ;
  while( getline( RLX, line ) ){ DMD << line << endl ; }
  RLX.clear( ios::goodbit ) ;  RLX.seekg( pos ) ;
}
/*=====================================================*/
void printLIST_PERM_BONDS( ifstream &RLX, go_atom *chain_go, ofstream &DMD ){
  cout << "LIST_PERM_BONDS ...\n" ;
  DMD << txtKeyWords[LIST_PERM_BONDS] << endl ;
  goto_txt_key( RLX, LIST_BONDS ) ;
  string line ;
  int *ids = NULL, n = 0 ;
  char cbuf[32] ;
  while( getline( RLX, line ) ){
    if(line.find("//") == string::npos){/*line is not a comment*/
      ids = split_line2ints( line, n );
      sprintf( cbuf, "%4d %4d %3d %3d\n",ids[0], ids[1], 
	       chain_go[ ids[ 0 ] ].go_type, chain_go[ ids[ 1 ] ].go_type ) ;
      DMD << cbuf ;
      delete [] ids ;  ids = NULL ; n = 0 ;
    }
  }
}
/*=====================================================*/
void printHBA_LIST( ifstream &DAT, PDBchain chain, const int &n_atoms, 
		    ofstream &DMD ){
  cout << "HB_LIST ...\n" ;
  DMD << txtKeyWords[HB_LIST] << endl ;
  goto_dat_key( DAT, HB_LIST_DAT ) ;
  string line, *associates=NULL, *atoms=NULL, sbuf, sbuf2 ;
  int *rel_index=NULL, n=0, ix ;
  bool flag ;
  char cbuf[64] ;
  getline( DAT, line ) ;
  while( !is_end( line, HB_LIST_DAT ) ){
    if( get_associates( line, associates, rel_index, atoms, n ) ){
      for( int i=1; i<=n_atoms; i++ ){/*go through all amino acids in chain*/
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
	    sprintf( cbuf, "%3d ", ix ) ; 
	    sbuf2.assign( cbuf ) ; sbuf += sbuf2 ;
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
  if( number_of_arguments != 8 ){
    system("clear");
    cout << "Usage: ./txt2txt_Go.x rlxd.txt nat.pdb frc_fld init.txt f_go e_go e_HB fr\n" ;
    cout << "rlxd.txt: DMD conformation with no hydrogen bond formed and with all\n" ;
    cout << "          permanent bonds present (outcome of a relaxation procedure).\n";
    cout << "nat.pdb : native conformation of the protein in PDB format\n";
    cout << "frc_fld : force field file with info on atoms (mass, radius, bonded\n" ;
    cout << "          and nonbonded interactions)\n";
    cout << "init.txt: output conformation where we include hydrogen bonds in the\n" ;
    cout << "f_go    : (1+go_co_f)*(VdW_i+VdW_j) will be Go range for i,j\n" ;
    cout << "e_go    : interaction strenght for Go model (must be negative number!)\n";
    cout << "e_HB    : interaction strenght for hydrogen bond (must be negative number!)\n";
    cout << " fr       : energy-rotamer factor (Er=Er(frc_fld)*f_rot).Default fr=0\n";
    cout << "          backbone and physical interactions.\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(  int argc, char ** argv ){

  /*check number of arguments*/
  if( test_input(argc, argv) == false ) { return 1 ; }

  /*retrieve input*/
  int n_args = argc -1 ;
  ifstream RLX( argv[ 1] ) ;
  ifstream PDB( argv[ 2 ] ) ;
  PDBchain chain ;   PDB >> chain ;  PDB.close( ) ;
  ifstream DAT( argv[ 3 ] ) ;
  ofstream DMD( argv[ 4 ] ) ;
  double go_co_f = 0.35;
  double e_go = -1;
  double factor = -0.1;
  double fr=0.0;
  if(n_args>4){
    go_co_f = atof( argv[ 5 ] ) ;
    e_go = atof( argv[ 6 ] );
    factor = atof( argv[ 7 ] ) ;
    fr = atof( argv[ 8 ] ) ;
    cout<<"Using customized parameters:\n";
    printf("f_go=%4f e_go=%4f e_HB=%4f fr=%4f\n",
	   go_co_f,e_go,factor,fr);
  }
  else{
    cout<<"Using default parameters:\n";
    cout<<" f_go=0.35 e_go=-1 e_HB=-0.1 fr=0\n";
  }
  if( (e_go>0) || (factor>0) ){ exit_error( _EPSILON_ ) ; }
  factor /= -1.0 ;

  chain.remove_all_hydrogens( ) ;
  chain.renumberFully( 1 ) ;
  chain.renumberFullyAtomSerialNumber( 1 ) ;
  
  int n_atoms = chain.numberOfAtoms( ) ;
  go_atom *chain_go = new go_atom[ n_atoms + 1 ] ;
  int n_go_types = out_n_go_types( chain,n_atoms ) ;
  int *go2hv  = new int[ n_go_types + 1 ] ;
  dmd_atom_t *go2dmd_atom_t = new dmd_atom_t[ n_go_types + 1 ] ;
  int *go2id  = new int[ n_go_types + 1 ] ;
  fill_arrays( chain_go, n_atoms, go2hv, go2dmd_atom_t, 
	       go2id, n_go_types, chain ) ;

  /*print_chain_go( chain_go, n_atoms ) ;*/

  double *hcr_list = new double[ n_go_types + 1] ;
  printSYS_SIZE( RLX, DMD ) ;
  printNUM_ATOMS( n_atoms, DMD ) ;
  printATOM_TYPE( DAT, DMD, n_go_types, go2dmd_atom_t, hcr_list ) ;
  printNONEL_COL( DAT, DMD, n_go_types, go2dmd_atom_t, go2id,
    chain, go_co_f, e_go, chain_go, hcr_list ) ;
  printLINKED_PAIRS( DAT, DMD, n_go_types, go2dmd_atom_t, factor,fr);
  printREACT( DAT, DMD, n_go_types, go2dmd_atom_t ) ;
  printLIST_ATOMS_LIST_OF_BONDS( RLX, DMD, chain_go ) ;
  printLIST_PERM_BONDS( RLX, chain_go, DMD ) ;
  printHBA_LIST( DAT, chain,  n_atoms, DMD ) ;


  return 0 ;
}
