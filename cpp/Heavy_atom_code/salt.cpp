#include<iostream>
#include<fstream>
using namespace std;
#include "salt.h"

/*=================================================================*/
bool is_SALT( string &theString ){
  string query=theString.substr(17,3) ;
  if( all_salts.find(query) == string::npos ) return false ;
  return true ;
}
/*=================================================================*/
salt_XXX string_to_salt_XXX( const string &salt_name_XXX  ){
  int n;
  if(salt_name_XXX == "ACT" ){ n=0;  return static_cast<salt_XXX>(n) ; }
}
/*=================================================================*/
string dmd_salt_t2string( dmd_salt_t &type ){
  switch( type ){
  case _ACT_CA_:     return char2string("ACT  CA " );
  case _ACT_CB_:     return char2string("ACT  CB " );
  case _ACT_OG1_:    return char2string("ACT  OG1" );
  case _ACT_BOG1_:   return char2string("ACT BOG1" );
  case _ACT_OG2_:    return char2string("ACT  OG2" );
  case _ACT_BOG2_:   return char2string("ACT BOG2" );
  default :          return char2string("ERROR"   ) ;
  }
}
/*=================================================================*/
dmd_salt_t string2dmd_salt_t( string &salt_and_atom ){
  if( salt_and_atom == "ACT  CA " ){ return  _ACT_CA_      ; }    
  if( salt_and_atom == "ACT  CB " ){ return  _ACT_CB_      ; }    
  if( salt_and_atom == "ACT  OG1" ){ return  _ACT_OG1_     ; }    
  if( salt_and_atom == "ACT BOG1" ){ return  _ACT_BOG1_    ; }    
  if( salt_and_atom == "ACT  OG2" ){ return  _ACT_OG2_     ; }    
  if( salt_and_atom == "ACT BOG2" ){ return  _ACT_BOG2_    ; }    
}
/*=================================================================*/
bool belong_to_same_XXX( dmd_salt_t &type1, dmd_salt_t &type2 ){
  string stype1 = dmd_salt_t2string( type1 );
  string stype2 = dmd_salt_t2string( type2 );
  if( stype1.substr( 0, 3 ) == stype2.substr( 0, 3 ) ){
    return true ;
  }
  return false ;
}
/*=================================================================*/
string dmd_salt_t2salt( dmd_salt_t &type ){
  string stype = dmd_salt_t2string( type );
  return stype.substr( 4, 4 ) ;
}
/*=================================================================*/
string dmd_salt_t2at( dmd_salt_t &type ){
  string stype = dmd_salt_t2string( type );
  return stype.substr( 4, 4 ) ;
}
/*=================================================================*/
/*=================================================================*/
/*=================================================================*/

/*------------ class PDBsalt  (PDBsalt : public listPDBatom)---*/

ifstream &operator>>(ifstream &PDB, PDBsalt &theSalt){
  int position ;
  string line ; //string junkLine ;
  PDBatom atomEntry ;  getline( PDB, line ) ;  theSalt.empty( ) ;
  while( !(is_ATOM(line)) || !(is_SALT(line)) ){/*find first ATOM entry*/
    //cout <<"EXIT MESSAGE: We are reading a non-ATOM line\n" ;
    if( PDB.eof() ) return PDB ; 
    if( is_TER(line) || is_END(line) ){ PDB.seekg( position ); return PDB ; }
    position = PDB.tellg( ) ; getline( PDB, line );
  }

  line >> atomEntry ; 
  theSalt.name = atomEntry.getResName( ) ;
  theSalt.saltSeq = atomEntry.getResSeq( ) ; //initialize
  while( atomEntry.getResSeq( ) == theSalt.saltSeq || is_TER(line) ){
    /*cout<<"atomEntry=\n"<<atomEntry ;
    cout<<"atomEntry.getResSeq( )="<<atomEntry.getResSeq( )<<endl;
    cout<<"theSalt.saltSeq="<<theSalt.saltSeq<<endl;*/
    if( is_ATOM(line) ){ 
      theSalt.insertAtBack( atomEntry ) ;  /*cout<<"is ATOM line\n" ;*/
    }
    position = PDB.tellg( ) ; getline( PDB, line ) ;
    if( is_ANISOU(line) || is_SIGATM(line) || is_HETATM(line) ){ 
      /*cout<<"we continue\n";*/
      continue ; 
    }
    if( PDB.eof( ) ){ 
      /*cout<<"is end of file\n";*/
      return PDB ; 
    }
    if( !(is_ATOM(line)) ){ 
      /*cout<<"is not atom line\n";*/
      PDB.seekg( position ); return PDB ; 
    }
    line >> atomEntry ;
  }

  PDB.seekg( position ); return PDB ;
}
/*=================================================================*/
PDBsalt::PDBsalt( ){
  firstPtr = lastPtr = NULL ; name=""; saltSeq=-1;
}
/*=================================================================*/
PDBsalt::PDBsalt( const PDBsalt &clone ){
  /*cout<<"PDBsalt copy constructor\n";*/
  nodePDBatom *ptr = clone.firstPtr ;
  while( ptr ){
    /*cout<<ptr->theAtom;*/
    this->insertAtBack( ptr->theAtom ) ;
    ptr = ptr->next ;
  }
  name = clone.name ;
  saltSeq = clone.saltSeq ;
  /*cout<<"finished PDBsalt copyconstructor\n";*/  
}
/*=================================================================*/
PDBsalt::PDBsalt( string &salt_name ){
  PDBatom sd ;
  for(int i=0; i<Nsalt_atoms; i++ ){
    if( salt_atoms[i].find( salt_name ) != string::npos ){      
      salt_atoms[i] >> sd ;
      this->insertAtBack( sd ) ;
    }
  }
}
/*=================================================================*/
PDBsalt &PDBsalt::operator=( PDBsalt &right  ){
  if( ! this->isEmpty() ){  this->empty( ) ; }
  if( ( &right != this )  && !( right.isEmpty() ) ) {
    firstPtr = getNewNode( right.firstPtr->theAtom ) ;
    nodePDBatom *currentRightPtr=right.firstPtr->next, *currentPtr=firstPtr ;
    while( currentRightPtr != NULL ){
      currentPtr->next = getNewNode( currentRightPtr->theAtom ) ;
      currentPtr = currentPtr->next ;
      currentRightPtr = currentRightPtr->next ;
    }
    lastPtr = currentPtr ;
  }
  name = right.name ;  saltSeq = right.saltSeq ;
  return *this ;
}
/*=================================================================*/
PDBsalt &PDBsalt::operator=( const PDBsalt &right  ){
  PDBsalt tmp( right ) ;
  *this = tmp ; /*previous operator assignment definition*/
  return *this ;
}
/*=================================================================*/
/*=================================================================*/
