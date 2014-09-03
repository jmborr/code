using namespace std;
#include<iostream>
#include<fstream>
#include<string>
#include "amino_acid_param.h"
#include "miscellanea.h"

aa_XXX string_to_aa_XXX( const string &aa_name_XXX ){
   int n;
   if(aa_name_XXX == "GLY" ){ n=0;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "ALA" ){ n=1;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "VAL" ){ n=2;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "LEU" ){ n=3;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "ILE" ){ n=4;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "SER" ){ n=5;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "THR" ){ n=6;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "CYS" ){ n=7;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "MET" ){ n=8;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "PRO" ){ n=9;  return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "ASP" ){ n=10; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "ASN" ){ n=11; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "GLU" ){ n=12; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "GLN" ){ n=13; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "LYS" ){ n=14; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "ARG" ){ n=15; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "HIS" ){ n=16; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "PHE" ){ n=17; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "TYR" ){ n=18; return static_cast<aa_XXX>(n) ; }
   if(aa_name_XXX == "TRP" ){ n=19; return static_cast<aa_XXX>(n) ; }
   n=-1; return static_cast<aa_XXX>(n) ;
}
/*=========================================================*/
string dmd_atom_t2string( dmd_atom_t &type ){
  switch( type ){
  case _BKB_N_:       return char2string("BKB  N  " );
  case _BKB_N_HB_:    return char2string("BKB BN  ") ;
  case _PRO_N_:       return char2string("PRO  N  ") ;
  case _BKB_CA_:      return char2string("BKB  CA ") ;
  case _BKB_C_:       return char2string("BKB  C  ") ;
  case _BKB_O_:       return char2string("BKB  O  ") ;
  case _BKB_O_HB_:    return char2string("BKB BO  ") ;
  case _ALA_CB_:      return char2string("ALA  CB ") ;
  case _ARG_CB_:      return char2string("ARG  CB ") ;
  case _ARG_CG_:      return char2string("ARG  CG ") ;
  case _ARG_CD_:      return char2string("ARG  CD ") ;
  case _ARG_NE_:      return char2string("ARG  NE ") ;
  case _ARG_CZ_:      return char2string("ARG  CZ ") ;
  case _ARG_NH1_:     return char2string("ARG  NH1") ;
  case _ARG_NH1_HB_:  return char2string("ARG BNH1") ;
  case _ARG_NH2_:     return char2string("ARG  NH2") ;
  case _ARG_NH2_HB_:  return char2string("ARG BNH2") ;
  case _ASN_CB_:      return char2string("ASN  CB ") ;
  case _ASN_CG_:      return char2string("ASN  CG ") ;
  case _ASN_OD1_:     return char2string("ASN  OD1") ;
  case _ASN_OD1_HB_:  return char2string("ASN BOD1") ;
  case _ASN_ND2_:     return char2string("ASN  ND2") ;
  case _ASN_ND2_HB_:  return char2string("ASN BND2") ;
  case _ASP_CB_:      return char2string("ASP  CB ") ;
  case _ASP_CG_:      return char2string("ASP  CG ") ;
  case _ASP_OD1_:     return char2string("ASP  OD1") ;
  case _ASP_OD1_HB_:  return char2string("ASP BOD1") ;
  case _ASP_OD2_:     return char2string("ASP  OD2") ;
  case _ASP_OD2_HB_:  return char2string("ASP BOD2") ;
  case _CYS_CB_:      return char2string("CYS  CB ") ;
  case _CYS_SG_:      return char2string("CYS  SG ") ;
  case _GLN_CB_:      return char2string("GLN  CB ") ;
  case _GLN_CG_:      return char2string("GLN  CG ") ;
  case _GLN_CD_:      return char2string("GLN  CD ") ;
  case _GLN_OE1_:     return char2string("GLN  OE1") ;
  case _GLN_OE1_HB_:  return char2string("GLN BOE1") ;
  case _GLN_NE2_:     return char2string("GLN  NE2") ;
  case _GLN_NE2_HB_:  return char2string("GLN BNE2") ;
  case _GLU_CB_:      return char2string("GLU  CB ") ;
  case _GLU_CG_:      return char2string("GLU  CG ") ;
  case _GLU_CD_:      return char2string("GLU  CD ") ;
  case _GLU_OE1_:     return char2string("GLU  OE1") ;
  case _GLU_OE1_HB_:  return char2string("GLU BOE1") ;
  case _GLU_OE2_:     return char2string("GLU  OE2") ;
  case _GLU_OE2_HB_:  return char2string("GLU BOE2") ;
  case _HIS_CB_:      return char2string("HIS  CB ") ;
  case _HIS_CG_:      return char2string("HIS  CG ") ;
  case _HIS_ND1_:     return char2string("HIS  ND1") ;
  case _HIS_ND1_HB_:  return char2string("HIS BND1") ;
  case _HIS_CD2_:     return char2string("HIS  CD2") ;
  case _HIS_CE1_:     return char2string("HIS  CE1") ;
  case _HIS_NE2_:     return char2string("HIS  NE2") ;
  case _HIS_NE2_HB_:  return char2string("HIS BNE2") ;
  case _ILE_CB_:      return char2string("ILE  CB ") ;
  case _ILE_CG1_:     return char2string("ILE  CG1") ;
  case _ILE_CG2_:     return char2string("ILE  CG2") ;
  case _ILE_CD1_:     return char2string("ILE  CD1") ;
  case _LEU_CB_:      return char2string("LEU  CB ") ;
  case _LEU_CG_:      return char2string("LEU  CG ") ;
  case _LEU_CD1_:     return char2string("LEU  CD1") ;
  case _LEU_CD2_:     return char2string("LEU  CD2") ;
  case _LYS_CB_:      return char2string("LYS  CB ") ;
  case _LYS_CG_:      return char2string("LYS  CG ") ;
  case _LYS_CD_:      return char2string("LYS  CD ") ;
  case _LYS_CE_:      return char2string("LYS  CE ") ;
  case _LYS_NZ_:      return char2string("LYS  NZ ") ;
  case _LYS_NZ_HB_:   return char2string("LYS BNZ ") ;
  case _MET_CB_:      return char2string("MET  CB ") ;
  case _MET_CG_:      return char2string("MET  CG ") ;
  case _MET_SD_:      return char2string("MET  SD ") ;
  case _MET_CE_:      return char2string("MET  CE ") ;
  case _PHE_CB_:      return char2string("PHE  CB ") ;
  case _PHE_CG_:      return char2string("PHE  CG ") ;
  case _PHE_CD1_:     return char2string("PHE  CD1") ;
  case _PHE_CD2_:     return char2string("PHE  CD2") ;
  case _PHE_CE1_:     return char2string("PHE  CE1") ;
  case _PHE_CE2_:     return char2string("PHE  CE2") ;
  case _PHE_CZ_:      return char2string("PHE  CZ ") ;
  case _PRO_CB_:      return char2string("PRO  CB ") ;
  case _PRO_CG_:      return char2string("PRO  CG ") ;
  case _PRO_CD_:      return char2string("PRO  CD ") ;
  case _SER_CB_:      return char2string("SER  CB ") ;
  case _SER_OG_:      return char2string("SER  OG ") ;
  case _SER_OG_HB_:   return char2string("SER BOG ") ;
  case _TRH_CB_:      return char2string("THR  CB ") ;
  case _TRH_OG1_:     return char2string("THR  OG1") ;
  case _TRH_OG1_HB_:  return char2string("THR BOG1") ;
  case _TRH_CG2_:     return char2string("THR  CG2") ;
  case _TRP_CB_:      return char2string("TRP  CB ") ;
  case _TRP_CG_:      return char2string("TRP  CG ") ;
  case _TRP_CD1_:     return char2string("TRP  CD1") ;
  case _TRP_CD2_:     return char2string("TRP  CD2") ;
  case _TRP_NE1_:     return char2string("TRP  NE1") ;
  case _TRP_CE2_:     return char2string("TRP  CE2") ;
  case _TRP_CE3_:     return char2string("TRP  CE3") ;
  case _TRP_CZ2_:     return char2string("TRP  CZ2") ;
  case _TRP_CZ3_:     return char2string("TRP  CZ3") ;
  case _TRP_CH2_:     return char2string("TRP  CH2") ;
  case _TYR_CB_:      return char2string("TYR  CB ") ;
  case _TYR_CG_:      return char2string("TYR  CG ") ;
  case _TYR_CD1_:     return char2string("TYR  CD1") ;
  case _TYR_CD2_:     return char2string("TYR  CD2") ;
  case _TYR_CE1_:     return char2string("TYR  CE1") ;
  case _TYR_CE2_:     return char2string("TYR  CE2") ;
  case _TYR_CZ_:      return char2string("TYR  CZ ") ;
  case _TYR_OH_:      return char2string("TYR  OH ") ;
  case _TYR_OH_HB_:   return char2string("TYR BOH ") ;
  case _VAL_CB_:      return char2string("VAL  CB ") ;
  case _VAL_CG1_:     return char2string("VAL  CG1") ;
  case _VAL_CG2_:     return char2string("VAL  CG2") ;
  default :           return char2string("ERROR"   ) ;
  }
}
/*============================================================*/
dmd_atom_t string2dmd_atom_t( string &amino_and_atom ){
  if( amino_and_atom == "BKB  N  " ){ return  _BKB_N_      ; }    
  if( amino_and_atom == "BKB BN  " ){ return  _BKB_N_HB_   ; } 
  if( amino_and_atom == "PRO  N  " ){ return  _PRO_N_      ; }    
  if( amino_and_atom == "BKB  CA " ){ return  _BKB_CA_     ; }   
  if( amino_and_atom == "BKB  C  " ){ return  _BKB_C_      ; }    
  if( amino_and_atom == "BKB  O  " ){ return  _BKB_O_      ; }    
  if( amino_and_atom == "BKB BO  " ){ return  _BKB_O_HB_   ; } 
  if( amino_and_atom == "ALA  CB " ){ return  _ALA_CB_     ; }   
  if( amino_and_atom == "ARG  CB " ){ return  _ARG_CB_     ; }   
  if( amino_and_atom == "ARG  CG " ){ return  _ARG_CG_     ; }   
  if( amino_and_atom == "ARG  CD " ){ return  _ARG_CD_     ; }   
  if( amino_and_atom == "ARG  NE " ){ return  _ARG_NE_     ; }   
  if( amino_and_atom == "ARG  CZ " ){ return  _ARG_CZ_     ; }   
  if( amino_and_atom == "ARG  NH1" ){ return  _ARG_NH1_    ; }  
  if( amino_and_atom == "ARG BNH1" ){ return  _ARG_NH1_HB_ ; }
  if( amino_and_atom == "ARG  NH2" ){ return  _ARG_NH2_    ; }  
  if( amino_and_atom == "ARG BNH2" ){ return  _ARG_NH2_HB_ ; }
  if( amino_and_atom == "ASN  CB " ){ return  _ASN_CB_     ; }   
  if( amino_and_atom == "ASN  CG " ){ return  _ASN_CG_     ; }   
  if( amino_and_atom == "ASN  OD1" ){ return  _ASN_OD1_    ; }  
  if( amino_and_atom == "ASN BOD1" ){ return  _ASN_OD1_HB_ ; }
  if( amino_and_atom == "ASN  ND2" ){ return  _ASN_ND2_    ; }  
  if( amino_and_atom == "ASN BND2" ){ return  _ASN_ND2_HB_ ; }
  if( amino_and_atom == "ASP  CB " ){ return  _ASP_CB_     ; }   
  if( amino_and_atom == "ASP  CG " ){ return  _ASP_CG_     ; }   
  if( amino_and_atom == "ASP  OD1" ){ return  _ASP_OD1_    ; }  
  if( amino_and_atom == "ASP BOD1" ){ return  _ASP_OD1_HB_ ; }
  if( amino_and_atom == "ASP  OD2" ){ return  _ASP_OD2_    ; }  
  if( amino_and_atom == "ASP BOD2" ){ return  _ASP_OD2_HB_ ; }
  if( amino_and_atom == "CYS  CB " ){ return  _CYS_CB_     ; }   
  if( amino_and_atom == "CYS  SG " ){ return  _CYS_SG_     ; }   
  if( amino_and_atom == "GLN  CB " ){ return  _GLN_CB_     ; }   
  if( amino_and_atom == "GLN  CG " ){ return  _GLN_CG_     ; }   
  if( amino_and_atom == "GLN  CD " ){ return  _GLN_CD_     ; }   
  if( amino_and_atom == "GLN  OE1" ){ return  _GLN_OE1_    ; }  
  if( amino_and_atom == "GLN BOE1" ){ return  _GLN_OE1_HB_ ; }
  if( amino_and_atom == "GLN  NE2" ){ return  _GLN_NE2_    ; }  
  if( amino_and_atom == "GLN BNE2" ){ return  _GLN_NE2_HB_ ; }
  if( amino_and_atom == "GLU  CB " ){ return  _GLU_CB_     ; }   
  if( amino_and_atom == "GLU  CG " ){ return  _GLU_CG_     ; }   
  if( amino_and_atom == "GLU  CD " ){ return  _GLU_CD_     ; }   
  if( amino_and_atom == "GLU  OE1" ){ return  _GLU_OE1_    ; }  
  if( amino_and_atom == "GLU BOE1" ){ return  _GLU_OE1_HB_ ; }
  if( amino_and_atom == "GLU  OE2" ){ return  _GLU_OE2_    ; }  
  if( amino_and_atom == "GLU BOE2" ){ return  _GLU_OE2_HB_ ; }
  if( amino_and_atom == "HIS  CB " ){ return  _HIS_CB_     ; }   
  if( amino_and_atom == "HIS  CG " ){ return  _HIS_CG_     ; }   
  if( amino_and_atom == "HIS  ND1" ){ return  _HIS_ND1_    ; }  
  if( amino_and_atom == "HIS BND1" ){ return  _HIS_ND1_HB_ ; }
  if( amino_and_atom == "HIS  CD2" ){ return  _HIS_CD2_    ; }  
  if( amino_and_atom == "HIS  CE1" ){ return  _HIS_CE1_    ; }  
  if( amino_and_atom == "HIS  NE2" ){ return  _HIS_NE2_    ; }  
  if( amino_and_atom == "HIS BNE2" ){ return  _HIS_NE2_HB_ ; }
  if( amino_and_atom == "ILE  CB " ){ return  _ILE_CB_     ; }   
  if( amino_and_atom == "ILE  CG1" ){ return  _ILE_CG1_    ; }  
  if( amino_and_atom == "ILE  CG2" ){ return  _ILE_CG2_    ; }  
  if( amino_and_atom == "ILE  CD1" ){ return  _ILE_CD1_    ; }  
  if( amino_and_atom == "LEU  CB " ){ return  _LEU_CB_     ; }   
  if( amino_and_atom == "LEU  CG " ){ return  _LEU_CG_     ; }   
  if( amino_and_atom == "LEU  CD1" ){ return  _LEU_CD1_    ; }  
  if( amino_and_atom == "LEU  CD2" ){ return  _LEU_CD2_    ; }  
  if( amino_and_atom == "LYS  CB " ){ return  _LYS_CB_     ; }   
  if( amino_and_atom == "LYS  CG " ){ return  _LYS_CG_     ; }   
  if( amino_and_atom == "LYS  CD " ){ return  _LYS_CD_     ; }   
  if( amino_and_atom == "LYS  CE " ){ return  _LYS_CE_     ; }   
  if( amino_and_atom == "LYS  NZ " ){ return  _LYS_NZ_     ; }   
  if( amino_and_atom == "LYS BNZ " ){ return  _LYS_NZ_HB_  ; }
  if( amino_and_atom == "MET  CB " ){ return  _MET_CB_     ; }   
  if( amino_and_atom == "MET  CG " ){ return  _MET_CG_     ; }   
  if( amino_and_atom == "MET  SD " ){ return  _MET_SD_     ; }   
  if( amino_and_atom == "MET  CE " ){ return  _MET_CE_     ; }   
  if( amino_and_atom == "PHE  CB " ){ return  _PHE_CB_     ; }   
  if( amino_and_atom == "PHE  CG " ){ return  _PHE_CG_     ; }   
  if( amino_and_atom == "PHE  CD1" ){ return  _PHE_CD1_    ; }  
  if( amino_and_atom == "PHE  CD2" ){ return  _PHE_CD2_    ; }  
  if( amino_and_atom == "PHE  CE1" ){ return  _PHE_CE1_    ; }  
  if( amino_and_atom == "PHE  CE2" ){ return  _PHE_CE2_    ; }  
  if( amino_and_atom == "PHE  CZ " ){ return  _PHE_CZ_     ; }   
  if( amino_and_atom == "PRO  CB " ){ return  _PRO_CB_     ; }   
  if( amino_and_atom == "PRO  CG " ){ return  _PRO_CG_     ; }   
  if( amino_and_atom == "PRO  CD " ){ return  _PRO_CD_     ; }   
  if( amino_and_atom == "SER  CB " ){ return  _SER_CB_     ; }   
  if( amino_and_atom == "SER  OG " ){ return  _SER_OG_     ; }   
  if( amino_and_atom == "SER BOG " ){ return  _SER_OG_HB_  ; }
  if( amino_and_atom == "THR  CB " ){ return  _TRH_CB_     ; }   
  if( amino_and_atom == "THR  OG1" ){ return  _TRH_OG1_    ; }  
  if( amino_and_atom == "THR BOG1" ){ return  _TRH_OG1_HB_ ; } 
  if( amino_and_atom == "THR  CG2" ){ return  _TRH_CG2_    ; }  
  if( amino_and_atom == "TRP  CB " ){ return  _TRP_CB_     ; }   
  if( amino_and_atom == "TRP  CG " ){ return  _TRP_CG_     ; }   
  if( amino_and_atom == "TRP  CD1" ){ return  _TRP_CD1_    ; }  
  if( amino_and_atom == "TRP  CD2" ){ return  _TRP_CD2_    ; }  
  if( amino_and_atom == "TRP  NE1" ){ return  _TRP_NE1_    ; }  
  if( amino_and_atom == "TRP  CE2" ){ return  _TRP_CE2_    ; }  
  if( amino_and_atom == "TRP  CE3" ){ return  _TRP_CE3_    ; }  
  if( amino_and_atom == "TRP  CZ2" ){ return  _TRP_CZ2_    ; }  
  if( amino_and_atom == "TRP  CZ3" ){ return  _TRP_CZ3_    ; }  
  if( amino_and_atom == "TRP  CH2" ){ return  _TRP_CH2_    ; }  
  if( amino_and_atom == "TYR  CB " ){ return  _TYR_CB_     ; }   
  if( amino_and_atom == "TYR  CG " ){ return  _TYR_CG_     ; }   
  if( amino_and_atom == "TYR  CD1" ){ return  _TYR_CD1_    ; }  
  if( amino_and_atom == "TYR  CD2" ){ return  _TYR_CD2_    ; }  
  if( amino_and_atom == "TYR  CE1" ){ return  _TYR_CE1_    ; }  
  if( amino_and_atom == "TYR  CE2" ){ return  _TYR_CE2_    ; }  
  if( amino_and_atom == "TYR  CZ " ){ return  _TYR_CZ_     ; }   
  if( amino_and_atom == "TYR  OH " ){ return  _TYR_OH_     ; }   
  if( amino_and_atom == "TYR BOH " ){ return  _TYR_OH_HB_  ; }
  if( amino_and_atom == "VAL  CB " ){ return  _VAL_CB_     ; }   
  if( amino_and_atom == "VAL  CG1" ){ return  _VAL_CG1_    ; }  
  if( amino_and_atom == "VAL  CG2" ){ return  _VAL_CG2_    ; }  
  return  _NON_ATOM_ ; 
}
/*============================================================*/
bool is_HB( dmd_atom_t &type ){
  switch( type ){
  case _BKB_N_HB_:    return true ;
  case _BKB_O_HB_:    return true ;
  case _ARG_NH1_HB_:  return true ;
  case _ARG_NH2_HB_:  return true ;
  case _ASN_OD1_HB_:  return true ;
  case _ASN_ND2_HB_:  return true ;
  case _ASP_OD1_HB_:  return true ;
  case _ASP_OD2_HB_:  return true ;
  case _GLN_OE1_HB_:  return true ;
  case _GLN_NE2_HB_:  return true ;
  case _GLU_OE1_HB_:  return true ;
  case _GLU_OE2_HB_:  return true ;
  case _HIS_ND1_HB_:  return true ;
  case _HIS_NE2_HB_:  return true ;
  case _LYS_NZ_HB_:   return true ;
  case _SER_OG_HB_:   return true ;
  case _TRH_OG1_HB_:  return true ;
  case _TYR_OH_HB_:   return true ;
  default : return false ;
  }
}
/*=========================================================*/
bool belong_to_same_XXX( dmd_atom_t &type1, dmd_atom_t &type2 ){
  string stype1 = dmd_atom_t2string( type1 );
  string stype2 = dmd_atom_t2string( type2 );
  if( stype1.substr( 0, 3 ) == stype2.substr( 0, 3 ) ){
    return true ;
  }
  return false ;
}
/*============================================================*/
string dmd_atom_t2amino( dmd_atom_t &type){
  string stype = dmd_atom_t2string( type );
  return stype.substr( 0, 3 ) ;
}
/*============================================================*/
string dmd_atom_t2at( dmd_atom_t &type ){
  string stype = dmd_atom_t2string( type );
  return stype.substr( 4, 4 ) ;
}
/*============================================================*/
string dmd_atom_t2stringII( dmd_atom_t &type1, dmd_atom_t &type2 ){
  string stype1 = dmd_atom_t2string( type1 );
  string stype2 = dmd_atom_t2string( type2 );
  string amino  = stype1.substr( 0, 3 ) ;
  string at1    = stype1.substr( 4, 4 ) ;
  string at2    = stype2.substr( 4, 4 ) ;
}
/*============================================================*/
void goto_txt_key( ifstream &DMD, const dmd_txt_key &key ){
  string signal = txtKeyWords[ key ], line ;
  DMD.seekg( 0 ) ;
  do{ 
    getline( DMD, line ) ;   /*cout<<"line="<<line<<endl;*/
  }while( line != signal ) ;
}
/*============================================================*/
