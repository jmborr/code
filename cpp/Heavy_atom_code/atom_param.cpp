using namespace std;
#include<iostream>
#include<string>
#include<fstream.h>
#include"atom_param.h"
#include"miscellanea.h"
#include"bibliography.h"

/*===================================================================*/
bool is_HB( string &amino_and_atom ){
    if( amino_and_atom == "BKB BN  " ){ return  true   ; } 
    if( amino_and_atom == "BKB BO  " ){ return  true   ; } 
    if( amino_and_atom == "ARG BNH1" ){ return  true ; }
    if( amino_and_atom == "ARG BNH2" ){ return  true ; }
    if( amino_and_atom == "ASN BOD1" ){ return  true ; }
    if( amino_and_atom == "ASN BND2" ){ return  true ; }
    if( amino_and_atom == "ASP BOD1" ){ return  true ; }
    if( amino_and_atom == "ASP BOD2" ){ return  true ; }
    if( amino_and_atom == "GLN BOE1" ){ return  true ; }
    if( amino_and_atom == "GLN BNE2" ){ return  true ; }
    if( amino_and_atom == "GLU BOE1" ){ return  true ; }
    if( amino_and_atom == "GLU BOE2" ){ return  true ; }
    if( amino_and_atom == "HIS BND1" ){ return  true ; }
    if( amino_and_atom == "HIS BNE2" ){ return  true ; }
    if( amino_and_atom == "LYS BNZ " ){ return  true  ; }
    if( amino_and_atom == "SER BOG " ){ return  true  ; }
    if( amino_and_atom == "THR BOG1" ){ return  true ; } 
    if( amino_and_atom == "TYR BOH " ){ return  true  ; }
    return false ;
}
/*===================================================================*/
bool is_comment( string &line ){
  if( line.substr( 0, 1 ) == "#" ){ return true ; }
  return false ;
}
/*===================================================================*/
bool has_comment( string &line ){
  if( line.find("#") != string::npos ){ return true ; }
  return false ;
}
/*===================================================================*/
bool goto_dat_key( ifstream &DAT, const dat_key &key ){
  int pos=DAT.tellg( ) ;
  string line, signal("#") ; 
  signal += datKeyWords[ key ] ; /*cout<<"signal=|"<<signal<<"|\n";*/
  DAT.seekg( 0 ) ;
  do{ 
    if( !getline( DAT, line ) ){
      cerr<<"\nERROR: no key \""<<datKeyWords[ key ]<<"\" found\n";
      DAT.clear( ios::goodbit ) ; DAT.seekg( pos ) ;
      return false ;
    }
  }while( (line != signal) ) ;
  return true ;
 }
/*===================================================================*/
bool remove_trailing_comment( string &line ){
  int x = static_cast<int>(line.find("#") );
  /*cout <<"x="<<x<<endl;*/
  if(x != -1){line.erase( x ) ;} /*will erase from "#" to the end of line*/
}
/*===================================================================*/
bool is_end( string &line ){
  if( line.substr( 0, 4 ) == "#END" ){ return true ; }
  return false ;
}
/*===================================================================*/
bool is_end( string &line, const dat_key &key  ){
  string end("#END ") ; end += datKeyWords[ key ] ;
  if( line == end ){ return true ; }
  return false ;
}
/*===================================================================*/
void print_atom_type_line( char *buf, ifstream &DAT, string &amino_and_atom ){
  string sbuf ;
  char cbuf[256] ;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  string line ;
  /*to_match =  dmd_atom_t2string( type ) ;*/
  do{
    if( !getline(DAT, line) ){
      cout << "ERROR: no \""<<amino_and_atom<<"\" found in the force field!\n";
      exit(1);
    }
  }while( line.find( amino_and_atom ) == string::npos ) ;
  double mass, hcr ;
  sbuf = line.substr(  9, 8 ) ;
  string2char( sbuf, cbuf ) ; mass = atof( cbuf ) ;
  sbuf = line.substr( 17, 6 ) ;
  string2char( sbuf, cbuf ) ; hcr  = atof( cbuf ) ;
  sprintf( buf, "%9.6f %9.6f %9.6f\n", mass, hcr, hcr*1.0001 ) ;
  DAT.seekg( pos ) ;
}
/*===================================================================*/
double out_hard_core_radius( ifstream &DAT, string &amino_and_atom ){
  string sbuf, line ;
  char cbuf[256] ;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  do{
    if( !getline(DAT, line) ){
      cout<<"error line="<<line<<endl;
      cout << "ERROR: no \""<<amino_and_atom<<"\" found in force field!\n" ;
      exit(1);
    }
  }while( line.find( amino_and_atom ) == string::npos ) ;
  double hcr ;
  sbuf = line.substr( 17, 6 ) ;
  string2char( sbuf, cbuf ) ; hcr  = atof( cbuf ) ;
  DAT.seekg( pos ) ;
  return hcr ;
}
/*===================================================================*/
double out_charge( ifstream &DAT, string &amino_and_atom ){
  string sbuf, line ;
  char cbuf[256] ;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  do{
    if( !getline(DAT, line) ){
      cout<<"error line="<<line<<endl;
      cout << "ERROR: no \""<<amino_and_atom<<"\" found in force field!\n" ;
      exit(1);
    }
  }while( line.find( amino_and_atom ) == string::npos ) ;
  double q ;
  sbuf = line.substr( 24, 7 ) ;
  string2char( sbuf, cbuf ) ; q  = atof( cbuf ) ;
  DAT.seekg( pos ) ;
  return q ;
}
/*===================================================================*/
double out_hydrophobicity( ifstream &DAT, string &amino_and_atom ){
  string sbuf, line ;
  char cbuf[256] ;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  do{
    if( !getline(DAT, line) ){
      cout<<"error line="<<line<<endl;
      cout << "ERROR: no \""<<amino_and_atom<<"\" found in force field!\n" ;
      exit(1);
    }
  }while( line.find( amino_and_atom ) == string::npos ) ;
  double hb ;
  sbuf = line.substr( 32, 7 ) ;
  string2char( sbuf, cbuf ) ; hb  = atof( cbuf ) ;
  DAT.seekg( pos ) ;
  return hb ;
}
/*===================================================================*/
bool is_HB_link( string &link, ifstream &DAT ){
  string line, link_clone = link ;
  remove_trailing_comment( link_clone ) ;
  remove_firsts( " ", link_clone ) ;
  remove_lasts( " ", link_clone ) ; 
  /*cout<<"clone=|"<<link<<"|\n";*/
  int pos = DAT.tellg( ) ; DAT.seekg( 0 ) ;
  if( goto_dat_key( DAT, LINK_HB_DAT) ){
    getline( DAT, line ) ;
    while( !is_end( line, LINK_HB_DAT ) ){
      if( !is_comment(line) ){
        remove_trailing_comment( line ) ;
	remove_firsts( " ", line ) ;
	remove_lasts( " ", line ) ; 
	if( line == link_clone ){	  
	  DAT.seekg( pos ) ;
	  return true ;
	}
      }
      getline( DAT, line ) ;
    }  
  }
  DAT.seekg( pos ) ;
  return false ;
}
/*===================================================================*/
bool is_there_linked_pair( ifstream &DAT, string &aa_at1, string &aa_at2,
			   string &link ){
  int pos = DAT.tellg( ) ; DAT.seekg( 0 ) ;
  string line, aa_atB1, aa_atB2 ;
  dat_key types_of_links[ ]={ LINK_PAIRS_DAT, LINK_HB_DAT } ;
  int num_types_of_links = 2 ;
  
  /*getline( DAT, line ) ; cout<<"line=|"<<line<<"|\n";exit(1);*/
  for( int z=0; z<num_types_of_links; z++ ){
    if( goto_dat_key( DAT, types_of_links[ z ] ) ){
      getline( DAT, line ) ; /*cout<<"line=|"<<line<<"|\n";exit(1);*/
      while( !is_end( line, types_of_links[ z ] ) ){
	if( !is_comment( line ) ){
	  /*cout<<"line=|"<<line<<"|\n";exit(1);*/
	  aa_atB1 = line.substr( 0, 8 ) ;
	  aa_atB2 = line.substr( 9, 8 ) ;
	  if( (aa_atB1==aa_at1 && aa_atB2==aa_at2) ||
	      (aa_atB1==aa_at2 && aa_atB2==aa_at1)   ){
	    remove_trailing_comment( line ) ;
	    link = line ;
	    return true ;
	  }
	}
	getline( DAT, line ) ;      
      }
    }
    DAT.seekg( 0 ) ;
  }

  DAT.seekg( pos ) ;
  return false ;
}
/*===================================================================*/
bool out_linked_pair( string &link, string &aa_at1, string &aa_at2, 
		       char *cbuf ){
  if( !is_comment(link) ){
    remove_trailing_comment(link) ;
    string sbuf ;
    aa_at1 = link.substr( 0, 8 ) ;
    aa_at2= link.substr( 9, 8 ) ;
    sbuf = link.substr( 17, link.length( )-17 ) ;
    string2char( sbuf, cbuf ) ;
    return true ;
  }
  return false ;
}
/*===================================================================*/
bool  out_linked_pairII( ifstream &DAT, string &aa_at1,
		       string &aa_at2, char *cbuf ){
  string link ;
  getline(DAT, link) ;
  return out_linked_pair( link, aa_at1, aa_at2, cbuf ) ;
}
/*===================================================================*/
bool out_prerelaxed_linked_pair( string &link, string &aa_at1, string &aa_at2, 
		       char *cbuf ){
  if( !is_comment(link) ){
    /*cout<<"link=|"<<link<<"|\n" ;*/
    remove_trailing_comment(link) ;  /*cout<<"link=|"<<link<<"|\n" ;*/
    char cbuf2[256] ;
    string sbuf, *snumbers ;
    char *cnumbers ;
    int n_numbers ;
    double low_hcd, high_hcd ; 
    aa_at1 = link.substr( 0, 8 ) ;
    aa_at2= link.substr( 9, 8 ) ;
    if( !is_HB( aa_at1 ) && !is_HB( aa_at1 ) ){
      sbuf=link.substr(17,link.length( )-17);/*cout<<"sbuf=|"<<sbuf<<"|\n" ;*/
      remove_firsts( " ", sbuf ) ;
      remove_lasts( " ", sbuf ) ; /*cout<<"sbuf=|"<<sbuf<<"|\n" ;*/
      snumbers = split_nonDegenerated( " ", sbuf, n_numbers ) ;
      sbuf.assign(" ") ;
      if( n_numbers > 2 ){/*There are internal barriers*/
	for( int i=2; i<n_numbers-1 ; i+=2 ){
	  if( snumbers[i].find( s_inf ) == string::npos ){
	    snumbers[i].assign( "0.000" ) ;/*remove finite energy steps*/
	  }
	  sbuf += snumbers[i-1] + " " + snumbers[i] + " " ;
	}
      }
      /*cout<<"sbuf=|"<<sbuf<<"|\n" ;*/
      string2char( snumbers[0], cbuf2 ) ; low_hcd = atof( cbuf2 ) ;
      string2char( snumbers[n_numbers-1], cbuf2 ) ; high_hcd = atof( cbuf2 ) ;
      /*cout<<"low_hdc="<<low_hcd<<" high_hcd="<<high_hcd<<endl;*/
      string2char( sbuf, cbuf2 ) ;

      sprintf( cbuf, "%7.5f %7.5f %7.5f %7.5f %7.5f %s %7.5f %7.5f %7.5f %7.5f %7.5f",
	       low_hcd*0.5, low_hcd*0.85, static_cast<double>(_AP_INF_),
	       low_hcd, static_cast<double>(_AP_INF_),
	       cbuf2,
	       high_hcd, -1*static_cast<double>(_AP_INF_),
	       high_hcd*1.15, -1*static_cast<double>(_AP_INF_),high_hcd*1.25 );
      return true ;
    }
  }
  return false ;
}
/*===================================================================
  Bonds between the atom and atoms belonging to same amino acid*/
string *out_all_same_permanent_bonds( string &aa_at, string &resName, 
					       ifstream &DAT, int &nBonds ){
  /*cout<<"aa_at="<<aa_at<<" resName="<<resName<<endl;*/
  string line, *friends = NULL, begin("#"), end("#END ") ;
  string *sbuf = new string[100] ;
  begin.append( resName ) ;
  end.append( resName ) ;  /*cout<<"begin="<<begin<<" end="<<end<<endl;*/
  nBonds = 0 ;
  int pos = DAT.tellg( ) ;
  DAT.seekg( 0 ) ;
  do{ 
    getline( DAT, line ) ;
  } while( line.find( datKeyWords[LIST_PERM_BONDS_DAT] ) == string::npos ) ;
  if( aa_at.substr( 0, 3 ) == "BKB" ){
    do{ 
      getline( DAT, line ) ;
    } while( line.find( "#BKB") == string::npos ) ;
    while(line.find( "#END BKB") == string::npos ){
      if( !is_comment(line) && line.substr( 0, 8 )==aa_at &&
	  line.substr( 18, 1 )=="S" ){
	sbuf[ nBonds ] = line.substr( 9, 8 ) ;
	nBonds++ ;
      }
      getline( DAT, line ) ;
    }
  }
  do{ 
    getline( DAT, line ) ;
  } while( line.find( begin ) == string::npos ) ;
  while(line.find( end ) == string::npos ){
    if( !is_comment(line) && line.substr( 0, 8 )==aa_at &&
	line.substr( 18, 1 )=="S" ){
      sbuf[ nBonds ] = line.substr( 9, 8 ) ;
      nBonds++ ;
    }
    getline( DAT, line ) ;
  }
  if( nBonds ){
    friends = new string[ nBonds ] ;
    for( int i=0; i<nBonds; i++ ){ friends[ i ] = sbuf[ i ] ; }
  }
  if( DAT.eof() ){ DAT.clear( ios::goodbit ) ; }
  DAT.seekg( pos ) ;  
  return friends ;
}
/*===================================================================
  Bonds between the atom and atoms belonging to next amino acid*/
string *out_all_next_permanent_bonds( string &aa_at, string &resName,
			 string &resName2, ifstream &DAT, int &nBonds ){
  string line, *friends = NULL, begin("#"), end("#END ") ;
  string *sbuf = new string[100] ;
  nBonds = 0 ;
  int pos = DAT.tellg( ) ;
  DAT.seekg( 0 ) ;
  do{ 
    getline( DAT, line ) ;
  } while( line.find( datKeyWords[LIST_PERM_BONDS_DAT] ) == string::npos ) ;
  if( aa_at.substr( 0, 3 ) == "BKB" ){
    /*cout<<"aa_at=|"<<aa_at<<"|\n";*/
    do{ 
      getline( DAT, line ) ;
    } while( line.find( "#BKB") == string::npos ) ;
    while(line.find( "#END BKB") == string::npos ){/*BKB_i--BKB_i+1*/
      if( !is_comment(line) && line.substr( 0, 8 )==aa_at &&
	  line.substr( 18, 1 )=="N" ){
	/*cout<<"line ="<<line<<endl;*/
	sbuf[ nBonds ] = line.substr( 9, 8 ) ;
	nBonds++ ;
      }
      getline( DAT, line ) ;
    }
    begin.append( resName2 ) ;
    end.append( resName2 ) ;/*cout<<"begin="<<begin<<" end="<<end<<endl;*/
    do{ 
      getline( DAT, line ) ;
    } while( line.find( begin ) == string::npos ) ;
    while(line.find( end ) == string::npos ){/*BKB_i---SIDECHAIN_i+1*/
      if( !is_comment(line) && line.substr( 0, 8 )==aa_at &&
	  line.substr( 18, 1 )=="N" ){
	/*cout<<"line2="<<line<<endl;*/
	sbuf[ nBonds ] = line.substr( 9, 8 ) ;
	nBonds++ ;
      }
      getline( DAT, line ) ;
    }
  }/*matches if( aa_at.substr( 0, 3 ) == "BKB" ) */
  else{/*atom belongs to the sidechain, not the backbone*/
    do{ 
      getline( DAT, line ) ;
    } while( line.find( begin ) == string::npos ) ;
    begin.append( resName ) ;
    end.append( resName ) ;
    while(line.find( end )==string::npos ){/*SIDECHAIN_i---BKB_i+1 (note: there are no bonds between one atom of sidechain i and one atom of sidechain i+1)*/
      if( !is_comment(line) && line.substr( 0, 8 )==aa_at &&
	  line.substr( 18, 2 )=="N" ){
	sbuf[ nBonds ] = line.substr( 9, 8 ) ;
	nBonds++ ;
      }
      getline( DAT, line ) ;
    }
  }  
  if( nBonds ){
    friends = new string[ nBonds ] ;
    for( int i=0; i<nBonds; i++ ){ 
      friends[ i ] = sbuf[ i ] ; 
      /*cout<<"friends["<<i<<"]="<<friends[ i ]<<endl;*/
    }
  }
  if( DAT.eof() ){ DAT.clear( ios::goodbit ) ; }
  DAT.seekg( pos ) ;  
  return friends ;
}
/*===================================================================*/
double get_inner_VW_factor( ifstream &DAT ){
  int pos = DAT.tellg( ), n=0 ; 
  DAT.seekg( 0 ) ;
  double f = 0.0, *dtmp=NULL ;
  string line ;

  goto_dat_key( DAT, VAN_DER_WAALS_DAT ) ;  getline( DAT, line ) ;
  while( !is_end( line, VAN_DER_WAALS_DAT ) ){
    if( !is_comment( line ) ){
      remove_trailing_comment( line ) ; 
      remove_firsts( " ", line ) ; remove_lasts( " ", line ) ;
      dtmp = split_line2doubles( line, n ) ;
      if( dtmp ){ f = dtmp[ 0 ] ;  break ; }
      else{ 
	cerr<<"ERROR with the default Van der Waals line\n" ; 
	exit( 1 ) ;
      }
    }
    getline( DAT, line ) ;
  }
  DAT.seekg( pos ) ;
  if( f <= 0.0 ){
    cout <<"ERROR: FORCE FIELD DOES NOT CONTAIN VAN DER WAALS INTERACTION\n" ;
    exit( 1 ) ;
  }
  DAT.seekg( pos ) ;
  return f ;
}
/*===================================================================*/
void get_default_Van_der_Waals( string &type1, string &type2,
				string &line, ifstream &DAT ){
  /* cout<<"void get_default_Van_der_Waals(string&,string&,string&,ifstream&)\n";*/
  double r1 = out_hard_core_radius( DAT, type1 ) ;
  double r2 = out_hard_core_radius( DAT, type2 ) ;
  double hcd = r1 + r2 ;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  string sbuf ;
  char cbuf[64] ;
  double *numbers = NULL ;
  int n_numbers = 0 ;
  goto_dat_key( DAT, VAN_DER_WAALS_DAT ) ;
  /*read the default Van der Waals line*/
  do{ getline( DAT, sbuf ) ; }while( is_comment( sbuf ) ) ;
  remove_trailing_comment( sbuf ) ;
  numbers = split_line2doubles( sbuf, n_numbers ) ;
  if( numbers ){
    sbuf.erase( 0 ) ;  line.erase( 0 ) ;  numbers[ 0 ] *= hcd ;
    for( int i=1; i<n_numbers; i+=2 ){ numbers[ i ] *= hcd ; }
    for( int i=0; i<n_numbers-1; i++ ){
      sprintf( cbuf, "%10.6f ", numbers[ i ] ) ;
      sbuf.assign( cbuf ) ; line += sbuf ;
    }
    sprintf( cbuf, "%10.6f", numbers[ n_numbers-1 ] ) ;
    sbuf.assign( cbuf ) ; line += sbuf ;
  }
  else{ cout<<"ERROR: no default Van der Waals line detected!\n"; exit(1) ; }
  DAT.seekg( pos ) ;
}
/*===================================================================*/
bool get_default_hydrophobicity( string &type1, string &type2,
				string &line, ifstream &DAT ){
  /*cout<<"void get_default_hydrophobicity(string&,string&,string&,ifstream&)\n";*/
  double q1 = out_hydrophobicity( DAT, type1 ) ;
  double q2 = out_hydrophobicity( DAT, type2 ) ;
  
  double r1 = out_hard_core_radius( DAT, type1 ) ;
  double r2 = out_hard_core_radius( DAT, type2 ) ;
  double q1q2 = q1+q2 ;
  if(fabs(q1)<0.01 || fabs(q2)<0.01){ return false ; }
  double hcd = r1 + r2 ;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  string sbuf ;
  char cbuf[64] ;
  double *numbers = NULL ;
  int n_numbers = 0 ;
  if( goto_dat_key( DAT, HYDROPHOBICITY_DAT ) ){
    do{ getline( DAT, sbuf ) ; }while( is_comment( sbuf ) ) ;
    remove_trailing_comment( sbuf ) ;
    numbers = split_line2doubles( sbuf, n_numbers ) ;
    if( numbers ){
      sbuf.erase( 0 ) ;  line.erase( 0 ) ;  numbers[ 0 ] *= hcd ;
      for( int i=1; i<n_numbers; i+=2 ){ numbers[ i ] *= hcd ; }
      for( int i=4; i<n_numbers; i+=2 ){ numbers[ i ] *= q1q2 ; }
      /*numbers[2] is potential barrier due to Van deer Walls repulsion, thus
	is independent of q1q1*/
      for( int i=0; i<n_numbers-1; i++ ){
	sprintf( cbuf, "%10.6f ", numbers[ i ] ) ;
	sbuf.assign( cbuf ) ; line += sbuf ;
      }
      sprintf( cbuf, "%10.6f", numbers[ n_numbers-1 ] ) ;
      sbuf.assign( cbuf ) ; line += sbuf ;
    }
    else{ cout<<"ERROR: no default hydrophobicity line detected!\n"; exit(1); }
    DAT.seekg( pos ) ;  return true ;
  }
  return false ;
}
/*===================================================================*/
bool get_default_hydrophobicity( const double &q1, const double &q2,
				 const double &r1, const double &r2,
				 string &line, ifstream &DAT        ){
  double q1q2 = q1+q2 ;
  if(fabs(q1)<0.01 || fabs(q2)<0.01){ return false ; }
  double hcd = r1 + r2 ;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  string sbuf ;
  char cbuf[64] ;
  double *numbers = NULL ;
  int n_numbers = 0 ;
  if( goto_dat_key( DAT, HYDROPHOBICITY_DAT ) ){
    do{ getline( DAT, sbuf ) ; }while( is_comment( sbuf ) ) ;
    remove_trailing_comment( sbuf ) ;
    numbers = split_line2doubles( sbuf, n_numbers ) ;
    if( numbers ){
      sbuf.erase( 0 ) ;  line.erase( 0 ) ;  numbers[ 0 ] *= hcd ;
      for( int i=1; i<n_numbers; i+=2 ){ numbers[ i ] *= hcd ; }
      for( int i=4; i<n_numbers; i+=2 ){ numbers[ i ] *= q1q2 ; }
      /*numbers[2] is potential barrier due to Van deer Walls repulsion, thus
	is independent of q1q1*/
      for( int i=0; i<n_numbers-1; i++ ){
	sprintf( cbuf, "%10.6f ", numbers[ i ] ) ;
	sbuf.assign( cbuf ) ; line += sbuf ;
      }
      sprintf( cbuf, "%10.6f", numbers[ n_numbers-1 ] ) ;
      sbuf.assign( cbuf ) ; line += sbuf ;
    }
    else{ cout<<"ERROR: no default hydrophobicity line detected!\n"; exit(1); }
    DAT.seekg( pos ) ;  return true ;
  }
  return false ;
}
/*===================================================================*/
void rescale_default_hydrophobicity_by_factor( string &line, const double &f ){
  double *numbers = NULL ;
  int n_numbers = 0 ;
  string sbuf ;
  char cbuf[64] ;
  remove_trailing_comment( line ) ;
  numbers = split_line2doubles( line, n_numbers ) ;
  if( numbers ){
    line.erase( 0 ) ;
      for( int i=4; i<n_numbers; i+=2 ){ numbers[ i ] *= f ; }
      /*numbers[2] is potential barrier due to Van deer Walls repulsion, thus
	is independent of q1q1*/
      for( int i=0; i<n_numbers-1; i++ ){
	sprintf( cbuf, "%10.6f ", numbers[ i ] ) ;
	sbuf.assign( cbuf ) ; line += sbuf ;
      }
      sprintf( cbuf, "%10.6f", numbers[ n_numbers-1 ] ) ;
      sbuf.assign( cbuf ) ; line += sbuf ;
  }
}
/*===================================================================*/
double get_HB_strength( ifstream &DAT ){
  double *numbers=NULL ;
  int n_numbers, pos = DAT.tellg( ) ; DAT.seekg( 0 ) ;
  string line ;
  if( ! goto_dat_key( DAT, HB_STRENGTH_DAT ) ){
    cerr<<"ERROR: no "<<datKeyWords[HB_STRENGTH_DAT]<<" section found!\n";
    exit(1);
  }
  do{ getline( DAT, line ) ; }while( is_comment( line ) ) ;
  remove_trailing_comment( line ) ;
  numbers = split_line2doubles( line, n_numbers ) ;
  if( numbers ){ DAT.seekg( pos ) ; return numbers[0] ; }
  else{ cerr<<"hydrogen bond strength entry is incorrect\n!" ; exit(1); }
}
/*===================================================================*/
void get_signal_Hydrogen_Bond( string &type1, string &type2,
			       string &line, ifstream &DAT ){
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  goto_dat_key( DAT, SIGNAL_HB_DAT ) ;
  string t1, t2, sbuf ;  getline( DAT, sbuf ) ;
  while( !is_end( sbuf, SIGNAL_HB_DAT ) ){
    if( !is_comment( sbuf ) ){
      t1 = sbuf.substr( 0, 8 ) ; t2 = sbuf.substr( 9, 8 ) ;
      if( t1==type1 && t2==type2 ){
	/*cout<<"sbuf="<<sbuf<<endl;*/
	line = sbuf.substr( 17, sbuf.size( )-17 ) ;
	return ;
      }
    }
    getline( DAT, sbuf ) ;
  }
  DAT.seekg( pos ) ;
}
/*===================================================================*/
bool do_they_react( ifstream &DAT, string &aa_at1, string &aa_at2,
		    string &reaction ){
  string line, aa_atB1, aa_atB2 ;

  if( goto_dat_key( DAT, REACT_DAT ) ){
    getline( DAT, line ) ;
    while( !is_end( line, REACT_DAT ) ){
      if( !is_comment( line ) ){
	aa_atB1 = line.substr( 0, 8 ) ;
	aa_atB2 = line.substr( 9, 8 ) ;
	if( (aa_atB1==aa_at1 && aa_atB2==aa_at2) ||
	    (aa_atB1==aa_at2 && aa_atB2==aa_at1)   ){
	    remove_trailing_comment( line ) ;
	    reaction = line ;
	    return true ;
	}
      }
      getline( DAT, line ) ;      
    }
  }
  return false ;
}
/*===================================================================*/
bool get_reactants_products(  string &line, string &sold1, string &sold2,
			      string &snew1, string &snew2 ){
  if( !is_comment( line ) ){
    sold1 = line.substr( 0,  8 ) ;
    sold2 = line.substr( 9,  8 ) ;
    snew1 = line.substr( 18, 8 ) ;
    snew2 = line.substr( 27, 8 ) ;
    return true ;
  }
  return false ;
}
/*===================================================================*/
bool get_associates( string &line, string *&associates, int *&rel_index,
		     string *&atoms, int &n ){
  int l_assoc = 12 ; /*12 is length of each associate plus indexes plus spaces*/
  char cbuf[4] ;
  string sbuf ;
  if( is_comment( line ) ){ return false ; }  /*cout<<"line="<<line<<endl;*/
  if( associates ){ delete [] associates; associates=NULL;}/*some clean up first*/
  if( rel_index ){ delete [] rel_index ; rel_index = NULL ; }
  if( atoms ){ delete [] atoms ; atoms = NULL ; }
  n = 0 ;
  remove_trailing_comment(line);  remove_firsts(" ",line);  remove_lasts(" ",line);
  /*cout<<"line=|"<<line<<"|\n";*/
  n = static_cast<int>((line.size( )+1.1)/l_assoc) ; /*cout<<"n="<<n<<endl;*/
  if( n ){
    associates = new string[n];   rel_index = new int[n];   atoms = new string[n];
    for( int i=0; i<n; i++ ){
      associates[ i ] = line.substr( i*l_assoc, 8 ) ;
      atoms[ i ] = associates[ i ].substr( 4, 4 ) ;
      sbuf = line.substr( i*l_assoc+9, 2 );   string2char( sbuf, cbuf ) ;
      rel_index[ i ] = atoi( cbuf ) ;
    }
    /*for(int i=0;i<n;i++){cout<<"associates["<<i<<"]=|"<<associates[ i ]<<"|, atoms["<<i<<"]=|"<<atoms[i]<<"|, rel_index["<<i<<"]="<<rel_index[ i ]<<endl;}*/
    return true ;
  }
  return false ;
}
/*===================================================================*/
bool have_rotamers( ifstream &DAT, string &type1, string &type2 ){
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  string sbuf, aa_at1, aa_at2 ;
  /*jump over the default rotamer barrier line*/
  if( !goto_dat_key( DAT, ROTAMER_BARRIERS_DAT ) ){ 
    DAT.seekg( pos ) ;return false ; 
  }
  do{ getline( DAT, sbuf ) ; }while( !is_comment( sbuf ) ) ;
  do{
    getline( DAT, sbuf ) ;
    if( !is_comment( sbuf ) ){
      aa_at1 = sbuf.substr(0, 8 ) ;
      aa_at2 = sbuf.substr(9, 8 ) ;
      if( (aa_at1==type1 && aa_at2==type2) ||
	  (aa_at2==type1 && aa_at1==type2)    ){ 
	DAT.seekg( pos ) ;  return true ; 
      }
    }
  }while( !is_end( sbuf, ROTAMER_BARRIERS_DAT ) ) ;
  DAT.seekg( pos ) ;  return false ;
}
/*===================================================================*/
double get_default_rotamer_barrier( ifstream &DAT ){
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  string sbuf ;
  char cbuf[64] ;
  goto_dat_key( DAT, ROTAMER_BARRIERS_DAT ) ;
  do{ getline( DAT, sbuf ) ; }while( is_comment( sbuf ) ) ;
  remove_trailing_comment(sbuf); remove_firsts(" ",sbuf); remove_lasts(" ",sbuf);
  string2char( sbuf, cbuf ) ;
  DAT.seekg( pos ) ;
  return atof( cbuf ) ;
}
/*===================================================================*/
void default_rotamer_barrier( ifstream &DAT, string &link, char * cbuf ){
  double bd = get_default_rotamer_barrier( DAT ) ; /*cout<<"bd="<<bd<<endl;*/
  double b ;
  string sbuf, *snumbers ;
  int n_numbers ;
  char cbuf2[256] ;
  if( !is_comment(link) ){
    remove_trailing_comment(link) ;
    sbuf=link.substr(17,link.length( )-17);
    remove_firsts( " ", sbuf ) ;
    remove_lasts( " ", sbuf ) ; /*cout<<"sbuf=|"<<sbuf<<"|\n" ;*/
    snumbers = split_nonDegenerated( " ", sbuf, n_numbers ) ;
    sbuf = "  " + snumbers[0] + "  " ;
    for( int i=2; i<n_numbers-1 ; i+=2 ){
      string2char( snumbers[i], cbuf2 ) ; b = atof( cbuf2 ) ; b *= bd ;
      sprintf( cbuf2, "%6.3f", b ) ; snumbers[i].assign( cbuf2 ) ;
      sbuf += snumbers[i-1] + "  " + snumbers[i] + "  " ;      
    }
    sbuf += snumbers[n_numbers-1] ; string2char( sbuf, cbuf ) ;
  }
}
/*===================================================================*/
double charge_of( string &stype, ifstream &DAT ){
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  char cbuf[64] ;
  string line, aa_at, q ;
  goto_dat_key( DAT, TYPE_ATOMS_DAT ) ;
  do{
    getline( DAT, line ) ;
    if( !is_comment( line ) ){
      aa_at = line.substr( 0, 8 ) ;
      if( aa_at == stype ) {
	q = line.substr( 24, 7 ) ; string2char( q, cbuf ) ;
	DAT.seekg( pos ) ;  return atof( cbuf ) ;	
      }
    }
  }while( !is_end( line, TYPE_ATOMS_DAT ) ) ;
  DAT.seekg( pos ) ;  return 0.0 ;
}
/*===================================================================*/
bool is_charged( string &stype, ifstream &DAT ){
  if( fabs( charge_of( stype, DAT ) )<0.001 ){ return false ; }
  return true ;
}
/*===================================================================*/
bool get_electrostatics( string &stype1, string &stype2, string &line, ifstream &DAT ){
  /*cout<<"void get_electrostatics(string&,string&,string&,ifstream&)\n";*/
  double q1 = charge_of( stype1, DAT ); if( fabs(q1)<0.001 ){return false;}
  double q2 = charge_of( stype2, DAT ); if( fabs(q2)<0.001 ){return false;}
  /*cout<<stype1<<"  "<<q1<<"  "<<stype2<<"  "<<q2<<endl;*/
  double q1q2 = q1 * q2 ;
  string sbuf, sbuf2, aa_at1, aa_at2 ;
  double *numbers = NULL ;
  int n_numbers = 0 ;
  char cbuf[ 64 ] ;
  int pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  getline(DAT,sbuf);  /*cout<<"sbuf="<<sbuf<<endl;*/
  if( !goto_dat_key( DAT, ELECTROSTATICS_DAT ) ){ 
    cerr<<"\nERROR: no key \""<<datKeyWords[ELECTROSTATICS_DAT]<<"\" found\n";
    DAT.seekg( pos ) ; return false ;
  } 
  /*jump over default electrostatic line*/
  do{getline(DAT,sbuf);/*cout<<"sbuf="<<sbuf<<endl;*/}while(is_comment(sbuf));
  
  do{
    getline( DAT, sbuf ) ;  
    if( !is_comment( sbuf ) ){
      aa_at1 = sbuf.substr( 0, 8 ) ; aa_at2 = sbuf.substr( 9, 8 ) ;
      if( (aa_at1==stype1 && aa_at2==stype2) ||
	  (aa_at1==stype2 && aa_at2==stype1)   ){
	remove_trailing_comment( sbuf ) ;
	remove_firsts( " ", sbuf ) ;  remove_lasts( " ", sbuf ) ;
	sbuf2 = sbuf.substr( 17, sbuf.length( )-17 ) ;
	numbers = split_line2doubles( sbuf2, n_numbers ) ;
	if( numbers ){
	  sbuf.erase( 0 ) ;  line.erase( 0 ) ;
	  for( int i=2; i<n_numbers; i+=2 ){
	    numbers[ i ] *= q1q2 ;
	    if(i==2){
	      numbers[2]=fabs(numbers[2]);
	    }/*make sure the barrier is positive*/
	  }
	  for( int i=0; i<n_numbers; i++ ){
	    sprintf( cbuf, "%10.6f  ", numbers[ i ] ) ;
	    sbuf.assign( cbuf ) ; line += sbuf ;	
	  }
	}
	DAT.seekg( pos ) ;  return true ;
      }
    }
  }while( !is_end( line, ELECTROSTATICS_DAT ) ) ;
  DAT.seekg( pos ) ;  return false ;
}
/*===================================================================*/
bool get_default_electrostatics( string &type1, string &type2,
				  string &line, ifstream &DAT ){
  double q1 = charge_of( type1, DAT ); if( fabs(q1)<0.001 ){return false;}
  double q2 = charge_of( type2, DAT ); if( fabs(q2)<0.001 ){return false;}
  double q1q2 = q1*q2,  *numbers = NULL ;
  double r1 = out_hard_core_radius( DAT, type1 ) ;
  double r2 = out_hard_core_radius( DAT, type2 ) ;
  double hcd = r1 + r2 ;
  int n_numbers, pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  string sbuf ;
  char cbuf[64] ;
  if( !goto_dat_key( DAT, ELECTROSTATICS_DAT ) ){
    cout<<"void get_default_electrostatics(string&,string&,string&,ifstream&)\n";
    cerr<<"\nERROR: no key \""<<datKeyWords[ELECTROSTATICS_DAT]<<"\" found\n";
    DAT.seekg( pos ) ; return false ;
  }
  do{getline(DAT,sbuf);/*cout<<"sbuf="<<sbuf<<endl;*/}while(is_comment(sbuf));
  remove_trailing_comment( sbuf ) ;
  numbers = split_line2doubles( sbuf, n_numbers ) ;
  if( numbers ){
    sbuf.erase( 0 ) ;  line.erase( 0 ) ;  numbers[ 0 ] *= hcd ;
    for( int i=1; i<n_numbers; i+=2 ){ numbers[ i ] *= hcd ; }
    for( int i=4; i<n_numbers; i+=2 ){ numbers[ i ] *= q1q2 ; }
    /*numbers[2] is potential barrier due to Van deer Walls repulsion, thus
      is independent of q1q1*/
    for( int i=0; i<n_numbers-1; i++ ){
      sprintf( cbuf, "%10.6f ", numbers[ i ] ) ;
      sbuf.assign( cbuf ) ; line += sbuf ;
    }
    sprintf( cbuf, "%10.6f", numbers[ n_numbers-1 ] ) ;
    sbuf.assign( cbuf ) ; line += sbuf ;
  }
  else{ cout<<"ERROR: no default electrostatics line detected!\n"; exit(1); }
  DAT.seekg( pos ) ;
  return true ;
}
/*===================================================================*/
bool get_default_electrostatics( const double &q1, const double &q2,
				 const double &r1, const double &r2,
				 string &line, ifstream &DAT        ){
  double q1q2 = q1*q2,  *numbers = NULL ;
  if( fabs(q1q2) < 0.001 ){ return false ; }
  double hcd = r1 + r2 ;
  int n_numbers, pos = DAT.tellg( ) ;  DAT.seekg( 0 ) ;
  string sbuf ;
  char cbuf[64] ;
  if( !goto_dat_key( DAT, ELECTROSTATICS_DAT ) ){
    cout<<"void get_default_electrostatics(string&,string&,string&,ifstream&)\n";
    cerr<<"\nERROR: no key \""<<datKeyWords[ELECTROSTATICS_DAT]<<"\" found\n";
    DAT.seekg( pos ) ; return false ;
  }
  do{getline(DAT,sbuf);/*cout<<"sbuf="<<sbuf<<endl;*/}while(is_comment(sbuf));
  remove_trailing_comment( sbuf ) ;
  numbers = split_line2doubles( sbuf, n_numbers ) ;
  if( numbers ){
    sbuf.erase( 0 ) ;  line.erase( 0 ) ;  numbers[ 0 ] *= hcd ;
    for( int i=1; i<n_numbers; i+=2 ){ numbers[ i ] *= hcd ; }
    for( int i=2; i<n_numbers; i+=2 ){ 
      numbers[ i ] *= q1q2 ;
      if(i==2){numbers[2]=fabs(numbers[2]);}/*assure interior barrier is positive*/
    }
    for( int i=0; i<n_numbers-1; i++ ){
      sprintf( cbuf, "%10.6f ", numbers[ i ] ) ;
      sbuf.assign( cbuf ) ; line += sbuf ;
    }
    sprintf( cbuf, "%10.6f", numbers[ n_numbers-1 ] ) ;
    sbuf.assign( cbuf ) ; line += sbuf ;
  }
  else{ cout<<"ERROR: no default electrostatics line detected!\n"; exit(1); }
  DAT.seekg( pos ) ;
  return true ;
}
/*===================================================================*/
void rescale_default_electrostatics_by_factor( string &line, const double &f ){
  rescale_default_hydrophobicity_by_factor( line, f ) ;
}
/*===================================================================*/
/*this version of rescale strips the link of the atom info*/
void rescale_barriers_of_link_by_factor( string &link, const double &factor ){
  string sbuf, *snumbers=NULL, sbuf2("") ;
  int n_numbers ;
  char cbuf[256] ;
  double e ;
  remove_trailing_comment( link ) ;
  remove_firsts( " ", link ) ;  remove_lasts( " ", link ) ;
  sbuf=link.substr(17,link.length( )-17) ;
  snumbers = split_nonDegenerated( " ", sbuf, n_numbers ) ;
  if( n_numbers > 2 ){/*There are internal barriers*/
    sbuf=link.substr( 0, 17 ) ;
    for( int i=2; i<n_numbers-1 ; i+=2 ){
      if( snumbers[i].find( s_inf ) == string::npos ){/*not infinite barrier*/
	string2char( snumbers[i], cbuf ) ; e = atof( cbuf ) ;
	e *= factor ; sprintf( cbuf, "%9.6f", e ) ; 
	snumbers[i].assign( cbuf ) ;
      }
    }
    for(int i=0; i<n_numbers; i++){ sbuf += "  " + snumbers[i]; }
    link = sbuf ;
  }
}
/*===================================================================*/
