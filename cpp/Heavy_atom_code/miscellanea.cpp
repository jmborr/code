using namespace std;
#include<iostream>
#include<fstream>
#include<string>
#include "pdbClasses2.h"
#include "miscellanea.h"
#include "bibliography.h"
/*===================================================================*/
string char2string( const char *the_char ){
  string tmp( the_char ) ;
  return tmp ;
}
/*===================================================================*/
void string2char( string &the_string, char *the_char ){
  the_string.copy( the_char, the_string.length( ), 0 ) ;
  the_char[ the_string.length( ) ] = 0 ;  
}
/*===================================================================*/
bool go_to( ifstream &IN, const string &start ) {
  string line ;
  do{/*go to linked pair section*/
    getline(IN, line) ;
  }while( line != start ) ;
  if( line == start ){ return true ; }
  return false ;
}
/*===================================================================*/
double *split_line2doubles( const string &line, int &n ){
  string sbuf = line, *stmp = NULL ;
  char div = ' ' ;
  remove_firsts( " ", sbuf ) ;  remove_lasts( " ", sbuf ) ;
  /*cout<<"sbuf=|"<<sbuf<<"|\n";*/
  stmp = split_nonDegenerated( &div , sbuf, n ) ;
  /*cout<<"n="<<n<<", sbuf[0]="<<sbuf[0]<<", sbuf["<<n-1<<"]="<<sbuf[n-1]<<endl;*/
  double *tmp = NULL ;
  if( stmp ){
    tmp = new double[ n ] ;
    for( int i=0; i<n; i++ ){ string2double( stmp[i], tmp[ i ] ) ; }
  }
  return tmp ;
}
/*===================================================================*/
int *split_line2ints( const string &line, int &n ){
  int *itmp = NULL ;
  double *dtmp = split_line2doubles( line, n ) ;
  if( dtmp ){
    itmp = new int[ n ] ;
    for( int i=0; i<n; i++ ){ itmp[ i ] = static_cast<int>( dtmp[ i ] ) ; }
  }
  return itmp ;
}
/*===================================================================*/
/*===================================================================*/
/*===================================================================*/

