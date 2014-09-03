#include<iostream.h>
using namespace std;
#include<string>
using std::string;
#include<fstream.h>
#include "bibliography.h"

/*===============================================================*/
bool is( string &line, string &the_class ){
  string query = "@" ; 
  query.append( the_class ) ;
  /*cout<<"query="<<query<<endl;*/
  if( line.find( the_class ) != string::npos ){ return true ; }
  return false ;
}
/*===============================================================*/
void remove_spaces( string &word){
  int x = word.find( " " ) ;
  while( x < string::npos ){
    word.replace( x, 1, "" ) ;
    x = word.find( " ", x ) ;
  }
  return;
}
/*===============================================================*/
void remove_last_colon( string &word){
  string last = word.substr( word.size( ) - 1, 1 ) ;
  if( last == "," ){ 
    string word2 = word.substr( 0, word.size( ) - 1 ) ;
    word = word2 ;
  }
}
/*===============================================================*/
int remove_lasts( const char *g, string &line ){
  string gS( g ), s ;
  int c = 0, begin = line.size( )-gS.size( ) ;
  do{
    s = line.substr( begin, line.size( ) ) ;
    if( s == gS ){ c++; line.erase( begin ) ; }
    else{ break ; }
    begin -= gS.size( ) ;
  }while(1) ;
  return c ;
}
/*===============================================================*/
int remove_firsts( const char *g, string &line ){
  string gS( g ), s ;
  int c = 0, begin = 0, gsize = gS.size( ) ;
  do{
    s = line.substr( 0, gsize ) ;
    if( s == gS ){ c++; line = line.substr( gsize, line.size( ) ) ; }
    else{ break ; }
  }while(1) ;
  return c ;
}
/*===============================================================*/
int replace( const char *oldC, const char *newC, string &line ){
  int c=0 ;
  string oldS( oldC ), newS( newC ) ;
  int s=oldS.size( ) ;
  int x = line.find( oldS ) ;
  while( x < string::npos ){
    line.replace( x, s, newS ) ;
    x = line.find( oldS , x ) ;
    c++ ;
  }
  return c ;
}
/*===============================================================*/
bool replace_last( const char *oldC, const char *newC, string &line ){
  string oldS( oldC ), newS( newC ) ;
  int s=oldS.size( ) ;
  int x = line.rfind( oldS ) ;
  line.replace( x, s, newS ) ;
  /*cout<<"line="<<line<<endl;*/
}
/*===============================================================*/
int remove( const char *garbage, string &line ){
  return replace( garbage, "", line ) ;
}
/*===============================================================*/
string char_to_string ( const char *line ){
  string tmp( line ) ;
  return tmp ;
}
/*===============================================================*/
int how_many( const char *div , string &line ){/*how many div's are there in line*/
  int c=0 ;
  string divS( div ) ;
  /*cout <<"divS="<<divS<<endl;*/
  int x = line.find( divS ) ;
  while( x < string::npos ){ 
    /*cout<<"x="<<x<<endl;*/
    c++ ; 
    x = line.find( divS , x+1 ) ; 
  }
  /*cout<<"c="<<c<<endl;*/
  return c ; 
}
/*===============================================================*/
int how_many_nonDegenerated( const char *div , string &line ){/*how many div's are there in line, but if two or more div's come together one after the other, then count them as only one instance*/
  int c=0 ;
  string divS( div ) ;  /*cout <<"divS="<<divS<<endl;*/
  int x = line.find( divS ), x_old=-1000000, l=divS.size( ) ;
  while( x < string::npos ){
    /*cout<<"x="<<x<<" ";*/
    if( x-x_old==l ){ x_old=x ; }
    if( x != x_old ){ 
      c++ ; x_old = x ; 
    }
    x = line.find( divS , x+1 ) ; 
  }
  /*cout<<"line=|"<<line<<"|, divS=|"<<divS<<"|, c="<<c<<endl;*/
  return c ; 
}
/*===============================================================*/
string * split( const char *div , string &line, int &n ){
  n =  how_many( div, line ) + 1; /*cout<<"n="<<n<<endl;*/
  string *tmp = NULL ;
  if( n>0 ){
    tmp = new string [ n ] ;
    string divS( div ) ;
    int c = 0 , x = 0, y = line.find( divS ), z ;
    while( y < string::npos ){
      z = y - x ;
      /*cout<<"c="<<c<<" x="<<x<<" y="<<y<<" z="<<z<<endl;*/
      tmp[c] = line.substr( x, z ) ;
      /*cout <<"c="<<c<<" tmp[0]="<<tmp[c]<<endl;*/
      c++ ;
      x = y + divS.size( ) ;
      y = line.find( divS, y+1 ) ;
    }
    z = line.size( ) - x ; tmp[c] = line.substr( x, z ) ;
  }
  return tmp ;
}
/*===============================================================*/
string * split_nonDegenerated( const char *div , string &line, int &n ){
  n =  how_many_nonDegenerated( div, line ) + 1;   /*cout<<"n="<<n<<endl;*/
  string *tmp = NULL ;
  if( n>0 ){
    tmp = new string [ n ] ;
    string divS( div ) ; /*cout<<"divS=|"<<divS<<"|\n" ;*/
    int c = 0 , x = 0, l= divS.size( ), y = line.find( divS ), z ;
    while( y < string::npos ){
      z = y - x ;
      if(z>0){
	/*cout<<"c="<<c<<" x="<<x<<" y="<<y<<" z="<<z<<endl;*/
	tmp[c] = line.substr( x, z ) ;
	/*cout <<"tmp["<<c<<"]=|"<<tmp[c]<<"|\n";*/
	c++ ; y += l ;
      }
      else{ y++ ; }
      x = y ;  y = line.find( divS, x ) ;
    }
    z = line.size( ) - x ; tmp[c] = line.substr( x, z ) ;
  }
  /*for(int i=0;i<n;i++){cout<<"tmp["<<i<<"]="<<tmp[i]<<endl;}*/
  return tmp ;
}
/*===============================================================*/
bool last_name_first( string &name ){
  if( how_many( ",", name ) == 0 ){
    string first_name, family_name, line ;
    int x = name.rfind( "." ) ;

    first_name = name.substr( 0, x+1 ) ;
    remove_lasts( "~", first_name ) ;

    family_name = name.substr( x+1, name.size( ) ) ;
    remove_spaces( family_name ) ;
    remove_firsts( "~", family_name ) ;

    name = family_name ; name.append( ", " ) ;
    name.append( first_name ) ;
    return true ;
  }
  return false ;
}
/*===============================================================*/
bool leave_only_last_name( string &name ){
  string first_name, family_name, line ;
  if( how_many( ",", name ) == 0 ){
    int x = name.rfind( "." ) ;
    family_name = name.substr( x+1, name.size( ) ) ;
    remove_spaces( family_name ) ;
    remove_firsts( "~", family_name ) ;
    name = family_name ;
    return true ;
  }
  else if( how_many( ",", name ) == 1 ){
    int x = name.find( "," ) ;
    family_name = name.substr( 0, x ) ;
    remove_spaces( family_name ) ;
    name = family_name ;
    return true ;
  }
  return false ;
}
/*===============================================================*/
/*returns true if query found, leave the pointer in at the 
  beggining of query. If not found, returns false and does not
  change the position of the pointer
*/
bool find_next( const char *query, ifstream &in ){
  int n=0, i=1, j ;
  int p = in.tellg( ) ;
  char c ;
  bool flag = false ;
  /*int junk=0 ;*/
  while( (c = in.get( )) != EOF ){
    /*while( query[junk] != '\0' ){
      cout << query[junk] <<endl;
      junk++ ;
      }*/
    /*printf("%c",c );*/
    if( c == query[0] ){/*the first character matches*/
      /*cout <<"found "<<c<<endl;*/
      flag = true ;
      j = in.tellg( ) ;
      while( query[i] != '\0' ){/*compare next characters*/
	c = in.get( ) ;
	/*cout <<"c="<<c<<" query["<<i<<"]="<<query[i]<<"--";*/
	if( c != query[i] ){
	  /*cout<<endl;*/
	  flag = false ; 
	  i = 1 ;
	  break ; 
	}
	i++ ;
      }
      if( flag == true ){
	in.seekg( j - 1 ) ;
	return true ;	
      }
    }      
  }
  if( in.eof() ){ in.clear( ios::goodbit ) ; }
  in.seekg( p ) ;/*return to origin if query not found*/
  return false ;
}
/*===============================================================*/
int find_next2( const char *query, FILE *in ){
  int i = 0, x = 0 ;
  long p = ftell( in ) ;
  char c ;
  bool flag = false ;
  while( (c = fgetc( in )) != EOF ){
    x++ ;
    if( c == query[0] ){
      flag = true ; i++ ;
      while( query[i] != '\0' ){
	c = fgetc( in ) ;
	x++ ; 
	if( c != query[i] ){
	  i = 0 ;
	  flag = false ;
	  break ;
	}
	i++ ;
      }
    }
    if( flag == true ){
      fseek( in, p, SEEK_SET ) ;
      return x ;
    }
  }
  if( feof( in ) ){ clearerr( in ) ; fseek( in, p, SEEK_SET ) ; }
  return 0 ;
}
/*===============================================================*/
int my_getc( const char *query, char *buff, FILE *in ){
  int x = find_next2( query, in ) - 1;
  if( x ){
    for( int i=0; i<x; i++ ){
      buff[ i ] = getc( in ) ;
    }
    buff[ x ] = '\0' ;
  }
  return x ;  
}
/*===============================================================*/
bool print_till_next( const char *query, FILE *in,
		      ofstream &out ){
  /*cout<<"FUNCTION bool print_till_next( const char *, ifstream &, ofstream & )\n";*/
  char c; char *d = &c + 1 ;
  char buff[1000] ;
  int i=0 ;
  long j, p = ftell( in ) ;
  bool flag = false ;
  while( (c = fgetc( in )) != EOF ){
    /*cout << c ;*/
    if( c != query[0] ){ 
      out << c ; 
    }
    else{
      /*cout <<"found "<<c<<endl;*/
      flag = true ;  buff[i] = c ;  i++ ;  j = ftell( in ) ;
      while( query[i] != '\0' ){/*compare next characters*/
	c = fgetc( in ) ;  buff[i] = c ;
	/*cout<<c;*/
	if( c != query[i] ){
	  for( int k=0; k<=i; k++ ){ out << buff[k] ; }
	  i = 0 ;  flag = false ;  break ;
	}
	i++ ;
      }
    }
    if( flag ){  fseek( in, j - 1, SEEK_SET ) ;  return true ;  }
  }
  if( feof( in ) ){ clearerr( in ) ; fseek( in, p, SEEK_SET ) ; }
  return false ;
}
/*===============================================================*/
string *get_handles_from_cite_command( string &line, int &nh ){
  string *tmp = NULL ;
  nh = 0 ;
  if( line.find("\\cite{") != string::npos ){
    int begin = line.find("\\cite{") + 6 ;
    /*cout<<"begin="<<begin<<endl;*/
    int lapse = line.find("}") - begin ;
    line = line.substr( begin, lapse ) ;
    tmp = split( ",", line, nh ) ;
    /*cout<<"nh="<<nh<<endl;*/
    for( int i=0; i<nh; i++ ){ 
      remove_spaces( tmp[i] ) ;
      /*cout << tmp[i] <<endl;*/
    }
    return tmp ;
  }
}
/*===============================================================*/
string *create_list_of_handles( int &n, ifstream &paper ){
  string *tmp = NULL ;
  list_alias all_handles ;
  all_handles.create_list_of_handles( paper ) ;
  /*cout << all_handles ;*/
  return all_handles.output_aliases( n );
}
/*===============================================================*/
string *sort_by_last_name_year_last_names( string *handles ,
					   const int &nh,
					   ifstream &bibl ){
  int n = nh ;
  string *tmp = handles ;
  list_entry citations ;
  citations.fill_with_list_handles( handles, nh, bibl ) ;
  citations.buble_sort_by_last_name_year_last_names( ) ;
  return citations.output_handles( n ) ;
}
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/


/*===============  C L A S S   A L I A S    =====================*/
void parse_alias( const string &line, string &alias, string &full ){
  int begin = line.find("{")+1, lapse = line.find("=") -  begin ;
  alias = line.substr( begin, lapse ) ;
  /*cout<<"alias="<<alias<<"|"<<endl;*/
  remove_spaces( alias ) ;
  /*cout<<"alias="<<alias<<"|"<<endl;*/
  begin = line.find("\"") + 1 ; lapse = line.find("\"", begin ) -  begin ;
  full = line.substr( begin, lapse ) ;
  /*cout<<"full="<<full<<"|"<<endl;*/
}
/*===============================================================*/
bool is_alias( const string &line ){
  if( line.find("@string") != string::npos ){ return true ; }
  return false ;
}
/*===============================================================*/
alias::alias( ){ a = "" ; full = "" ; }
/*===============================================================*/
alias::alias( const string &line){
  /*cout<<"line="<<line<<endl;*/
  if( line.find("string") ){
    string the_alias, the_full ;
    parse_alias( line, the_alias, the_full) ;
    a = the_alias ; full = the_full ;
  }
  else{ a = ""; full = ""; }
}
/*===============================================================*/
void alias::fill_members( string &line ){
  if( line.find("@string") ){
    string the_alias, the_full ;
    parse_alias( line, the_alias, the_full) ;
    a = the_alias ; full = the_full ;
  }
}
/*===============================================================*/
string alias::print( ){
  string tmp ="@string {" ;  
  tmp.append(a);
  tmp.append("=\"");
  tmp.append(full);
  tmp.append("\"}");
  return tmp;
}
/*===============================================================*/
void alias::assign( const string &the_alias, const string &the_full ){
  a = the_alias ; full = the_full ;
}
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/



/*===============  C L A S S   N O D E _ A L I A S  =================*/
node_alias::node_alias( const alias &input ){
  the_alias = input ; next = NULL ;
}
/*===============================================================*/



/*===============  C L A S S   L I S T _ A L I A S ================*/
ifstream &operator>>(ifstream &in, list_alias &receiver ){
  string line ;
  /*record actual position of file listOfContacts*/
  int position=in.tellg( ) ;
  while( getline( in, line ) ){
    /*cout <<"line="<<line<<endl;*/
    if( is_alias( line ) ){ 
      /*cout<<"success!\n";*/
      receiver.insertAtBack( line ) ; 
    }
  }
  /*rewind to original position, so that file can be used more than once*/
  if( in.eof() ){ in.clear( ios::goodbit ) ; }
  in.seekg( 0 ) ;
  /*getline( in, line ) ; cout<<"line="<<line<<"|\n";*/
  return in ;
}
/*===============================================================*/
ostream &operator<<(ostream &out, list_alias &output ){
  if( !output.isEmpty( ) ){
    string line ;
    node_alias *currentPtr = output.firstPtr ;
    while( currentPtr != NULL){
      line = currentPtr->the_alias.print( ) ;
      out<<line<<endl;
      currentPtr = currentPtr->next ;
    }
  }
  return out ;
}
/*===============================================================*/
list_alias::list_alias( ){
  firstPtr = lastPtr = NULL ;
}
/*===============================================================*/
list_alias::~list_alias( ){
  /*cout<<"FUNCTION list_alias::~list_alias( )\n";*/
  /*cout <<*this ;*/
  this->empty( ) ; 
}
/*===============================================================*/
node_alias *list_alias::getNewNode( const alias &input ){
  node_alias *ptr = new node_alias( input ) ;
  assert( ptr != 0 ) ;
  return ptr ;
}
/*===============================================================*/
bool list_alias::isEmpty( ){ return firstPtr == NULL ; }
/*===============================================================*/
void list_alias::rmFromFront( ){
  node_alias *currentPtr = firstPtr ;
  firstPtr = ( *currentPtr ).next ; delete currentPtr ;
}
/*===============================================================*/
void list_alias::empty( ){
  if( !isEmpty( ) ){
    while( firstPtr != NULL ){ this->rmFromFront( ); }
  }
}
/*===============================================================*/
void list_alias::insertAtFront( const alias &input ){
  node_alias *newPtr = getNewNode( input ) ;
  if( isEmpty() ){ firstPtr = lastPtr = newPtr ;}
  else{ ( *newPtr ).next = firstPtr  ; firstPtr = newPtr ; }
}
/*===============================================================*/
void list_alias::insertAtBack( const alias &input ){
  /*cout<<"FUNCTION void list_alias::insertAtBack( const alias & )\n";*/
  node_alias *newPtr = getNewNode( input ) ;
  if( isEmpty( ) ){ firstPtr = lastPtr = newPtr ; }
  else{ lastPtr->next = newPtr ; lastPtr = newPtr ; }
  /*cout<<"FINISHED void list_alias::insertAtBack( const alias & )\n";*/
}
/*===============================================================*/
void list_alias::insertAtBack( const string &line ){
  alias input( line ) ;
  node_alias *newPtr = getNewNode( input ) ;
  if( isEmpty() ){ firstPtr = lastPtr = newPtr ; }
  else{ ( *lastPtr ).next = newPtr ; lastPtr = newPtr ; }

}
/*===============================================================*/
int list_alias::size( ){
  int size = 0 ;
  if( isEmpty() ){ size = 0 ; }
  else{
    node_alias *currentPtr = firstPtr ;
    while( currentPtr != NULL ){
      size++ ;
      currentPtr = currentPtr->next ;
    }
  }
  return size ;
}
/*===============================================================*/
int list_alias::create_list_of_handles( ifstream &paper ){
  /*cout<<"int list_alias::create_list_of_handles( ifstream &)\n" ;*/
  int p = paper.tellg( ) ;
  if( !(this->isEmpty( )) ){ this->empty( ) ; }
  int nh  ;
  int position = paper.tellg( ) ;
  string line, empty_str( "" ) ;
  char cline[1000] ;
  string *tmp ;
  alias new_alias ;
  node_alias *currentPtr ;
  while( find_next("\\cite", paper) ){
    /*getline( paper, line) ; cout<<"line="<<line<<endl;*/
    paper.get( cline, 999, '}' ) ;
    line.assign( cline ) ; line.append( "}" ) ;
    /*cout<<line<<endl;*/
    tmp = get_handles_from_cite_command( line, nh );
    for( int i=0; i<nh; i++ ){
      if( !(this->is_there_alias(tmp[i])) ){
	new_alias.assign( tmp[i], empty_str ) ;
	this->insertAtBack( new_alias ) ;
      }
    }
  }
  paper.clear( ios::goodbit ) ; paper.seekg( p ) ; 
  paper.seekg( p ) ;
  return this->size( );
}
/*===============================================================*/
bool list_alias::is_there_alias( const string &line ){
  node_alias *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    if( currentPtr->the_alias.a == line ){ return true ; }
    currentPtr = currentPtr->next ;
  }
  return false ;
}
/*===============================================================*/
string *list_alias::output_aliases( int &n ){
  string *tmp = NULL ;
  n = this->size( ) ;
  /*cout <<"n="<<n<<endl;*/
  if( n > 0 ){
    tmp = new string[ n ] ;
    int i = 0 ;
    node_alias *currentPtr = firstPtr ;
    while( currentPtr != NULL ){
      tmp[i] = currentPtr->the_alias.a ;  i++ ;
      /*cout <<"i="<<i<<endl;*/
      currentPtr = currentPtr->next ;
    }
  }
  return tmp ;
}
/*===============================================================*/
string list_alias::get_the_full( const string &the_a ){
  /*cout<<"FUNCTION string list_alias::get_the_full( const string & )\n" ;*/
  string tmp = "" ;
  /*cout<<"the_a="<<the_a<<endl;*/
  /*cout<<*this<<endl;*/
  node_alias *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    if( currentPtr->the_alias.a == the_a ){
      tmp = currentPtr->the_alias.full ;
      return tmp ;
    }
    currentPtr = currentPtr->next ;
  }
  return tmp ;
}
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/



/*===============  C L A S S   F I E L D    =================*/
bool is_value_in_brakets( const string &line ){
  int position = line.find("{") ;
  if( position < line.find("\"") && position < line.size( ) ){
    return true ;
  }
  return false;
}
/*===============================================================*/
string in_brakets( const string &line ){
  int begin = line.find("{") + 1 ;
  int lapse = line.rfind("}") - begin ;
  return line.substr( begin, lapse ) ;
}
/*===============================================================*/
bool is_value_in_quotes( const string &line ){
  int position = line.find("\"") ;
  if( position < line.find("{") && position < line.size( ) ){
    return true ;
  }
  return false;
}
/*===============================================================*/
string in_quotes( const string &line ){
  int begin = line.find("\"") + 1 ;
  int lapse = line.rfind("\"") - begin ;
  return line.substr( begin, lapse ) ;
}
/*===============================================================*/
string in_naked( const string &line ){
  int begin = line.find("=") + 1 ;
  int lapse = line.size( ) - begin ;
  string tmp = line.substr( begin, lapse ) ;
  remove_spaces( tmp ) ;
  remove_last_colon( tmp ) ;
  /*cout <<"tmp="<<tmp<<"|\n" ;*/
  return tmp ;
}
/*===============================================================*/
field_kw determine_field_kw( const string &line ){
  int end = line.find("=") ; 
  string keyword = line.substr( 0, end ) ;
  remove_spaces( keyword ) ;
  /*cout<<"keyword="<<keyword<<endl;*/
  field_kw value=_NON_FIELD_ ;
  for( int i=0; i<nun_field_keywords; i++ ){
    if( keyword == field_keywords[i] ){
      value = static_cast<field_kw>( i ) ;
      break ;
    }
  }
  return value ;
}
/*===============================================================*/
bool is_field( const string &line ){
  /*cout <<"line="<<line<<endl;*/
  if( static_cast<int>( determine_field_kw( line ) ) ){ return true ; }
  return false ;
}
/*===============================================================*/
string determine_field_value( const string &line ){
  /*cout<<"THIS IS FUNCTION:string determine_field_value( const string & )\n";*/
  int index = static_cast<int>( determine_field_kw(line) ) ;
  /*cout <<"index="<<index<<endl;*/
  if( index ){
    if( is_value_in_brakets(line) ){ 
      /*cout<<"in brakets\n";*/
      return in_brakets(line) ; 
    }
    else if( is_value_in_quotes(line) ){ 
      /*cout<<"in quotes\n";*/
      return in_quotes(line) ; 
    }
    else{ 
      /*cout<<"naked line="<<line<<"\n";*/
      return in_naked(line) ; 
    }
  }
  return "" ;
}
/*===============================================================*/
ostream &operator<<(ostream &out, const field &input ){
  field tmp = input ;
  out << tmp.print( ) ;
  return out ;
}
/*===============================================================*/
ifstream &operator>>(ifstream &in, field &receiver){
  string line ;
  getline( in, line ) ;
  receiver.assign( line ) ;
  return in ;
}
/*===============================================================*/
field::field( ){
  kw = _NON_FIELD_ ;
  value = "" ;
}
/*===============================================================*/
field::field( const string &line ){
  if( is_field( line ) ){
    kw =  determine_field_kw( line ) ;
    value = determine_field_value( line ) ;
    /*cout<<"kw="<<field_keywords[static_cast<int>(kw)]<<endl;
      cout <<"value="<<value<<"|\n";*/
  }
  else{
    kw = _NON_FIELD_ ;
    value = "" ;
  }
}
/*===============================================================*/
bool field::assign( const string &line ){
  if( is_field( line ) ){
    kw =  determine_field_kw( line ) ;
    value = determine_field_value( line ) ;
    /*cout <<"value="<<value<<endl;*/
    return true;
  }
  return false;
}
/*===============================================================*/
string field::get_value( ){ return value ; }
/*===============================================================*/
field_kw field::get_keyword2( ){ return kw ; }
/*===============================================================*/
int  field::get_keyword3( ){ return static_cast<int>( kw ) ; }
/*===============================================================*/
string field::get_keyword( ){
  int index = static_cast<int>( kw ) ;
  /*cout <<"index="<<index<<endl;  exit(1);*/
  return field_keywords[ index ] ;
}
/*===============================================================*/
string field::print( ){
  int index = static_cast<int>( kw ) ;
  string tmp = "  " ;
  tmp.append( field_keywords[ index ] ) ;
  /*cout<<"kw="<<field_keywords[ index ]<<endl;*/
  tmp.append(" = {" ) ;
  tmp.append( value ) ;
  tmp.append( "}" ) ;
  /*cout<<"tmp="<<tmp<<endl;*/
  return tmp ;
}
/*===============================================================*/
bool field::last_names_first( ){
  if( kw == _AUTH_ ){
    int n = 0 ;
    string line = this->get_value( ) ;
    string *names = split( " and ", line, n ) ;
    value = "" ;
    for( int i=0; i<n; i++ ){
      last_name_first( names[i] ) ;
      value.append( names[i] ) ;
      value.append( " and " ) ;
    }
    remove_lasts( " and ", value ) ;
  }
  return false ;
}
/*===============================================================*/
string *field::last_names( int &n ){
  string *names = NULL ;
  if( kw == _AUTH_ ){
    n = 0 ;
    string line = this->get_value( ) ;
    names = split( " and ", line, n );
    for( int i=0; i<n; i++ ){
      leave_only_last_name(  names[i] ) ;
    }
  }
  return names ;
}
/*===============================================================*/
bool field::format_author_for_BJ( ){
  if( kw == _AUTH_ ){
    int n = 0 ;
    string *names = split( " and ", value, n ) ;
    /*cout<<"n="<<n<<endl;*/
    /*cout<<"names[0]="<<names[0]<<endl;*/
    last_name_first( names[0] ) ;
    value = names[ 0 ] ;
    /*cout<<"value="<<value<<endl;*/
    for( int i=1; i<n; i++ ){
      value.append( ", " ) ;
      value.append( names[ i ] ) ;
    }
    /*cout<<"value="<<value<<endl;*/
    if( n>1 ){ replace_last( ", ", " and ", value ) ; }
    return true ;
  }
  return false ;
}
/*===============================================================*/
bool field::unalias( const list_alias &all_alias ){
  node_alias *currentPtr = all_alias.firstPtr ;
  while( currentPtr != NULL ){
    /*cout<<"value="<<value<<" alias="<<currentPtr->the_alias.a<<endl;*/
    if( value == currentPtr->the_alias.a ){ 
      value = currentPtr->the_alias.full ;
      return true ;
    }
    currentPtr = currentPtr->next ;
  }
  return false ;
}
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/


/*===============  C L A S S   N O D E _ F I E L D    =================*/
node_field::node_field( const field &input ){
  the_field = input ;  next = NULL ;
}
/*===============================================================*/


/*===============  C L A S S   E N T R Y    =================*/
ifstream &operator>>( ifstream &in, entry &receiver ){
  if( !receiver.isEmpty( ) ){ receiver.empty( ) ; }
  string line ;
  do{ getline( in, line ) ; }while( !is_entry( line ) ) ;
  if( is_entry( line ) ){/*in case we previously reach eof*/
    if( has_handle( line ) ){
      receiver.kw = determine_entry( line ) ;
      /*cout <<"keyword="<<entry_keywords[static_cast<int>(receiver.kw)]<<endl;*/
      receiver.handle = determine_handle( line ) ;
      /*cout<<"handle="<<determine_handle( line )<<"|\n";*/
      while( getline( in, line ) ){
	/*cout<<"line="<<line<<"|\n";*/
	if( is_end_of_entry( line ) ){ break ; }
	receiver.insertAtBack( line ) ;
      }
    }
  }
  return in ;
}
/*===============================================================*/
ostream &operator<<( ostream &out, entry &output ){
  /*cout<<"FUNCTION ostream &operator<<(ostream &,entry &)\n";*/
  if( !output.isEmpty( ) ){
    int index = static_cast<int>( output.kw ) ;
    out <<"@"<<entry_keywords[ index ]<<"{"<<output.handle<<",\n" ;
    string tmp ;
    node_field *currentPtr = output.firstPtr ;
    while( currentPtr != NULL){
      tmp = currentPtr->the_field.print( ) ;
      out << tmp ;
      currentPtr = currentPtr->next ;
      if( currentPtr != NULL ){ out << ",\n" ; }
    }
    out << "\n}\n" ;
  }
  /*cout<<"FINISHED: ostream &operator << ( ostream &, entry & )\n";*/
  return out ;
}
/*===============================================================*/
entry_kw determine_entry( const string &line ){
  int begin = line.find("@") + 1 ; 
  int lapse = line.find("{") - begin ;
  string keyword = line.substr( begin, lapse ) ;
  remove_spaces( keyword ) ;
  /*cout <<"keyword="<<keyword<<endl;*/

  entry_kw value=_NON_ENTRY_ ;
  for( int i=0; i<nun_entry_keywords; i++ ){
    if( keyword == entry_keywords[i] ){
      value = static_cast<entry_kw>( i ) ;
      break ;
    }
  }
  return value ;
}
/*===============================================================*/
string determine_handle( const string &line ){
  string handle = "" ;
  if( is_entry( line ) ){
    int begin = line.find( "{" ) + 1 ;
    handle = line.substr( begin, line.size() - begin - 1 ) ;
  }
  return handle ;
}
/*===============================================================*/
bool is_entry( const string &line ){
  entry_kw kw = determine_entry( line );
  if( kw != _NON_ENTRY_ ){ return true ; }
  return false ;
}
/*===============================================================*/
bool is_end_of_entry( const string &line ){
  string tmp = line ;
  remove_spaces( tmp ) ;
  if( tmp == "}" ){ return true ; }
  return false ;
}
/*===============================================================*/
bool has_handle( const string &line ){
  if( is_entry( line ) ){
    int begin = line.find( "{" ) + 1 ;
    if( line.size( ) - begin > 1 ){ return true ; }
  }
  return false ;
}
/*===============================================================*/
entry::entry( ){
  firstPtr = lastPtr = NULL ; 
  kw = _NON_ENTRY_ ;
  handle = "" ;
}
/*===============================================================*/
entry::~entry( ){ 
  this->empty( ) ; 
  kw = _NON_ENTRY_ ;
  handle = "" ;
  /*cout<<"entry::~entry finished\n" ;*/
}
/*===============================================================*/
void entry::empty( ){
 if( !isEmpty( ) ){ while( firstPtr != NULL ){ this ->rmFromFront( ); } }
 kw = _NON_ENTRY_ ;
 handle.assign( "" ) ;
}
/*===============================================================*/
node_field *entry::getNewNode( const field &input ){
  node_field *ptr = new node_field( input ) ;
  assert( ptr != 0 ) ;
  return ptr ;
}
/*===============================================================*/
bool entry::isEmpty( ){ return firstPtr == NULL ; }
/*===============================================================*/
void entry::rmFromFront( ){
  node_field *currentPtr = firstPtr ;
  firstPtr = ( *currentPtr ).next ; delete currentPtr ;
}
/*===============================================================*/
void entry::insertAtBack( const field &input ){
  node_field *newPtr = getNewNode( input ) ;
  if( isEmpty() ){ firstPtr = lastPtr = newPtr ; }
  else{ ( *lastPtr ).next = newPtr ; lastPtr = newPtr ; }
}
/*===============================================================*/
void entry::insertAtBack( const string  &line ){
  /*cout <<"line="<<line<<endl;*/
  field item( line ) ;
  /*string tmp = item.print( ) ;   cout<<"item="<<tmp<<endl;*/
  this->insertAtBack( item ) ;
}
/*===============================================================*/
bool entry::fill( ifstream &in, const string &the_handle ){
  /*cout<<"bool entry::fill( ifstream &, const string &)\n";*/
  string line, query ;
  int position, old_position= in.tellg( ) ;
  do{
    position = in.tellg( ) ;
    getline( in, line ) ; 
    if( is_entry( line ) && has_handle( line ) ){
      query = determine_handle( line ) ;
      /*cout <<"query="<<query<<"|\n" ;*/
      if ( the_handle == query ){
	if( !this->isEmpty( ) ){ this->empty( ) ; }
	in.seekg( position ) ;	
	/*getline( in, line ) ; cout <<"line="<<line<<"|\n";*/
	in >> *this ;
	/*cout <<"this="<< *this ;*/
	/*rewind to original position*/
	in.clear( ios::goodbit ) ; in.seekg( old_position ) ; 
	in.seekg( old_position ) ;
	return true ;
      } 
    }
  }while( !in.eof( ) ) ;
  in.clear( ios::goodbit ) ; in.seekg( old_position ) ;
}
/*===============================================================*/
bool entry::is_there_field( const field_kw &query ){
  node_field *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    if( currentPtr->the_field.kw == query ){ return true ; }
    currentPtr = currentPtr->next ;
  }
  return false ;
}
/*===============================================================*/
bool entry::get_field( const field_kw &query, field &the_field ){
  node_field *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    if( currentPtr->the_field.kw == query ){
      the_field = currentPtr->the_field ;
      return true ;
    }
    currentPtr = currentPtr->next ;
  }
  return false ;
}
/*===============================================================*/
string entry::get_field( const field_kw &query ){
  field the_field ;
  this->get_field( query, the_field ) ;
  return the_field.get_value( ) ;
}
/*===============================================================*/
string entry::get_handle( ){ return handle ; }
/*===============================================================*/
bool entry::last_names_first( ){
  if( this->is_there_field( _AUTH_ ) ){
    node_field *currentPtr = firstPtr ;
    while( currentPtr != NULL ){
      if( currentPtr->the_field.kw == _AUTH_ ){
	return currentPtr->the_field.last_names_first( ) ;
      }
      currentPtr = currentPtr->next ;
    }
  }
  return false ;
}
/*===============================================================*/
int entry::unalias( const list_alias &all_alias ){
  int n = 0 ;
  node_field *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    if( currentPtr->the_field.unalias( all_alias ) ){ n++ ; }
    currentPtr = currentPtr->next ;
  }
  return n ;
}
/*===============================================================*/
bool entry::format_for_BJ( ){
  node_field *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    if( currentPtr->the_field.kw == _AUTH_ ){
      return currentPtr->the_field.format_author_for_BJ( ) ;
    }
    currentPtr = currentPtr->next ;
  }
  return false ;
}
/*===============================================================*/
entry &entry::operator=( entry &right_side ){
  if( ! this->isEmpty( ) ){  this->empty( ) ; }
  if( ( &right_side != this )  && !( right_side.isEmpty() ) ) {
    kw = right_side.kw ;
    handle = right_side.handle ;
    firstPtr = getNewNode( right_side.firstPtr->the_field ) ;
    node_field *currentRightPtr = right_side.firstPtr->next ;
    node_field *currentPtr = firstPtr ;
    while( currentRightPtr != NULL ){
      currentPtr->next = getNewNode( currentRightPtr->the_field ) ;
      currentPtr = currentPtr->next ;
      currentRightPtr = currentRightPtr->next ;
    }
  }
  return *this ;
}
/*===============================================================*/
string *entry::last_names( int &n ){
  string *tmp = NULL ;
  if( this->is_there_field( _AUTH_ ) ){
    node_field *currentPtr = firstPtr ;
    while( currentPtr != NULL ){
      if( currentPtr->the_field.kw == _AUTH_ ){
	tmp = currentPtr->the_field.last_names( n ) ;
	break ;
      }
      currentPtr = currentPtr->next ;
    }
  }
  return tmp ;
}
/*===============================================================*/
string *entry::last_name_year_last_names( int &n ){
  string *tmp = NULL, *last_names = NULL ;
  field the_authors, the_year ;
  if( this->get_field( _AUTH_, the_authors) &&
      this->get_field( _YEAR_, the_year )      ){
    last_names = the_authors.last_names( n ) ;
    n++; tmp = new string[ n ] ;
    tmp[0] = last_names[0];
    tmp[1] = the_year.get_value( ) ;
    for( int i=2; i<n; i++ ){ tmp[i] = last_names[ i-1 ] ; }
  }
  return tmp ;
  }
/*===============================================================*/
int entry::compare_by_last_name_year_last_names( entry &right ){
  int c = 0 ;
  field the_authors, the_authors2, the_year, the_year2 ;
  if ( this->get_field( _AUTH_, the_authors)  &&
       this->get_field( _YEAR_, the_year)     &&
       right.get_field( _AUTH_, the_authors2) &&
       right.get_field( _YEAR_, the_year2)      ){
    int n1, n2, n ;
    string *s1 = this->last_name_year_last_names( n1 ) ;
    string *s2 = right.last_name_year_last_names( n2 ) ;
    if( n1 < n2 ){ n = n1 ; }
    else{ n = n2 ; }
    string str1, str2 ;
    for( int i=0; i<n; i++ ){
      str1 = s1[i], str2 = s2[i] ;
      remove_firsts ( "{", str1 ) ;  remove_lasts( "}", str1 ) ;
      remove_firsts ( "{", str2 ) ;  remove_lasts( "}", str2 ) ;
      if( str1 < str2 ){ c = -1 ; break ; }
      else if( str1 > str2 ){ c = 1 ; break ; }
    }
    if( c==0 ){
      if( n1 < n2 ){ c=-1 ; }
      else if( n1 > n2 ){ c=1 ; }
    }
  }
  return c ;
}
/*===============================================================*/
alias entry::BJ_style_handle( ){
  alias pair ;
  string a = this->get_handle( ), full = "" ;
  int nh ;
  string *the_list = this->last_name_year_last_names( nh ) ;
  for( int i=0; i<nh; i++ ){
    remove_firsts( "{", the_list[i] ) ;
    remove_lasts( "}", the_list[i] ) ;
  }
  full.append( the_list[0] ) ;
  if( nh == 2 ){
    full.append( ", " ) ;
    full.append( the_list[1] ) ; 
  }
  else if( nh == 3 ){
    full.append( " \\& " ) ; full.append( the_list[2] ) ;
    full.append( ", " ) ; full.append( the_list[1] ) ;
  }
  else{
    full.append( " \\emph{et al.}, " ) ;
    full.append( the_list[1] ) ;
  }
  pair.assign( a, full ) ;
  return pair ;
}
/*===============================================================*/
void entry::print_BJ_style( ofstream &out ){
  /*cout<<*this<<endl;*/
  field f_author ; 
  this->get_field( _AUTH_, f_author ) ;
  /*cout<<"author="<<f_author<<endl;*/
  f_author.format_author_for_BJ( ) ;
  string s_author, s_title, s_year ;
  s_author = f_author.get_value( ) ;
  s_title = this->get_field( _TITLE_ ) ;
  s_year = this->get_field( _YEAR_ ) ;
  if( kw == _BOOK_ ){
    string publisher = this->get_field( _PUBLISHER_ ) ;
    string address = this->get_field( _ADDR_ ) ;
    out<<s_author<<" "<<s_year<<", "<<s_title<<". "
	<<publisher<<", "<<address<<"\n\n" ;
  }
  else if( kw == _ARTICLE_ ){
    string  s_journal = this->get_field( _JOURN_ ) ;
    string  s_volume = this->get_field( _VOL_ ) ;
    string  s_pages= this->get_field( _PAGES_ ) ;
    out<<s_author<<" "<<s_year<<", "<<s_title<<". \\emph{"
	<<s_journal<<"} "<<s_volume<<":"<<s_pages<<"\n\n" ;
  }
}
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/



/*===============  C L A S S   N O D E _ E N T R Y    =================*/
node_entry::node_entry( entry &input ){
  the_entry = input ;  next = NULL ;  prev = NULL ;
}
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/



/*===============  C L A S S   L I S T _ E N T R Y    ===========*/

ostream &operator<<( ostream &out, list_entry &output ){
  node_entry *currentPtr = output.firstPtr ;
  while( currentPtr != NULL ){
    out << currentPtr->the_entry ;
    currentPtr = currentPtr->next ;
    /*cout<<"currentPtr="<<currentPtr<<endl;*/
  }
  return out ;
}
/*===============================================================*/
node_entry *list_entry::getNewNode( entry &input ){
  node_entry *ptr = new node_entry( input ) ;
  assert( ptr != 0 ) ;
  return ptr ;
}
/*===============================================================*/
list_entry::list_entry( ){
  firstPtr = lastPtr = NULL ;
}
/*===============================================================*/
list_entry::~list_entry( ){ 
  /*cout<<"FUNCTION list_entry::~list_entry( )\n";*/
  this->empty( ) ; 
}
/*===============================================================*/
void list_entry::empty( ){
  /*cout<<"FUNCTION void list_entry::empty( )\n";*/
  if( !isEmpty( ) ){
    while( firstPtr != NULL ){ 
      /*cout<<"size="<<this->size( )<<endl;*/
      this ->rmFromFront( ); 
      /*cout<<"firstPtr="<<firstPtr<<endl;*/
    }
  }
}
/*===============================================================*/
bool list_entry::isEmpty( ){ return firstPtr == NULL ; }
/*===============================================================*/
void list_entry::rmFromFront( ){
  node_entry  *currentPtr = firstPtr ;
  firstPtr = currentPtr->next ; 
  if( firstPtr != NULL ){ firstPtr->prev = NULL ; }
  delete currentPtr ;
}
/*===============================================================*/
void list_entry::insertAtBack( entry &input ){
  /*cout<<"void list_entry::insertAtBack( entry & )";*/
  node_entry  *newPtr = getNewNode( input ) ;
  if( isEmpty() ){ firstPtr = lastPtr = newPtr ; }
  else{ 
    newPtr->prev = lastPtr ;
    lastPtr->next = newPtr ;     
    lastPtr = newPtr ; 
  }
}
/*===============================================================*/
void list_entry::fill_with_list_handles( string *list, const int &n, 
					 ifstream &bibl ){
  /*cout<<"void list_entry::fill_with_list_handles( string *, const int &, ifstream &bibl )\n";*/
  entry citation ;
  for( int i=0; i<n; i++){
    citation.fill( bibl, list[i] ) ;
    /*cout<<"list["<<i<<"]"<<"="<<list[i]<<endl;*/
    /*cout<<"citation=\n"<<citation<<endl ;*/
    this->insertAtBack( citation ) ;
  }
}
/*===============================================================*/
int list_entry::size( ){
  int size = 0 ;
  node_entry  *currentPtr = firstPtr ;
  while( currentPtr != NULL ){ 
    size++ ; 
    currentPtr = currentPtr->next ;
  }
  return size ;
}
/*===============================================================*/
int list_entry::unalias( const list_alias &all_alias ){
  int n = 0 ;
  node_entry *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    n += currentPtr->the_entry.unalias( all_alias ) ;
    currentPtr = currentPtr->next ;
  }
  return n ;
}
/*===============================================================*/
void list_entry::buble_sort_by_last_name_year_last_names( ){
  /*cout<<"FUNCTION void buble_sort_by_last_name_year_last_names( )\n" ;*/
  int size = this->size( ) ;
  if( size > 1 ){
    node_entry *ptr, *ptrp, *ptrpp, *ptrn, *c ;
    c = firstPtr->next ;
    do{
      ptr = c ; ptrn = ptr->next ; ptrp = ptr->prev ;
      ptrpp = ptrp->prev ; c = c->next ;
      while( ptr->the_entry.compare_by_last_name_year_last_names( ptrp->the_entry ) < 0 ){
	ptr->prev = ptrp->prev;
	if( ptrpp != NULL ){ ptrpp->next = ptr ; }
	ptrp->prev = ptr ;  ptr->next = ptrp ; ptrp->next = ptrn ;
	if( ptrn != NULL ){ ptrn->prev = ptrp ; }
	ptrp = ptr->prev;
	if( ptrp == NULL ){ break ; }
	ptrpp = ptrp->prev ; ptrn = ptr->next ;
      }
    }while( c != NULL ) ;
    /*after ordering, put firstPtr and lastPtr at the begining
      and end or the ordered chain*/
    c = firstPtr->prev ;
    while( c != NULL ){ firstPtr = c ; c = firstPtr->prev ; }
    c = lastPtr->next ;
    while( c != NULL ){ lastPtr = c ;  c = lastPtr->next ;  }
  }
}
/*===============================================================*/
string *list_entry::output_handles( int &n ){
  string *tmp = NULL ;
  n = this->size( ) ;
  if( n>0 ){
    tmp = new string[ n ] ;
    int i = 0 ;
    node_entry *currentPtr = firstPtr ;
    while( currentPtr != NULL ){
      tmp[ i ] = currentPtr->the_entry.get_handle( ) ; i++ ;
      currentPtr = currentPtr->next ;
    }
  }
  return tmp ;
}
/*===============================================================*/
int list_entry::BJ_style_handles( list_alias &BJ_handles ){
  int n = 0 ;
  if( !this->isEmpty( ) ){
    BJ_handles.empty( ) ;
    alias pair ;
    node_entry *currentPtr = firstPtr ;
    while( currentPtr != NULL ){
      pair = currentPtr->the_entry.BJ_style_handle( ) ;
      /*cout<<pair.print( )<<endl;*/
      BJ_handles.insertAtBack( pair ) ; n++ ;
      currentPtr = currentPtr->next ;
    }
  }
  return n ;
}
/*===============================================================*/
void list_entry::print_BJ_style( ofstream &out ){
  if( !this->isEmpty( ) ){
    node_entry *currentPtr = firstPtr ;
    while( currentPtr != NULL ){
      currentPtr->the_entry.print_BJ_style( out ) ;
      currentPtr = currentPtr->next ;
    }
  }
}
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/
/*===============================================================*/
