#ifndef _BIBL_
#define _BIBL_
#include<iostream.h>
using namespace std;
#include<string>

bool is( string &, string & ) ;
void remove_spaces( string & ) ;
void remove_last_colon( string & ) ;
int remove_lasts( const char *, string & ) ;
int remove_firsts( const char *, string & ) ;
int replace( const char *, const char *, string & ) ;
bool replace_last( const char *, const char *, string & ) ;
int remove( const char *, string & ) ;
string char_to_string ( const char * ) ;
int how_many( const char * , string & ) ;
int how_many_nonDegenerated( const char * , string & ) ;
string * split( const char * , string &, int & ) ;
string * split_nonDegenerated( const char * , string &, int & ) ;
bool last_name_first( string & ) ;
bool leave_only_last_name( string & ) ;
bool find_next( const char *, ifstream & ) ;
int find_nex2( const char *, FILE * ) ;
int my_getc( const char *, char *, FILE * ) ;
bool print_till_next( const char *, FILE *, ofstream & ) ;
string *get_handles_from_cite_command( string &, int & ) ;
string *create_list_of_handles( int &, ifstream & ) ;
string *sort_by_last_name_year_last_names( string * , const int &, ifstream & ) ;

/*===============  C L A S S   A L I A S    =====================*/
void parse_alias( const string &, string &, string & ) ;
bool is_alias( const string & ) ;

class alias{

  friend class node_alias ;
  friend class list_alias ;
  friend class field ;

 public:
  alias( ) ;
  alias( const string & ) ;
  void fill_members( string & ) ;
  string print( ) ;
  void assign( const string &, const string & );

 private:
  string a ;
  string full  ;
};
/*===========================================================*/


/*===============  C L A S S   N O D E _ A L I A S =================*/
class node_alias{

  friend ostream &operator<<(ostream &, list_alias & ) ;

  friend class list_alias ;
  friend class field ;

 public:
  node_alias( const alias& ) ;

 private:
  alias the_alias ;
  node_alias *next ;
};
/*===========================================================*/



/*===============  C L A S S   L I S T _ A L I A S =================*/
class list_alias{

  friend ifstream &operator>>(ifstream &, list_alias & ) ;
  friend ostream &operator<<(ostream &, list_alias & ) ;

  friend class alias ;
  friend class node_alias ;
  friend class field ;

 public:
  list_alias( ) ;
  ~list_alias( ) ;
  bool isEmpty( ) ;
  void rmFromFront( ) ;
  void empty( ) ;
  void insertAtFront( const alias & );
  void insertAtBack( const alias & );
  void insertAtBack( const string & );
  int size( ) ;
  int create_list_of_handles( ifstream & ) ;
  bool is_there_alias( const string & ) ;
  string *output_aliases( int & ) ;
  string get_the_full( const string & ) ;

 private:
  node_alias *firstPtr, *lastPtr ;
  node_alias *getNewNode( const alias & ) ;
};
/*===========================================================*/



/*===============  C L A S S   F I E L D    =================*/

typedef enum { 
  _NON_FIELD_=0,
  _AUTH_=1 , _TITLE_, _JOURN_, _VOL_, _PAGES_, _YEAR_,
  _PUBLISHER_, _ADDR_,  _CROSSREF_
}field_kw ;

const static int nun_field_keywords= 10;
const static string field_keywords[] = {
  "non_field",
  "author",
  "title",
  "journal",
  "volume",
  "pages",
  "year",
  "publisher",
  "address",
  "crossref"
};

field_kw determine_field_kw( const string & ) ;
bool is_field( const string & ) ;

class field{

  friend class node_field ;
  friend class entry ;

  friend ostream &operator<<( ostream &, const field & ) ;
  friend ifstream &operator>>(ifstream &, field &);

 public:
  field( ) ;
  field( const string & ) ;
  bool assign( const string & ) ;
  string get_value( ) ;
  string get_keyword( ) ;
  field_kw get_keyword2( ) ;
  int get_keyword3( ) ;
  string print( ) ;
  bool last_names_first( ) ;    /*only for author field*/
  string *last_names( int & ) ;/*only for author field*/
  bool format_author_for_BJ( ) ;/*only for author field*/
  bool unalias( const list_alias & ) ;


 private:
  field_kw kw ;
  string value ;
};

/*===============  C L A S S   N O D E _ F I E L D    =================*/
class node_field{

  friend ostream &operator << ( ostream &, entry & ) ;

  friend class entry ;

 public:
  node_field( const field & ) ;

 private:
  field the_field ;
  node_field *next ;
};
/*===============================================================*/


/*===============  C L A S S   E N T R Y    =================*/

typedef enum { 
  _NON_ENTRY_=0,
  _BOOK_=1 , _ARTICLE_, _INPROCEEDINGS_
}entry_kw ;

const static int nun_entry_keywords= 4;
const static string entry_keywords[] = {
  "non_entry",
  "book",
  "article",
  "inproceedings"
};

entry_kw determine_entry( const string & );
string determine_handle( const string & );
bool is_entry( const string & ) ;
bool is_end_of_entry( const string & ) ;
bool has_handle( const string & ) ;

class entry{

  friend ifstream &operator >> ( ifstream &, entry & ) ;
  friend ostream &operator << ( ostream &, entry & ) ;

  friend class field;
  friend class node_field;

 public:
  entry( ) ;
  ~entry( ) ;
  void empty( ) ;
  bool isEmpty( ) ;
  void rmFromFront( ) ;
  void insertAtBack( const field & );
  void insertAtBack( const string  & );
  bool fill( ifstream &, const string & ) ;
  bool is_there_field( const field_kw & ) ;
  bool get_field( const field_kw &, field & ) ;
  string get_field( const field_kw & ) ;
  string get_handle( ) ;
  bool last_names_first( ) ;
  int unalias( const list_alias & ) ;
  bool format_for_BJ( ) ;
  entry &operator=( entry & ) ;
  string *last_names( int & ) ;
  string *last_name_year_last_names( int & ) ;
  int compare_by_last_name_year_last_names( entry & ) ;
  alias BJ_style_handle( ) ;
  void print_BJ_style( ofstream & ) ;

 private:
  entry_kw kw ;
  string handle ;
  node_field *firstPtr, *lastPtr ;
  node_field *getNewNode( const field & ) ;
};
/*===============================================================*/


/*===============  C L A S S   N O D E _ E N T R Y    =================*/

class node_entry{

  friend class list_entry ;
  friend ostream &operator << ( ostream &, list_entry & ) ;

 public:
  node_entry( entry & ) ;

 private:
  entry the_entry ;
  node_entry *next ;
  node_entry *prev ;
};
/*===============================================================*/


/*===============  C L A S S   L I S T _ E N T R Y    =================*/

class list_entry{

  friend ostream &operator<<( ostream &, list_entry & ) ;

 public:
  list_entry( ) ;
  ~list_entry( ) ;
  void empty( ) ;
  bool isEmpty( ) ;
  void rmFromFront( ) ;
  void insertAtBack( entry & );
  void fill_with_list_handles( string *, const int &, ifstream & ) ;
  int size( ) ;
  int unalias( const list_alias & ) ;
  void buble_sort_by_last_name_year_last_names( ) ;
  string *output_handles( int & ) ;
  int BJ_style_handles( list_alias & ) ;
  void print_BJ_style( ofstream & ) ;

 private:
  node_entry *firstPtr ;
  node_entry *lastPtr ;
  node_entry *getNewNode( entry & ) ;
};

#endif/*_BIBL_*/
