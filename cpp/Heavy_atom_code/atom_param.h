#ifndef _ATP_
#define _ATP_
#define _AP_INF_ 999.999

static const string s_inf("999.999") ;

static const string datKeyWords[]={
  "C.TYPES OF ATOMS"                        ,
  "F.LINKED PAIRS"                          ,
  "LIST ROTAMER BARRIERS"                   ,
  "HIDROGEN BOND LINKS"                     ,
  "IJ.LIST OF PERMANENT BONDS"              ,
  "VAN DER WAALS"                           ,
  "ELECTROSTATICS"                          ,
  "HYDROPHOBICITY"                          ,
  "HYDROGEN BOND STRENGTH"                  ,
  "NON BONDED SIGNAL HYDROGEN BOND"         ,
  "G.REACTIONS"                             ,
  "L.LIST OF HYDROGEN BONDING ASSOCIATIONS" ,
};

typedef enum {  
  NON_KEY_DAT=-1        ,
  TYPE_ATOMS_DAT        , 
  LINK_PAIRS_DAT        ,
  ROTAMER_BARRIERS_DAT  ,
  LINK_HB_DAT           ,
  LIST_PERM_BONDS_DAT   ,
  VAN_DER_WAALS_DAT     ,
  ELECTROSTATICS_DAT    ,
  HYDROPHOBICITY_DAT    ,
  HB_STRENGTH_DAT       ,
  SIGNAL_HB_DAT         ,
  REACT_DAT             ,
  HB_LIST_DAT
} dat_key;


bool is_comment( string & ) ;
bool has_comment( string & ) ;
bool goto_dat_key( ifstream &, const dat_key & ) ;
bool remove_trailing_comment( string & ) ;
bool is_end( string & ) ;
bool is_end( string &, const dat_key &  ) ;
void print_atom_type_line( char *, ifstream &, string & ) ;
double out_hard_core_radius( ifstream &, string & ) ;

double out_charge( ifstream &, string & ) ;
double charge_of( string &stype, ifstream &DAT ) ;
bool is_charged( string &, ifstream & ) ;
bool get_electrostatics( string &, string &, string &, ifstream & ) ;
bool get_default_electrostatics( string &, string &, string &, ifstream & ) ;
bool get_default_electrostatics( const double&, const double&, const double&,
				 const double&, string &, ifstream &       ) ;
void rescale_default_electrostatics_by_factor( string &, const double & ) ;

double out_hydrophobicity( ifstream &, string & ) ;
bool get_default_hydrophobicity( string &, string &, string &, ifstream & ) ;
bool get_default_hydrophobicity( const double&, const double&, const double&,
				 const double&, string &, ifstream &       ) ;
void rescale_default_hydrophobicity_by_factor( string &, const double & ) ;

bool is_HB( string & ) ;
bool is_HB_link( string &, ifstream & ) ;
double get_HB_strength( ifstream & ) ;
void get_signal_Hydrogen_Bond( string &, string &, string &, ifstream & ) ;

bool out_linked_pair( string &, string &, string &, char * ) ;
bool out_linked_pairII( ifstream &, string &, string &, char * ) ;
bool is_there_linked_pair( ifstream &, string &, string &, string & ) ;
bool out_prerelaxed_linked_pair( string &, string &, string &, char * ) ;
string *out_all_same_permanent_bonds(string &,string &,ifstream &,int &);
string *out_all_next_permanent_bonds(string&,string&,string&,ifstream&,int&);
void  filter_friends_same( string *, string &, int & ) ;

double get_inner_VW_factor( ifstream & ) ;
void get_default_Van_der_Waals( string &, string &, string &, ifstream & ) ;

bool do_they_react( ifstream &, string &, string &, string & ) ;
bool get_reactants_products( string &,string &,string &, string &, string & );
bool get_associates( string &, string *&, int *&, string *&, int & ) ;

bool have_rotamers( ifstream &, string &, string & ) ;
double get_default_rotamer_barrier( ifstream & ) ;
void default_rotamer_barrier( ifstream &, string &, char * ) ;

void rescale_barriers_of_link_by_factor( string &, const double & ) ;
#endif/*ATP*/
