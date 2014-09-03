#ifndef _SALT_
#define _SALT_

#include<string>
#include "pdbClasses2.h"
#include "miscellanea.h"

static const string all_salts="ACT" ;

typedef enum {
  NON_SALT=-1,
  ACT=0
}salt_XXX ;

const static int N_salt_X=20 ;

static const int n_dmd_salt_t = 6 ;

typedef enum {
  _NON_SALT_ATOM_=0,
  /*acetate anion 2          3           4           5         6   */
  _ACT_CA_=1, _ACT_CB_, _ACT_OG1_, _ACT_BOG1_, _ACT_OG2_, _ACT_BOG2_
}dmd_salt_t ;

static const int Nsalt_atoms = 4 ;
static string salt_atoms[ ]={
"ATOM    107  CA  ACT C  15       0.735   1.331  -0.186  1.00  0.53           C  ",
"ATOM    108  CB  ACT C  15       0.000   0.000   0.000  1.00  0.90           C  ",
"ATOM    109  OG1 ACT C  15       0.661  -0.989   0.270  1.00  1.79           O  ",
"ATOM    110  OG2 ACT C  15      -1.213  -0.005  -0.126  1.00  1.17           O  ",
} ;

salt_XXX string_to_salt_XXX( const string &  ) ;
string dmd_salt_t2string( dmd_salt_t & ) ;
dmd_salt_t string2dmd_salt_t( string & ) ;
bool belong_to_same_XXX( dmd_salt_t &, dmd_salt_t & ) ;
string dmd_salt_t2salt( dmd_salt_t & ) ;
string dmd_salt_t2at( dmd_salt_t & ) ;

/*class SALT derived from class listPDBatom*/

class PDBsalt : public listPDBatom {

  friend ifstream &operator>>(ifstream &, PDBsalt &);

 public:
  PDBsalt( ) ;
  PDBsalt( string &) ;
  PDBsalt( const PDBsalt & ) ;
  PDBsalt &operator=( PDBsalt &  ) ;
  PDBsalt &operator=( const PDBsalt &  ) ;

 protected:
  string name ;
  int saltSeq ;

} ;

#endif /*_SALT_*/
