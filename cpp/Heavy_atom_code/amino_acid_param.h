#ifndef _AAP_
#define _AAP_

const static int N_aa_XXX=20;

typedef enum {
  NON=-1,
  GLY=0, ALA , VAL, LEU, ILE, SER, THR, CYS, MET, PRO,
  ASP, ASN, GLU, GLN, LYS, ARG, HIS, PHE, TYR, TRP
}aa_XXX;

const static int N_aa_X=20;

typedef enum {
  Z=-1,
  G=0, A, V, L, I, S, T, C, M, P, D, N, E, Q, K, R, H, F, Y, W
}aa_X;


static const string txtKeyWords[]={
  "A.SYSTEM SIZE",
  "B.NUMBER OF ATOMS",
  "C.TYPES OF ATOMS",
  "D.NON-ELASTIC COLLISIONS",
  "E.ELASTIC COLLISIONS",
  "F.LINKED PAIRS",
  "G.REACTIONS",
  "H.LIST OF ATOMS",
  "I.LIST OF BONDS",
  "IJ.LIST OF PERMANENT BONDS",
  "J.BOND TABLE LENGTH",
  "K.LIST OF PARAMETERS",
  "COLLISION TABLE LENGTH",
  "L.LIST OF HYDROGEN BONDING ASSOCIATIONS",
  "M.LIST OF GID"
};

typedef enum {  
  NON_KEY=-1,
  SYS_SIZE, 
  NUM_ATOMS, 
  TYPE_ATOMS, 
  NONEL_COL, 
  EL_COL, 
  LINK_PAIRS, 
  REACT, 
  LIST_ATOMS, 
  LIST_BONDS, 
  LIST_PERM_BONDS, 
  NUM_BONDS, 
  LIST_PARAM, 
  COL_TABLE,
  HB_LIST,
  GID_LIST
} dmd_txt_key;


static const int n_dmd_atom_t = 110 ;

typedef enum {
  _NON_ATOM_=0,

  /*backbone (BKB)          3         4        5        6         7*/
  _BKB_N_=1, _BKB_N_HB_, _PRO_N_, _BKB_CA_, _BKB_C_, _BKB_O_, _BKB_O_HB_,

  /*ALA*/
  _ALA_CB_,

  /*ARG         10       11        12        13   */
  _ARG_CB_, _ARG_CG_, _ARG_CD_, _ARG_NE_, _ARG_CZ_,
  /*  14        15            16          17      */ 
  _ARG_NH1_, _ARG_NH1_HB_, _ARG_NH2_, _ARG_NH2_HB_,

  /*ASN        19         20          21          22           23     */
  _ASN_CB_, _ASN_CG_, _ASN_OD1_, _ASN_OD1_HB_, _ASN_ND2_, _ASN_ND2_HB_,

  /*ASP        25         26          27          28           29     */
  _ASP_CB_, _ASP_CG_, _ASP_OD1_, _ASP_OD1_HB_, _ASP_OD2_, _ASP_OD2_HB_,

  /*CYS        31    */
  _CYS_CB_, _CYS_SG_,

  /*GLN        33        34         35    */
  _GLN_CB_, _GLN_CG_, _GLN_CD_, _GLN_OE1_,
  /*   36           37          38        */
  _GLN_OE1_HB_, _GLN_NE2_, _GLN_NE2_HB_,

  /*GLU        40         41       42     */
  _GLU_CB_, _GLU_CG_, _GLU_CD_, _GLU_OE1_,
  /*    43          44          45        */
  _GLU_OE1_HB_, _GLU_OE2_, _GLU_OE2_HB_,

  /*GLY*/

  /*HIS        47         48          49           50    */
  _HIS_CB_, _HIS_CG_, _HIS_ND1_, _HIS_ND1_HB_, _HIS_CD2_,
  /*  51         52          53      */
  _HIS_CE1_, _HIS_NE2_, _HIS_NE2_HB_,

  /*ILE         55         56         57    */
  _ILE_CB_, _ILE_CG1_, _ILE_CG2_, _ILE_CD1_,

  /*LEU         59        60         61    */
  _LEU_CB_, _LEU_CG_, _LEU_CD1_, _LEU_CD2_,

  /*LYS        63        64        65        66          67     */
  _LYS_CB_, _LYS_CG_, _LYS_CD_, _LYS_CE_, _LYS_NZ_, _LYS_NZ_HB_,

  /*MET        69        70        71    */
  _MET_CB_, _MET_CG_, _MET_SD_, _MET_CE_,

  /*PHE        73        74         75         76         77         78    */
  _PHE_CB_, _PHE_CG_, _PHE_CD1_, _PHE_CD2_, _PHE_CE1_, _PHE_CE2_, _PHE_CZ_,

  /*PRO        80        81     */
  _PRO_CB_, _PRO_CG_, _PRO_CD_,

  /*SER        83         84      */
  _SER_CB_, _SER_OG_, _SER_OG_HB_,

  /*THR        86           87           88   */
  _TRH_CB_, _TRH_OG1_, _TRH_OG1_HB_, _TRH_CG2_, 

  /*TRP        90         91         92        93         94     */
  _TRP_CB_, _TRP_CG_, _TRP_CD1_, _TRP_CD2_, _TRP_NE1_, _TRP_CE2_,
  /*  95         96         97         98    */
  _TRP_CE3_, _TRP_CZ2_, _TRP_CZ3_, _TRP_CH2_,

  /*TYR        100       101         102       103     */
  _TYR_CB_, _TYR_CG_, _TYR_CD1_, _TYR_CD2_, _TYR_CE1_,
  /*  104       105       106        107     */
  _TYR_CE2_, _TYR_CZ_, _TYR_OH_, _TYR_OH_HB_,

  /*VAL        109        110   */
  _VAL_CB_, _VAL_CG1_, _VAL_CG2_
} dmd_atom_t;

aa_XXX string_to_aa_XXX( const string &  ) ;
string dmd_atom_t2string( dmd_atom_t & ) ;
dmd_atom_t string2dmd_atom_t( string & ) ;
bool is_HB( dmd_atom_t & ) ;
bool belong_to_same_XXX( dmd_atom_t &, dmd_atom_t & ) ;
string dmd_atom_t2amino( dmd_atom_t & ) ;
string dmd_atom_t2at( dmd_atom_t & ) ;
void goto_txt_key( ifstream &, const dmd_txt_key & ) ;
#endif /*_AAP_*/
