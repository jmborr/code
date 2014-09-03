#ifndef _PSEUDOTXT_H_
#define _PSEUDOTXT_H_

/*ENERGY FUNCTION*/

/*define the hydrophobic interaction*/
static double HP_E = 0.25;

/*define the salt bridge interaction*/
static double SB_E = 0.5;

/*ANOTHER important interaction is the AROMATIC-type interactions*/
/*First, we just introduce the interaction bwteen PRO and TRP*/
/*LATER, we will add all other interactions*/
static double AROMATIC_E = 0.25;

/*define the mainchain mainchain interaction*/
static double EHB_M_M = 1.0;

/*define the mainchain-sidechain interaction*/
static double EHB_S_M = 1.0;

/*define the Dihedral rotation constraint interaction*/
static double Dihedral_E  = 0.5;

void setHP_E(double e){
  HP_E = e;
}

void setSB_E(double e){
  SB_E = e;
}

void setAROMATIC_E(double e){
  AROMATIC_E = e;
}

void setMMHB_E(double e){
  EHB_M_M = e;
}

void setSMHB_E(double e){
  EHB_S_M = e;
}

void setDI_E(double e){
  Dihedral_E = e;
}
/*END of ENERGY DEFINITION*/

static double BOND_DEV   = 0.03;
static double BOND_DEV_O = 0.03;
static double BOND_DEV_G = 0.050;

void setBOND_DEV(double inpt){
  BOND_DEV = inpt;
}

void setBOND_DEV_G(double inpt){
  BOND_DEV_G = inpt;
}

static double ir_ext = 1.00;

const static double N_HDR  = 1.35;
const static double C_HDR  = 1.50;
const static double O_HDR  = 1.35;
const static double CA_HDR = 1.65;
const static double CB_HDR = 1.65;
const static double CB_IR  = 2.00;

const static double N_M  = 1.0;
const static double C_M  = 1.0;
const static double O_M  = 1.0;
const static double CA_M = 1.0;
const static double CB_M = 1.0;

typedef enum {
  _NON_ATOM_=0,
  /*generic type*/
  _N_=1, _C_, _O_, _CA_, _CB_, _PRO_N_, 
  /*derived type*/
  _N_HB_, _O_HB_,
  /*CB types*/
  _CYS_CB_, _MET_CB_, _PHE_CB_, _ILE_CB_, _LEU_CB_, 
  _VAL_CB_, _TRP_CB_, _TYR_CB_, _ALA_CB_, _GLY_CB_,
  _THR_CB_, _SER_CB_, _GLN_CB_, _ASN_CB_, _GLU_CB_, 
  _ASP_CB_, _HIS_CB_, _ARG_CB_, _LYS_CB_, _PRO_CB_,   
  /*CG types*/
  _CYS_CG_, _MET_CG_, _PHE_CG_, _ILE_CG_, _LEU_CG_, 
  _VAL_CG_, _TRP_CG_, _TYR_CG_, _ALA_CG_, _GLY_CG_,
  _THR_CG_, _SER_CG_, _GLN_CG_, _ASN_CG_, _GLU_CG_, 
  _ASP_CG_, _HIS_CG_, _ARG_CG_, _LYS_CG_, _PRO_CG_, 
  /*CG2 types*/
  _ILE_CG2_, _THR_CG2_, _VAL_CG2_,
  /*derived type*/
  _ASP_HBA_, _ASN_HBA_, _GLN_HBA_, 
  _GLU_HBA_, _SER_HBA_, _THR_HBA_
} dmd_atom_t;


/*define hydrophobicity*/
typedef enum {
  APOLOR=1, AMPHI, POLAR, IRRELEVANT
} hp_t;

const static hp_t dmd_atom_hp[]={
  /*CB*/
  IRRELEVANT, IRRELEVANT, IRRELEVANT, IRRELEVANT, IRRELEVANT,
  IRRELEVANT, IRRELEVANT, IRRELEVANT, APOLOR,     IRRELEVANT,
  IRRELEVANT, IRRELEVANT, IRRELEVANT, IRRELEVANT, IRRELEVANT,
  IRRELEVANT, IRRELEVANT, IRRELEVANT, IRRELEVANT, IRRELEVANT,
  /*CG1*/
  APOLOR, APOLOR, APOLOR, APOLOR, APOLOR,
  APOLOR, APOLOR, APOLOR, IRRELEVANT, IRRELEVANT,
  POLAR, POLAR, POLAR, POLAR, POLAR, 
  POLAR, POLAR, POLAR, POLAR, POLAR,
  /*CG2*/
  APOLOR, POLAR, APOLOR
};

const static int is_aromatic[]={
  /*CB*/
  0,     0,     0,     0,     0,
  0,     0,     0,     0,     0,
  0,     0,     0,     0,     0,
  0,     0,     0,     0,     0,
  /*CG1*/
  0,     0,     1,     0,     0,
  0,     1,     1,     0,     0,
  0,     0,     0,     0,     0,
  0,     0,     0,     0,     1,
  /*CG2*/
  0,     0,     0
};

const static double dmd_atom_hp_scale[]={
  /*CB*/
  0.0,     0.0,     0.0,     0.0,     0.0,
  0.0,     0.0,     0.0,    -0.40,    0.0,
  0.0,     0.0,     0.0,     0.0,     0.0,
  0.0,     0.0,     0.0,     0.0,     0.0,
  /*CG1*/
 -0.56,   -0.42,   -0.62,   -1.00,   -0.84,
 -0.93,   -0.20,   -0.20,    0.00,    0.00,
  0.16,    0.18,    0.78,    0.78,    0.78,
  0.78,    0.02,    1.00,    0.87,    0.36,
  /*CG2*/
 -1.00,    0.16,   -0.93
};
/*end of hydrophobicity definition*/

/*end of hydrophobicity defintion*/


/***********************************************************************
  Define the Dihedral (actually correspoding to Chi1 rotamer
  using the distance between Ci and CGi to implement the constraints
  chi1: 0, 240 -- [d1, d2]
  chi1: 120    -- [d3, ..)
***********************************************************************/
const static double Dihedral_MIN = 2.0;
const static double Dihedral_MAX = 6.0;
const static double Dihedral_C_CG[20][3]={
  //CYS
  3.00, 3.32, 4.06,
  //MET: no-definition
  INF, INF, INF,
  //PHE
  3.73, 4.22, 5.10,
  //ILE: 
  2.80, 3.05, 3.79,
  //LEU
  3.22, 3.55, 4.25,
  //3.42, 3.78, 4.23,
  //VAL
  2.80, 3.05, 3.79,
  //TRP: no-def
  4.12, 4.58, 5.53,
  //TYR
  4.00, 4.54, 5.47,
  //ALA: no-def
  INF, INF, INF,
  //GLY: no-def
  INF, INF, INF,
  //THR
  2.82, 3.10, 3.77,
  //SER
  2.68, 3.06, 3.68,
  //GLN: no-def
  INF, INF, INF,
  //ASN
  3.12, 3.48, 4.16,
  //GLU: no-def
  INF, INF, INF,
  //ASP
  3.12, 3.48, 4.16,
  //HIS
  3.57, 4.05, 4.87,
  //ARG: no-def
  INF, INF, INF,
  //LYS: no-def
  INF, INF, INF,
  //PRO: no-def
  INF, INF, INF
};
/*end of Dihedral definition*/

/*define the geometry of hydrogen bonding*/
/*First, only the backbone HB*/
const static double HB_N_O[2] ={2.80, 3.12};
const static double HB_N_C[4] ={3.70, 3.80, 3.91, 4.23};
const static double HB_O_C[3] ={3.48, 3.60, 4.00};
const static double HB_O_CA[3]={3.47, 3.60, 4.04};
const static double HB_MAX = 6.0;
const static double max_height = 1.0E8;
const static double eps = 1.0E-6;
/*Second, it is believed that the Side-Main chain HB are important.
  We will introduce it later*/
const static double HB_ASP_N[2] = {3.52, 4.04};
const static double HB_ASP_C[2] = {4.42, 4.94};
const static double HB_ASP_CA[2]= {4.08, 4.76};

const static double HB_ASN_N[2] = {3.52, 4.04};
const static double HB_ASN_C[2] = {4.42, 4.94};
const static double HB_ASN_CA[2]= {4.08, 4.76};

const static double HB_GLN_N[2] = {3.60, 4.04};
const static double HB_GLN_C[2] = {4.38, 4.86};
const static double HB_GLN_CA[2]= {4.30, 4.90};

const static double HB_GLU_N[2] = {3.60, 4.04};
const static double HB_GLU_C[2] = {4.38, 4.86};
const static double HB_GLU_CA[2]= {4.30, 4.90};

const static double HB_SER_N[2] = {2.87, 3.27};
const static double HB_SER_C[2] = {3.77, 4.23};
const static double HB_SER_CA[2]= {3.64, 4.08};

const static double HB_THR_N[2] = {2.87, 3.27};
const static double HB_THR_C[2] = {3.77, 4.23};
const static double HB_THR_CA[2]= {3.64, 4.08};
/*endof geometry of hydrogen bonding*/

/*define interaction min-dist-matrix: 5x5*/
const static double min_NN   = 2.75;
const static double min_NC   = 3.15;
const static double min_NO   = 2.80;
const static double min_NCa  = 3.36;//not used
const static double min_NCb  = 2.85;
const static double min_CC   = 3.52;//not used
const static double min_CO   = 2.85;//not used
const static double min_CCa  = 3.45;//not used
const static double min_CCb  = 3.05;
const static double min_OO   = 2.70;//not used
const static double min_OCa  = 2.70;
const static double min_OCb  = 2.70;
const static double min_CaCa = 3.41;//not used
const static double min_CaCb = 3.24;//not used
const static double min_CbCb = 3.07;//not used


string txtKeyWords[]={
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

#endif
