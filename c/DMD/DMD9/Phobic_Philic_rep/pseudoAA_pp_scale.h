/**************************************************************************************
 ****  The next step after 4-bead model is to include the side-chains which will    ***
 ****  introduce the side-chain entropy to the model. It is believed that the	    ***
 ****  side-chain entropy have a strong contribution to the helix propensity. It is ***
 ****  also believed that the etropic contribution to protein folding is equally    ***
 ****  improtant to enthalpic contribution. Therefore, here we include the          ***
 ****  side-chains our previous model. For simplicity we only add one CG atom to    ***
 ****  each amino acids, and for those b-branched amino acids we add two CG atoms:  ***
 ****  VAL, THR, ILE.								    ***
 ****  										    ***
 ****   ( CG2 )------------ CG1							    ***
 ****      \  *            * /							    ***
 ****       \   *       *   /							    ***
 ****        \	  *  *     /							    ***
 ****         \     CB    /							    ***
 ****          \    |    /							    ***
 ****           \   |   /							    ***
 ****            \  |  /							    ***
 ****             \ | /								    ***
 ****               CA								    ***
 ****            /     \							    ***
 ****  	       /         \							    ***
 ****        /             \							    ***
 ****      /                 \							    ***
 ****     N                    C'						    ***
 ****                          |						    ***
 ****                          |						    ***
 ****                         (O)						    ***
 ****  										    ***
 ****  The rotational freedom can characterized by chi, theta(CA-CB-CG1).	    ***
 **************************************************************************************
 ****   Feng Ding
 ****   fding@bu.edu
 ***********************/
#ifndef _PSEUDOAA_H_
#define _PSEUDOAA_H_

#include "PDBLib.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
using namespace std;

const static int N_aa_t=20;
typedef enum {
  NON=-1, 
  CYS=0, MET, PHE, ILE, LEU, VAL, TRP, TYR, ALA, GLY, 
  THR, SER, GLN, ASN, GLU, ASP, HIS, ARG, LYS, PRO
}aa_t;

/*define the GEOMETRY for each amino acid*/
const static double N_CA_CB    = 110.4;
const static double CB_N_CA_C  = 122.0;
const static double N_CA_CBXY  = 180.0-acos(cos(N_CA_CB*rpi)/cos(CB_N_CA_C*rpi))/rpi;
const static double PRO_N_CA_CB  = 103.0;
const static double PRO_N_CA_CBXY= 180.0-acos(cos(PRO_N_CA_CB*rpi)/cos(CB_N_CA_C*rpi))/rpi;

const static double N_CA_C   = 111.2;
const static double CA_C_O   = 120.8;
const static double C_N_CA   = 121.7;
const static double CA_C_N   = 116.2;

const static double C_CA_CBXY = 360.0-N_CA_CBXY-N_CA_C;
const static double C_CA_CB   = 180.0-acos(cos(CB_N_CA_C*rpi)*cos(C_CA_CBXY*rpi))/rpi;
const static double PRO_C_CA_CBXY = 360.0-PRO_N_CA_CBXY-N_CA_C;
const static double PRO_C_CA_CB   = 180.0-acos(cos(CB_N_CA_C*rpi)*cos(PRO_C_CA_CBXY*rpi))/rpi;

const static double N_CA     = 1.458;
const static double CA_C     = 1.525;
const static double C_O      = 1.231;
const static double N_H      = 1.020;
const static double C_N      = 1.329;
const static double CA_CB    = 1.521;

/*next is about the gamma atom GEOMETRY*/
//Gamma1 atom
const static int N_G1_data_t = 7;
typedef enum{G1_CB=0, G1_CB_D, G1_CA, G1_CA_D, G1_HDR, G1_IR, G1_M} G1_data_t;
/*----------------------------------------------------------------------------------------------
  0: CG-CB; 1: DELTA(CG-CB); 2: CG-CA; 3: DELTA(CG-CA); 4: Hardcore Radis; 5: Interaction Range
  ----------------------------------------------------------------------------------------------
  A: If the Delta is zero, means a rigid bond and will be assigned a common deviation coef.
  otherwise, the bonds is flexible and it is the deivation of the "flexible-bonds"
  
  B: If the Interaction Range is set zero, the GAMMA atom is not interting by
  hydrophobic interactoin (HP) and salt bridge (SB) interactions. But it may
  still interact with the backbone to form the m-s hb.
  ---------------------------------------------------------------------------------------------*/
const static double G1_CONST[N_aa_t][N_G1_data_t]={
  //CYS: disulf bond length 2\AA(S-S): H (hydrophobic)(S)
  1.83, 0.00, 2.80, 0.00, 1.80, 2.20, 1.0,//2.20
  //MET: H(S)
  2.76, 0.00, 3.71, 0.30, 1.65, 2.90, 3.0,//2.90
  //PHE: H
  2.91, 0.00, 3.79, 0.00, 1.65, 3.15, 6.0,//3.15
  //ILE: (CG1) H --- beta-branched
  1.52, 0.00, 2.52, 0.00, 1.65, 2.20, 1.0,//2.20
  //LEU: H
  1.94, 0.00, 3.04, 0.00, 1.65, 3.00, 3.0,//3.00
  //VAL: --- beta-branched
  1.52, 0.00, 2.50, 0.00, 1.65, 2.20, 1.0,//2.20
  //TRP: H(treated as a rigid object) AR+AMPHI
  3.33, 0.00, 4.25, 0.00, 1.65, 3.50, 9.0,//3.50

  //TYR: Amphi + AR
  3.32, 0.00, 4.16, 0.00, 1.65, 3.25, 7.0,//3.25
  //ALA: no GAMMA atom
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0,
  //GLY: not either
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.0,
  //THR: AMPHI --- beta-branched
  1.52, 0.00, 2.49, 0.00, 1.65, 2.20, 1.0,//2.20
  //SER: P(O)
  1.45, 0.00, 2.43, 0.00, 1.35, 2.00, 1.0,//G1_IR modified by Jose (prev=0.00)
  //GLN: HBA + P
  2.47, 0.00, 3.40, 0.30, 1.65, 2.00, 4.0,//G1_IR modified by Jose (prev=0.00)
  //ASN: HBA (hydrogen bond accepter) & P (polar)
  1.94, 0.00, 2.88, 0.00, 1.65, 2.00, 3.0,//G1_IR modified by Jose (prev=0.00)
  //GLU: HBA + SBN
  2.47, 0.00, 3.40, 0.30, 1.65, 2.00, 4.0,
  //ASP: HBA & SBN (salt-bridge netative-) 
  1.94, 0.00, 2.88, 0.00, 1.65, 2.00, 3.0,
  //HIS: rigid ring structure + P
  2.65, 0.00, 3.55, 0.00, 1.65, 2.00, 5.0,//G1_IR modified by Jose (prev=0.00)
  //ARG: SBP (hydrophobic)
  4.45, 0.45, 5.55, 0.55, 1.65, 3.00, 6.0,
  //LYS: SBP (salt-bridge positive +)
  3.40, 0.35, 4.60, 0.40, 1.65, 3.00, 4.0,
  //PRO: "P"
  1.83, 0.00, 2.28, 0.00, 1.85, 2.50, 2.0//2.50
};
//Gamma2 atoms
const static int N_G2_data_t=7;
typedef enum{G2_CB=0, G2_CA, G2_HDR, G2_IR, G1_G2, G2_CA_CB_G1, G2_M} G2_data_t;

const static double ILE_G2[N_G2_data_t]={1.94, 2.87, 1.65, 2.75, 2.85, 124.7, 2};
/* G2_IR of THR_G2 modified by Jose (prev.=0.00)*/
const static double THR_G2[N_G2_data_t]={1.43, 2.41, 1.35, 2.00, 2.41, 121.6, 1};//O
const static double VAL_G2[N_G2_data_t]={1.52, 2.50, 1.65, 2.00, 2.49, 125.1, 1};

/*for PROLINE, it is different from others*/
const static double PRO_N_CG   = 1.80;
const static double PRO_N_CG_D = 0.09;
const static double PRO_C_CG   = 3.01;
const static double PRO_C_CG_D = 0.50;
const static double PRO_C_C1_MIN = 2.90;
const static double PRO_C_C1_MAX = 3.50;
const static double PRO_C_CB   = 3.63;
const static double PRO_C_CB_D = 0.05;

/*end of GEOMETRY definition*/

/*Determine the intrinsic axis of the amino acids
  FOR INITIALIZING AMIN ACID STRUCTURE
  
  Two Streched amino acids(psi = -180, and phi = 180)
                     (ix)
       Ca--(iy)   N   |    C
     / a|  \    /  \ a|  /
    N   |    C       Ca--(iy)
       (ix)   
  angle : a is such that the two intrinsic x axis are parallel
  ------------------------------*/
double C_CA     = sqrt(C_N*C_N + N_CA*N_CA - 2.0*C_N*N_CA*cos(C_N_CA*rpi));
double N_C0_CA  = asin(N_CA*sin(C_N_CA*rpi)/C_CA)/rpi;
double C_CA_N0  = 180.0 - N_C0_CA - C_N_CA;
const static double CA_C_CA  = CA_C_N + N_C0_CA;

const static double CA_CA = sqrt(CA_C*CA_C + C_CA*C_CA - 2.0*CA_C*C_CA*cos(CA_C_CA*rpi));

double CA_CA_C0 = asin(CA_C*sin(CA_C_CA*rpi)/CA_CA)/rpi;
const static double C_CA0_CA = 180.0 - CA_C_CA - CA_CA_C0;
double CA_CA_N0 = C_CA_N0 - CA_CA_C0;

const static double init_N_CA_x   = (N_CA_C + C_CA0_CA - CA_CA_N0)/2.0;
const static double init_C_CA_x   = N_CA_C - init_N_CA_x;

const char* cstr_aa3[]={
  "CYS", "MET", "PHE", "ILE", "LEU", "VAL", "TRP", "TYR", "ALA", "GLY", 
  "THR", "SER", "GLN", "ASN", "GLU", "ASP", "HIS", "ARG", "LYS", "PRO"
};

const char* cstr_aa1[]={
  "C",   "M",   "F",   "I",   "L",   "V",   "W",   "Y",   "A",   "G",
  "T",   "S",   "Q",   "N",   "E",   "D",   "H",   "R",   "K",   "P"
};

aa_t aa_Name2Type(const char* aa_name){
  if(strlen(aa_name)==3){
    if(!strcmp(aa_name, cstr_aa3[CYS])) return CYS;
    if(!strcmp(aa_name, cstr_aa3[MET])) return MET;
    if(!strcmp(aa_name, cstr_aa3[PHE])) return PHE;
    if(!strcmp(aa_name, cstr_aa3[ILE])) return ILE;
    if(!strcmp(aa_name, cstr_aa3[LEU])) return LEU;
    if(!strcmp(aa_name, cstr_aa3[VAL])) return VAL;
    if(!strcmp(aa_name, cstr_aa3[TRP])) return TRP;
    if(!strcmp(aa_name, cstr_aa3[TYR])) return TYR;
    if(!strcmp(aa_name, cstr_aa3[ALA])) return ALA;
    if(!strcmp(aa_name, cstr_aa3[GLY])) return GLY;
    if(!strcmp(aa_name, cstr_aa3[THR])) return THR;
    if(!strcmp(aa_name, cstr_aa3[SER])) return SER;
    if(!strcmp(aa_name, cstr_aa3[GLN])) return GLN;
    if(!strcmp(aa_name, cstr_aa3[ASN])) return ASN;
    if(!strcmp(aa_name, cstr_aa3[GLU])) return GLU;
    if(!strcmp(aa_name, cstr_aa3[ASP])) return ASP;
    if(!strcmp(aa_name, cstr_aa3[HIS])) return HIS;
    if(!strcmp(aa_name, cstr_aa3[ARG])) return ARG;
    if(!strcmp(aa_name, cstr_aa3[LYS])) return LYS;
    if(!strcmp(aa_name, cstr_aa3[PRO])) return PRO;
  }
  else if(strlen(aa_name)==1){
    if(!strcmp(aa_name, cstr_aa1[CYS])) return CYS;
    if(!strcmp(aa_name, cstr_aa1[MET])) return MET;
    if(!strcmp(aa_name, cstr_aa1[PHE])) return PHE;
    if(!strcmp(aa_name, cstr_aa1[ILE])) return ILE;
    if(!strcmp(aa_name, cstr_aa1[LEU])) return LEU;
    if(!strcmp(aa_name, cstr_aa1[VAL])) return VAL;
    if(!strcmp(aa_name, cstr_aa1[TRP])) return TRP;
    if(!strcmp(aa_name, cstr_aa1[TYR])) return TYR;
    if(!strcmp(aa_name, cstr_aa1[ALA])) return ALA;
    if(!strcmp(aa_name, cstr_aa1[GLY])) return GLY;
    if(!strcmp(aa_name, cstr_aa1[THR])) return THR;
    if(!strcmp(aa_name, cstr_aa1[SER])) return SER;
    if(!strcmp(aa_name, cstr_aa1[GLN])) return GLN;
    if(!strcmp(aa_name, cstr_aa1[ASN])) return ASN;
    if(!strcmp(aa_name, cstr_aa1[GLU])) return GLU;
    if(!strcmp(aa_name, cstr_aa1[ASP])) return ASP;
    if(!strcmp(aa_name, cstr_aa1[HIS])) return HIS;
    if(!strcmp(aa_name, cstr_aa1[ARG])) return ARG;
    if(!strcmp(aa_name, cstr_aa1[LYS])) return LYS;
    if(!strcmp(aa_name, cstr_aa1[PRO])) return PRO;
  }
  else return NON;
}

/*return the Euler Matrix by(theta, phi, psi), fixing the axis*/
void getCompEularMatrix(double** m, 
			double theta, double phi, double psi){
  double cth=cos(theta);
  double sth=sin(theta);
  double cphi = cos(phi);
  double sphi = sin(phi);
  double cpsi = cos(psi);
  double spsi = sin(psi);
  m[0][0] = cphi*cpsi - sphi*spsi*cth;
  m[0][1] = -cphi*spsi-sphi*cpsi*cth;
  m[0][2] = sphi*sth;
  m[1][0] = sphi*cpsi+cphi*spsi*cth;
  m[1][1] = -sphi*spsi+cphi*cpsi*cth;
  m[1][2] = -cphi*sth;
  m[2][0] = spsi*sth;
  m[2][1] = cpsi*sth;
  m[2][2] = cth;
}

/*return the Euler Matrix by(theta, phi, psi), fixing the space*/
void getAxisEularMatrix(double** m, 
			double theta, double phi, double psi){
  getCompEularMatrix(m, -theta, -psi, -phi);
}

double** newRMatrix(){
  double** rotate = new double*[3];
  rotate[0] = new double[9];
  rotate[1] = rotate[0]+3;
  rotate[2] = rotate[1]+3;
  return rotate;
}

const static double PHE_ARO_HDR = 2.30;
const static double TRP_ARO_HDR = 2.30;
const static double TYR_ARO_HDR = 2.30;

class pseudoAA{
 protected:
  atom* N;
  atom* CA;
  atom* C;
  atom* O;
  atom* H;
  atom* SideChains;
  int nSideChains;
  /*-----------------------------------------
    0, no side chain atoms (GLY)
    1, one CB side chain atom (ALA)
    2, CB + CG (15)
    3, CB + CG1 + CG2 (3) beta-branched
    -----------------------------------------*/
  int isBetaBranched;
  /*------------------------------
    1(true) for ILE, THR, VAL
    0 for others
    ------------------------------*/
  aa_t id;
  double phi;
  double psi;
  /*------------------------------------
    non INF value if having definitions
    range: [-PI, PI)
    ------------------------------------*/
  double chi; /*range: [0, 2PI)*/
 public:
  pseudoAA(const char* name){
    N  = new atom('N', N_CA*cos(init_N_CA_x*rpi), -N_CA*sin(init_N_CA_x*rpi), 0);
    H  = NULL;
    CA = new atom('C', 0, 0 ,0);
    C  = new atom('C', CA_C*cos(init_C_CA_x*rpi), CA_C*sin(init_C_CA_x*rpi), 0);
    double theta = init_C_CA_x+CA_C_O-180.0;
    O  = new atom('O', C->getR()->getX() + C_O*cos(theta*rpi), 
		  C->getR()->getY() + C_O*sin(theta*rpi), 0);
    id = aa_Name2Type(name);
    const double* G1=G1_CONST[id];
    const double* G2=NULL;
    char G1_name='C', G2_name='C';
    nSideChains=2;
    isBetaBranched=0;
    switch(id){
    case GLY:
      nSideChains=0;
      G1 = NULL;
      break;
    case ALA:
      nSideChains=1;
      G1 = NULL;
      break;
    case ILE:
      G2 = ILE_G2;
      nSideChains=3;
      isBetaBranched=1;
      break;
    case THR:
      G2 = THR_G2;
      G2_name = 'O';
      nSideChains=3;
      isBetaBranched=1;
      break;
    case VAL:
      G2 = VAL_G2;
      nSideChains = 3;
      isBetaBranched=1;
      break;
    case CYS:
    case MET:
      G1_name = 'S';
      break;
    case SER:
      G1_name = 'O';
      break;
    default:
      break;
    }
    if(nSideChains>0){
      SideChains = new atom[nSideChains]('C',0,0,0);
      double theta, phi, x, y, z;
      //CB -- for non-GLY amino acids
      if(id==PRO){
	theta = (180-CB_N_CA_C)*rpi;
	phi = (360.0-PRO_N_CA_CBXY-init_N_CA_x)*rpi;
	x = CA_CB*cos(theta)*cos(phi);
	y = CA_CB*cos(theta)*sin(phi);
	z = CA_CB*sin(theta);
	SideChains[0]=atom('C',x,y,z);
      }
      else{
	theta = (180-CB_N_CA_C)*rpi;
	phi = (360.0-N_CA_CBXY-init_N_CA_x)*rpi;
	x = CA_CB*cos(theta)*cos(phi);
	y = CA_CB*cos(theta)*sin(phi);
	z = CA_CB*sin(theta);
	SideChains[0]=atom('C',x,y,z);
      }
      //CG1 -- for non-ALA amino acids
      if(nSideChains>1 && G1){
	double G_CA = G1[G1_CA];
	double G_CB = G1[G1_CB];
	double ftmp = CA_CB*CA_CB + G_CB*G_CB - G_CA*G_CA;
	ftmp /= (2.0*CA_CB*G_CB);
	double ca_cb_g = acos(ftmp);
	double n_ca_cb = N_CA_CB*rpi;
	if(id==PRO) n_ca_cb = PRO_N_CA_CB*rpi;
	vec n_ca = *CA->getR() - *N->getR();
	vec ca_cb = *SideChains[0].getR() - *CA->getR();
	n_ca.Normalize();
	ca_cb.Normalize();
	double d_ca_cb;
	double d_n_ca;
	if(ca_cb_g>n_ca_cb){
	  d_ca_cb = G_CB*sin(ca_cb_g-n_ca_cb)/sin(n_ca_cb);
	  d_n_ca  = G_CB*sin(PI-ca_cb_g)/sin(n_ca_cb);
	}
	else if(ca_cb_g==n_ca_cb){
	  d_ca_cb = 0;
	  d_n_ca  = G_CB;
	}
	else{
	  d_ca_cb = -G_CB*sin(n_ca_cb-ca_cb_g)/sin(PI-n_ca_cb);
	  d_n_ca  =  G_CB*sin(ca_cb_g)/sin(PI-n_ca_cb);
	}
	vec vtmp = d_n_ca*n_ca + d_ca_cb*ca_cb + *SideChains[0].getR();
	SideChains[1] = atom(G1_name, vtmp.getX(), vtmp.getY(), vtmp.getZ());
	//PROLINE is special because the N_CG is fixed!
	if(id==PRO){
	  double a = N_CA*sin(n_ca_cb);
	  double b = G_CB*sin(ca_cb_g);
	  double c = CA_CB - N_CA*cos(n_ca_cb) - G_CB*cos(ca_cb_g);
	  double d = PRO_N_CG;
	  ftmp = d*d - a*a -b*b -c*c;
	  ftmp /= (2.0*a*b);
	  ftmp = acos(ftmp);
	  SideChains[1].getR()->Rotate(*CA->getR(), *SideChains[0].getR(), ftmp);
	  chi = ftmp;
	}
	//CG2 -- for ILE, THR, VAL
	if(nSideChains==3 && G2){
	  G_CA = G2[G2_CA];
	  G_CB = G2[G2_CB];
	  ftmp = CA_CB*CA_CB + G_CB*G_CB - G_CA*G_CA;
	  ftmp /= (2.0*CA_CB*G_CB);
	  ca_cb_g = acos(ftmp);
	  if(ca_cb_g>n_ca_cb){
	    d_ca_cb = G_CB*sin(ca_cb_g-n_ca_cb)/sin(n_ca_cb);
	    d_n_ca  = G_CB*sin(PI-ca_cb_g)/sin(n_ca_cb);
	  }
	  else if(ca_cb_g==n_ca_cb){
	    d_ca_cb = 0;
	    d_n_ca  = G_CB;
	  }
	  else{
	    d_ca_cb = -G_CB*sin(n_ca_cb-ca_cb_g)/sin(PI-n_ca_cb);
	    d_n_ca  =  G_CB*sin(ca_cb_g)/sin(PI-n_ca_cb);
	  }
	  vtmp = d_n_ca*n_ca + d_ca_cb*ca_cb + *SideChains[0].getR();
	  SideChains[2] = atom(G2_name, vtmp.getX(), vtmp.getY(), vtmp.getZ());
	  //Rotate as CA to CB axis
	  double theta = G2[G2_CA_CB_G1]*rpi;
	  SideChains[2].getR()->Rotate(*CA->getR(), *SideChains[0].getR(), theta);
	}
      }
    }
    else {//GLY -- no side chains
      SideChains = NULL;
    }
    phi = INF;
    psi = INF;
    if(id!=PRO) chi = 0;
  }

  pseudoAA(aa& theAA){
    N = new atom(*theAA.getN());
    H = NULL;
    CA= new atom(*theAA.getAtomAt(0));
    C = new atom(*theAA.getC());
    O = new atom(*theAA.getO());
    id = aa_Name2Type(cstr(theAA.getName()).getStr());
    const double* G1=G1_CONST[id];
    const double* G2=NULL;
    char G1_name='C', G2_name='C';
    isBetaBranched=0;
    nSideChains=2;
    switch(id){
    case GLY:
      nSideChains=0;
      G1 = NULL;
      break;
    case ALA:
      nSideChains=1;
      G1 = NULL;
      break;
    case ILE:
      G2 = ILE_G2;
      nSideChains=3;
      isBetaBranched=1;
      break;
    case THR:
      G2 = THR_G2;
      G2_name = 'O';
      nSideChains=3;
      isBetaBranched=1;
      break;
    case VAL:
      G2 = VAL_G2;
      nSideChains = 3;
      isBetaBranched=1;
      break;
    case CYS:
    case MET:
      G1_name = 'S';
      break;
    case SER:
      G1_name = 'O';
      break;
    default:
      break;
    }
    if(nSideChains>0){
      SideChains = new atom[nSideChains]('C',0,0,0);
      double theta, phi, x, y, z;
      //CB -- for non-GLY amino acids
      SideChains[0]=atom(*theAA.getAtomAt(1));
      //CG1 -- for non-ALA amino acids
      if(nSideChains>1 && G1){
	double G_CA = G1[G1_CA];
	double G_CB = G1[G1_CB];
	double ftmp = CA_CB*CA_CB + G_CB*G_CB - G_CA*G_CA;
	ftmp /= (2.0*CA_CB*G_CB);
	double ca_cb_g = acos(ftmp);
	double n_ca_cb = N_CA_CB*rpi;
	if(id==PRO) n_ca_cb = PRO_N_CA_CB*rpi;
	vec n_ca = *CA->getR() - *N->getR();
	vec ca_cb = *SideChains[0].getR() - *CA->getR();
	n_ca.Normalize();
	ca_cb.Normalize();
	double d_ca_cb;
	double d_n_ca;
	if(ca_cb_g>n_ca_cb){
	  d_ca_cb = G_CB*sin(ca_cb_g-n_ca_cb)/sin(n_ca_cb);
	  d_n_ca  = G_CB*sin(PI-ca_cb_g)/sin(n_ca_cb);
	}
	else if(ca_cb_g==n_ca_cb){
	  d_ca_cb = 0;
	  d_n_ca  = G_CB;
	}
	else{
	  d_ca_cb = -G_CB*sin(n_ca_cb-ca_cb_g)/sin(PI-n_ca_cb);
	  d_n_ca  =  G_CB*sin(ca_cb_g)/sin(PI-n_ca_cb);
	}
	vec vtmp = d_n_ca*n_ca + d_ca_cb*ca_cb + *SideChains[0].getR();
	SideChains[1] = atom(G1_name, vtmp.getX(), vtmp.getY(), vtmp.getZ());
	//PROLINE is special because the N_CG is fixed!
	if(id==PRO){
	  double a = N_CA*sin(n_ca_cb);
	  double b = G_CB*sin(ca_cb_g);
	  double c = CA_CB - N_CA*cos(n_ca_cb) - G_CB*cos(ca_cb_g);
	  double d = PRO_N_CG;
	  ftmp = d*d - a*a -b*b -c*c;
	  ftmp /= (2.0*a*b);
	  ftmp = acos(ftmp);
	  SideChains[1].getR()->Rotate(*CA->getR(), *SideChains[0].getR(), ftmp);
	  chi = ftmp;
	}
	//CG2 -- for ILE, THR, VAL
	if(nSideChains==3 && G2){
	  G_CA = G2[G2_CA];
	  G_CB = G2[G2_CB];
	  ftmp = CA_CB*CA_CB + G_CB*G_CB - G_CA*G_CA;
	  ftmp /= (2.0*CA_CB*G_CB);
	  ca_cb_g = acos(ftmp);
	  if(ca_cb_g>n_ca_cb){
	    d_ca_cb = G_CB*sin(ca_cb_g-n_ca_cb)/sin(n_ca_cb);
	    d_n_ca  = G_CB*sin(PI-ca_cb_g)/sin(n_ca_cb);
	  }
	  else if(ca_cb_g==n_ca_cb){
	    d_ca_cb = 0;
	    d_n_ca  = G_CB;
	  }
	  else{
	    d_ca_cb = -G_CB*sin(n_ca_cb-ca_cb_g)/sin(PI-n_ca_cb);
	    d_n_ca  =  G_CB*sin(ca_cb_g)/sin(PI-n_ca_cb);
	  }
	  vtmp = d_n_ca*n_ca + d_ca_cb*ca_cb + *SideChains[0].getR();
	  SideChains[2] = atom(G2_name, vtmp.getX(), vtmp.getY(), vtmp.getZ());
	  //Rotate as CA to CB axis
	  double theta = G2[G2_CA_CB_G1]*rpi;
	  SideChains[2].getR()->Rotate(*CA->getR(), *SideChains[0].getR(), theta);
	}
      }
    }
    else {//GLY -- no side chains
      SideChains = NULL;
    }
    phi = INF;
    psi = INF;
    if(id!=PRO) chi = 0;
  }
  
  pseudoAA(aa& theAA, pseudo_aa& thePseudoAA){
    if(theAA.getN()) N = new atom(*theAA.getN());
    else N = NULL;
    H = NULL;
    if(theAA.getC()) C = new atom(*theAA.getC());
    else C = NULL;
    if(theAA.getO()) O = new atom(*theAA.getO());
    else O = NULL;
    CA= new atom(*theAA.getAtomAt(0));
    id = aa_Name2Type(cstr(theAA.getName()).getStr());
    char G1_name='C', G2_name='C';
    isBetaBranched=0;
    vec G1;
    vec G2;
    switch(id){
    case GLY:
      nSideChains=0;
      break;
    case ALA:
      nSideChains=1;
      break;
    case ILE:
      nSideChains=3;
      isBetaBranched=1;
      G1 = *theAA.getAtomAt(3)->getR();
      G2 = *theAA.getAtomAt(1)->getR() + (*theAA.getAtomAt(2)->getR() - 
					  *theAA.getAtomAt(1)->getR()).norm()*ILE_G2[G2_CB];
      break;
    case THR:
      nSideChains=3;
      isBetaBranched=1;
      G1 = *theAA.getAtomAt(3)->getR(); //CG2
      G2_name = 'O';
      G2 = *theAA.getAtomAt(2)->getR(); //OG1
      break;
    case VAL:
      nSideChains=3;
      isBetaBranched=1;
      G1 = *theAA.getAtomAt(2)->getR();
      G2 = *theAA.getAtomAt(3)->getR();
      break;
    case PRO:
      nSideChains=2;
      G1 = (*theAA.getAtomAt(2)->getR() + 
	    *theAA.getAtomAt(3)->getR())/2.0;
      break;
    case LYS:
      nSideChains=2;
      G1 = *theAA.getAtomAt(4)->getR();
      break;
    case ARG:
      nSideChains=2;
      G1 = *theAA.getAtomAt(5)->getR();
      break;
    case MET:
      nSideChains=2;
      G1_name = 'S';
      G1 = *theAA.getAtomAt(3)->getR();
      break;
    case CYS:
      nSideChains=2;
      G1_name = 'S';
      G1 = *theAA.getAtomAt(2)->getR();
      break;
    case SER:
      nSideChains=2;
      G1_name = 'O';
      G1 = *theAA.getAtomAt(2)->getR();
      break;
    default:
      nSideChains=2;
      G1 = *thePseudoAA.getSide();
      break;
    }
    if(nSideChains>0){
      SideChains = new atom[nSideChains]('C',0,0,0);
      SideChains[0] = atom(*theAA.getAtomAt(1));
      if(nSideChains>1){
	SideChains[1] = atom(G1_name, G1.getX(), G1.getY(), G1.getZ());
	if(isBetaBranched){
	  SideChains[2] = atom(G2_name, G2.getX(), G2.getY(), G2.getZ());
	}
      }
    }
    else {//GLY -- no side chains
      SideChains = NULL;
    }
    phi=INF;
    psi=INF;
    if(N && CA && nSideChains>=2){
      vec* n  = N->getR();
      vec* ca = CA->getR();
      vec* cb = SideChains[0].getR();
      vec* cg = SideChains[1].getR();
      vec cg_cb = *cb - *cg;
      vec cb_ca = *ca - *cb;
      vec ca_n  = *n  - *ca;
      double mod2 = cb_ca*cb_ca;
      vec r1 = cg_cb - (cg_cb*cb_ca)*cb_ca/mod2;
      vec r2 = ca_n  - (cb_ca*ca_n )*cb_ca/mod2;
      double cos = r1*r2;
      cos /= sqrt((r1*r1)*(r2*r2));
      double coef = (cb_ca^r1)*r2;
      if(coef>0) chi = acos(cos);
      else chi=PI*2 - acos(cos);
    }
    else{
      chi=INF;
    }
  }

  ~pseudoAA(){
    delete N;
    delete CA;
    delete C;
    delete O;
    delete H;
    delete [] SideChains;
  }
  
  void intrinsicAxis(vec& x, vec& y, vec& z){
    //return the intric axises
    vec ca_n = *N->getR() - *CA->getR();
    ca_n.Normalize();
    vec ca_c = *C->getR() - *CA->getR();
    ca_c.Normalize();
    double n_ca_c = N_CA_C*rpi;
    double n_ca_x = init_N_CA_x*rpi;
    double c_ca_x = init_C_CA_x*rpi;
    //intrinsic x axis
    double dist_ca_n = sin(c_ca_x)/sin(n_ca_c);
    double dist_ca_c = sin(n_ca_x)/sin(n_ca_c);
    x = ca_n*dist_ca_n + ca_c*dist_ca_c;
    x.Normalize();
    //intrinsic z axis
    z = ca_n^ca_c;
    z.Normalize();
    //derive the intrinsic y axis = z^x
    y = z^x;
  }
  
  const vec nextCA(vec& xbar, vec& ybar, vec& zbar){
    this->intrinsicAxis(xbar, ybar, zbar);
    double theta = (C_CA0_CA + init_C_CA_x)*rpi;
    return *this->getCA()->getR()
      +CA_CA*cos(theta)*xbar
      +CA_CA*sin(theta)*ybar;
  }
  
  //rotate the whole amino acids
  void rotate(const vec& p1, const vec& p2, const double theta){
    N->getR()->Rotate(p1, p2, theta);
    CA->getR()->Rotate(p1, p2, theta);
    C->getR()->Rotate(p1, p2, theta);
    O->getR()->Rotate(p1, p2, theta);
    for(int i=0; i<nSideChains; i++) 
      SideChains[i].getR()->Rotate(p1, p2, theta);
  }
  
  //fix N, CA; rotate C,O,SideChains
  void rotatePhi(double dphi){
    if(phi!=INF){
      C->getR()->Rotate(*N->getR(), *CA->getR(), dphi);
      O->getR()->Rotate(*N->getR(), *CA->getR(), dphi);
      for(int i=0; i<nSideChains; i++)
	SideChains[i].getR()->Rotate(*N->getR(), *CA->getR(), dphi);
      phi += dphi/rpi;
      if(phi >= PI)phi -= PI*2.0;
    }
  }
  
  //fix CA, C, O; rotate N,SideChains
  void rotatePsi(double dpsi){
    if(psi!=INF){
      N->getR()->Rotate(*C->getR(), *CA->getR(), dpsi);
      for(int i=0; i<nSideChains; i++)
	SideChains[i].getR()->Rotate(*C->getR(), *CA->getR(), dpsi);
      psi += dpsi/rpi;
      if(psi >= PI)psi -= PI*2.0;
    }
  }
  
  //fix CA, CB, rotate CG or G2
  void rotateChi(double dchi){
    if(nSideChains>1){
      SideChains[1].getR()->Rotate(*getCA()->getR(), *getCB()->getR(), dchi);
      if(isBetaBranched){
	SideChains[2].getR()->Rotate(*getCA()->getR(), *getCB()->getR(), dchi);
      }
      chi += dchi;
      if(chi >= 2.0*PI) chi -= PI*2.0;
    }
  }
  
  int updateChi(double& out){
    if(id!=GLY && id!=ALA && 
       getN() && getCA() && getCB() && getG()){
      vec* Nr = getN()->getR();
      vec* CAr= getCA()->getR();
      vec* CBr= getCB()->getR();
      vec* CGr= getG()->getR();
      vec CG_CB = *CBr - *CGr;
      vec CB_CA = *CAr - *CBr;
      vec CA_N  = *Nr  - *CAr;
      double mod2 = CB_CA*CB_CA;
      vec r1 = CG_CB - (CG_CB*CB_CA)*CB_CA/mod2;
      vec r2 = CA_N  - (CB_CA*CA_N )*CB_CA/mod2;
      double cos = r1*r2;
      cos /= sqrt((r1*r1)*(r2*r2));
      double coef = (CB_CA^r1)*r2;
      if(coef>0) chi = acos(cos);
      else chi = PI*2 - acos(cos);
      out = chi;
      return 1;
    }
    return 0;
  }

  void joinTo(pseudoAA& prev, double prev_psi=-180, double phi=-180){
    static double** rotate = newRMatrix();
    static vec xbar, ybar, zbar;
    vec shift = prev.nextCA(xbar, ybar, zbar);
    
    xbar *= -1.0;
    zbar *= -1.0;
    rotate[0][0] = xbar.getX();
    rotate[0][1] = ybar.getX();
    rotate[0][2] = zbar.getX();

    rotate[1][0] = xbar.getY();
    rotate[1][1] = ybar.getY();
    rotate[1][2] = zbar.getY();

    rotate[2][0] = xbar.getZ();
    rotate[2][1] = ybar.getZ();
    rotate[2][2] = zbar.getZ();
    N->getR()->EularRotate(rotate); 
    *N->getR() += shift;
    CA->getR()->EularRotate(rotate);
    *CA->getR() += shift;
    C->getR()->EularRotate(rotate);
    *C->getR() += shift;
    O->getR()->EularRotate(rotate);
    *O->getR() += shift;
    for(int i=0; i<nSideChains; i++) {
      SideChains[i].getR()->EularRotate(rotate);
      *SideChains[i].getR() += shift;
    }
    prev.setPsi(-180);
    this->setPhi(-180);
    if(prev_psi>-180){
      double dpsi = prev_psi+180;
      dpsi *= rpi;
      this->rotate(*prev.getCA()->getR(),
		   *prev.getC()->getR(), dpsi);
      prev.getO()->getR()->Rotate(*prev.getCA()->getR(),
				  *prev.getC()->getR(), dpsi);
      prev.setPsi(prev_psi);
    }
    if(phi>-180){
      double dphi = phi+180;
      dphi *= rpi;
      rotatePhi(dphi);
    }
  }
  
  atom* getN(){return N;}
  
  atom* getCA(){return CA;}
  
  atom* getC(){return C;}
  
  atom* getO(){return O;}

  atom* getH(){return H;}

  void createH(vec& hvec){
    H = new atom('H',hvec);
  }

  atom* getSideChains(){return SideChains;}
  
  int getNSideChains(){return nSideChains;}
  
  atom* getCB(){
    if(id!=GLY)return SideChains;
    else return NULL;
  }
  
  atom* getG(){
    if(id!=GLY && id!=ALA) return SideChains+1;
    else return NULL;
  }
  
  atom* getG2(){
    if(isBetaBranched) return SideChains+2;
    else return NULL;
  }

  void shift(vec& r){
    if(getN()) *N->getR() += r;
    if(getC()) *C->getR() += r;
    if(getO()) *O->getR() += r;
    if(getCA()) *CA->getR() += r;
    for(int i=0; i<nSideChains; i++)
      *SideChains[i].getR() += r;
  }
  
  double getPhi(){return phi;}
  
  void setPhi(double x){ phi = x;}
  
  double getPsi(){return psi;}
  
  void setPsi(double x){psi = x;}
  
  int getID(){return id;}
  
  const char* getName(){
    return cstr_aa3[id];
  }

  void printPDB(ostream& out, int resIndex, int& atomIndex){
    double x,y,z;
    const char* name;
    char field[100];
    name = cstr_aa3[id];
    x = N->getR()->getX();
    y = N->getR()->getY();
    z = N->getR()->getZ();
    sprintf(&field[0],"ATOM  %5d  N   %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00              ",
	    atomIndex++, name, resIndex, x, y, z);
    out << field << endl;
    x = C->getR()->getX();
    y = C->getR()->getY();
    z = C->getR()->getZ();
    sprintf(&field[0],"ATOM  %5d  C   %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00              ",
	    atomIndex++, name, resIndex, x, y, z);
    out << field << endl;
    if(O){
      x = O->getR()->getX();
      y = O->getR()->getY();
      z = O->getR()->getZ();
      sprintf(&field[0],"ATOM  %5d  O   %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00              ",
	      atomIndex++, name, resIndex, x, y, z);
      out << field << endl;
    }
    x = CA->getR()->getX();
    y = CA->getR()->getY();
    z = CA->getR()->getZ();
    sprintf(&field[0],"ATOM  %5d  CA  %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00              ",
	    atomIndex++, name, resIndex, x, y, z);
    out << field << endl;
    if(nSideChains>=1){
      x = getCB()->getR()->getX();
      y = getCB()->getR()->getY();
      z = getCB()->getR()->getZ();      
      sprintf(&field[0],"ATOM  %5d  CB  %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00              ",
	      atomIndex++, name, resIndex, x, y, z);
      out << field << endl;
      if(nSideChains>=2){
	x = getG()->getR()->getX();
	y = getG()->getR()->getY();
	z = getG()->getR()->getZ();  
	char* atomName = new char[10];
	atomName="CG";
	if(id==MET || id==CYS) atomName="SG";
	else if(id==SER) atomName="OG";
	if(isBetaBranched) atomName="CG1";
	sprintf(&field[0],"ATOM  %5d  %-3s %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00              ",
		atomIndex++, atomName, name, resIndex, x, y, z);
	out << field << endl;
	if(nSideChains==3){
	  x = getG2()->getR()->getX();
	  y = getG2()->getR()->getY();
	  z = getG2()->getR()->getZ(); 
	  atomName="CG2";
	  if(id==THR) atomName="OG2";
	  sprintf(&field[0],"ATOM  %5d  %-3s %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00              ",
		  atomIndex++, atomName, name, resIndex, x, y, z);
	  out << field << endl;
	}
      }
    }
    if(getH()){
      x = getH()->getR()->getX();
      y = getH()->getR()->getY();
      z = getH()->getR()->getZ();
      sprintf(&field[0],"ATOM  %5d  H   %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00              ",
	      atomIndex++, name, resIndex, x, y, z);
      out << field << endl;
    }
  }
};

#endif /*_PSEUDO_AA_H_*/
