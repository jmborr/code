#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include "PDBLib.h"
#include "relax_txt.h"
#include "relax_aa.h"
#include "relax_pep.h"
#include "random.h"

using namespace std;

double estep=10.0;
double egap=100.0;
void type_sprint(char* buf, dmd_atom_t type, double mass, double hdr, double ir){
  if(ir==0) ir = hdr + eps;
  sprintf(buf, "%ld %lf %lf %lf", type, mass, hdr, ir);
}

void el_col_sprint(char* buf, dmd_atom_t a, dmd_atom_t b, double min_d){
  double min=0.1;
  int nstep=10;
  double step=(min_d-min)/nstep;
  int pt=sprintf(buf,"%ld %ld %lf ", a, b, min);
  for(int i=1; i<nstep; i++){
    pt+=sprintf(&buf[pt], "%lf %lf ", min+i*step, estep);
  }
  sprintf(&buf[pt],"%lf %lf", min_d, egap);
}

void single_bond_sprint(char* buf, double ave, double dev){
  double min=0.1;
  double step=(ave*(1.0-dev)-min)/6;
  int pt=sprintf(buf,"%lf ", min);
  for(int i=1; i<5; i++){
    pt+=sprintf(&buf[pt], "%lf %lf ", min+i*dev*ave, estep);
  }
  pt+=sprintf(&buf[pt], "%lf %lf %lf %lf ", 
	      ave*(1.0-dev), egap, ave*(1.0+dev), -egap);
  double max=ave*(1.0+6*dev);
  for(int i=2; i<6; i++){
    pt+=sprintf(&buf[pt], "%lf %lf ", ave*(1.0+dev*i), -estep);
  }
  sprintf(&buf[pt], "%lf %lf %lf %lf %lf", 
	  max, -egap, max+0.05,  2*egap+4.0*estep, max+0.1);
}

void single_bond_sprint_len(char* buf, double ave, double len_dev){
  single_bond_sprint(buf, ave, len_dev/ave);
}

void double_bond_sprint(char* buf, double l_ave, double l_dev, 
			double r_ave, double r_dev){
  if(l_ave*(1.0+l_dev) > r_ave*(1.0-r_dev)){
    cerr << "Error: overlaping bonds" << endl;
    exit(1);
  }
  double min=l_ave*(1.0-6*l_dev);
  int pt=sprintf(buf,"%lf ", min);
  for(int i=1; i<5; i++){
    pt+=sprintf(&buf[pt], "%lf %lf ", min+i*l_dev*l_ave, estep);
  }
  pt+=sprintf(&buf[pt], "%lf %lf %lf %lf ", l_ave*(1.0-l_dev), egap,
	      l_ave*(1.0+l_dev), -egap);
  
  double step=(r_ave*(1.0-r_dev)-l_ave*(1.0+l_dev))/11.0;
  for(int i=1; i<5; i++){
    pt+=sprintf(&buf[pt], "%lf %lf ", 
		l_ave*(1.0+l_dev)+i*step, -estep);
  }
  pt+=sprintf(&buf[pt], "%lf %lf %lf %lf ",
	      l_ave*(1.0+l_dev)+5*step, -10*egap,
	      l_ave*(1.0+l_dev)+6*step, 10*egap);
  for(int i=7; i<11; i++){
    pt+=sprintf(&buf[pt], "%lf %lf ", 
		l_ave*(1.0+l_dev)+i*step, estep);
  }
  pt+=sprintf(&buf[pt], "%lf %lf %lf %lf ", r_ave*(1.0-r_dev), egap,
	      r_ave*(1.0+r_dev),-egap);
  double max=r_ave*(1.0+6.0*r_dev);
  for(int i=2; i<6; i++){
    pt+=sprintf(&buf[pt], "%lf %lf ", r_ave*(1.0+i*r_dev), -estep);
  }
  sprintf(&buf[pt], "%lf %lf %lf %lf %lf", 
	  max, -egap, max+0.05, 2*egap+4.0*estep, max+0.1);
}

int side_sprint(char* buf, atom* cb1, int t1, atom* cb2, int t2, double hdr, double ir){
  double dist = cb1->getDist(*cb2);
  if(dist<hdr){
    int nstep=5;
    int hstep=10;
    double min=dist-0.1;
    double step = (hdr-min)/nstep;
    int pt=sprintf(buf, "%ld %ld %lf ", t1, t2, min);
    for(int i=1; i<nstep; i++){
      pt+=sprintf(&buf[pt], "%lf %lf ", min+i*step, hstep);
    }
    pt+=sprintf(&buf[pt], "%lf 100.0 %lf", hdr, ir);
  }
  else if(dist<ir){
    sprintf(buf, "%ld %ld %lf %lf", t1, t2, hdr, ir);
  }
  else return 0;
  return 1;
}

double getTriDev_angle(double a, double b, double dev, double theta, double& c){
  c = sqrt(a*a+b*b-2.0*a*b*cos(theta));
  double dummy = (a-b*cos(theta))/c*a;
  double dummy1= (b-a*cos(theta))/c*b;
  dummy = sqrt(dummy*dummy + dummy1*dummy1)*dev;
  dummy1 = dummy/c;
  return dummy1;
}

double getTriDev_length(double a, double b, double dev, double c){
  double theta = acos((a*a+b*b-c*c)/(2.0*a*b));
  double dummy = (a-b*cos(theta))/c*a;
  double dummy1= (b-a*cos(theta))/c*b;
  dummy = sqrt(dummy*dummy + dummy1*dummy1)*dev;
  dummy1 = dummy/c;
  return dummy1;
}

void list_atom_sprint(char* buf, int index, dmd_atom_t type, vec& r, vec& v){
  sprintf(buf,"%-4ld %-4ld %18.13lf %18.13lf %18.13lf %18.13lf %18.13lf %18.13lf", 
	  index, type, r.getX(), r.getY(), r.getZ(), 
	  v.getX(), v.getY(), v.getZ());
}

vec getCM(pseudoPep& p){
  vec tmp(0,0,0);
  int natom=0;
  for(int i=0; i<p.getLength(); i++){
    pseudoAA* theAA = p.getResidue(i);
    tmp += *theAA->getN()->getR(); natom++;
    tmp += *theAA->getC()->getR(); natom++;
    tmp += *theAA->getO()->getR(); natom++;
    tmp += *theAA->getCA()->getR(); natom++;
    if(theAA->getCB()) {
      tmp += *theAA->getCB()->getR(); 
      natom++;
      if(theAA->getG()){
	tmp += *theAA->getG()->getR();
	natom++;
	if(theAA->getG2()){
	  tmp += *theAA->getG2()->getR();
	  natom++;
	}
      }
    }
  }
  tmp /= static_cast<double>(natom);
  return tmp;
}

void printSYS_SIZE(ostream& out, double size, int dimension=3){
  out << txtKeyWords[SYS_SIZE] << endl;
  for(int i=0; i<dimension; i++){
    out << size << " ";
  }
  out << endl;
}

void printNUM_ATOMS(ostream& out, pseudoPep& p, 
		    int np=1, int ng=0){
  out << txtKeyWords[NUM_ATOMS] << endl;
  int nTotal = ng;
  pseudoAA* lastAA = p.getResidue(p.getLength()-1);
  atom* lastAtom = NULL;
  if(lastAtom=lastAA->getG2());
  else if(lastAtom=lastAA->getG());
  else if(lastAtom=lastAA->getCB());
  else if(lastAtom=lastAA->getCA());
  else{
    cerr << "last amino acids is not a protein!!" << endl;
    exit(1);
  }
  nTotal += np*lastAtom->getIndex();
  out << nTotal << endl;
}

int printATOM_TYPE(ostream& out, pseudoPep& p){
  out << txtKeyWords[TYPE_ATOMS] << endl;
  char buf[1024];
  /*generic atoms*/
  /*N*/
  type_sprint(buf, _N_, N_M, N_HDR, 0.0);
  out << buf << endl;
  /*C*/
  type_sprint(buf, _C_, C_M, C_HDR, 0.0);
  out << buf << endl;
  /*O*/
  type_sprint(buf, _O_, O_M, O_HDR, 0.0);
  out << buf << endl;
  /*CA*/
  type_sprint(buf, _CA_, CA_M, CA_HDR, 0.0);
  out << buf << endl;
  /*CB*/
  type_sprint(buf, _CB_, CB_M, CB_HDR, 0.0);
  out << buf << endl;
  /*PRO_N*/
  type_sprint(buf, _PRO_N_, N_M, N_HDR, 0.0);
  out << buf << endl;
  /*NP*/
  type_sprint(buf, _N_HB_, N_M, N_HDR, 0.0);
  out << buf << endl;
  /*OP*/
  type_sprint(buf, _O_HB_, O_M, O_HDR, 0.0);
  out << buf << endl;
  
  int type_index= _O_HB_;
  for(int i=0;i<p.getLength();i++){
    aa_t type=static_cast<aa_t>(p.getResidue(i)->getID());
    if(type!=GLY){
      //CB
      type_index++;
      type_sprint(buf, static_cast<dmd_atom_t>(type_index), CB_M, CB_HDR, 0.0);
      out << buf << endl;
      //CG1
      if(type!=ALA){
	type_index++;
	const double* cg_para = G1_CONST[type];
	type_sprint(buf, static_cast<dmd_atom_t>(type_index), 
		    cg_para[G1_M], cg_para[G1_HDR], cg_para[G1_IR]);
	out << buf << endl;
	//CG2
	if(type==ILE){
	  type_index++;
	  type_sprint(buf, static_cast<dmd_atom_t>(type_index), ILE_G2[G2_M], ILE_G2[G2_HDR], ILE_G2[G2_IR]);
	  out << buf << endl;
	}
	else if(type==THR){
	  type_index++;
	  type_sprint(buf, static_cast<dmd_atom_t>(type_index), THR_G2[G2_M], THR_G2[G2_HDR], THR_G2[G2_IR]);
	  out << buf << endl;
	}
	else if(type==VAL){
	  type_index++;
	  type_sprint(buf, static_cast<dmd_atom_t>(type_index), VAL_G2[G2_M], VAL_G2[G2_HDR], VAL_G2[G2_IR]);
	  out << buf << endl;
	}
      }
    }
  }
  return type_index;
}

void printNONEL_COL(ostream& out, pseudoPep& p, int* cbIndex){
  out << txtKeyWords[NONEL_COL] << endl;
  char buf[1024];
  /*mainchain hydogen related*/
  /*N-O*/
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _N_, _O_, min_NO, HB_N_O[1]);
  out << buf << endl;
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _N_, _O_HB_, min_NO, HB_N_O[1]);
  out << buf << endl; 
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _N_HB_, _O_HB_, min_NO, HB_N_O[1]);
  out << buf << endl;
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _N_HB_, _O_, min_NO, HB_N_O[1]);
  out << buf << endl;
  /*N-C*/
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _N_, _C_, min_NC, HB_N_C[3]);
  out << buf << endl;
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _N_HB_, _C_, min_NC, HB_N_C[3]);
  out << buf << endl;

  /*C-O*/
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _O_, _C_, min_CO, HB_O_C[2]);
  out << buf << endl;
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _O_HB_, _C_, min_CO, HB_O_C[2]);
  out << buf << endl;  

  /*Ca-O*/
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _O_, _CA_, min_OCa, HB_O_CA[2]);
  out << buf << endl;
  sprintf(buf, "%ld %ld 1.0 %lf 10.0 %lf 0.0000", _O_HB_, _CA_, min_OCa, HB_O_CA[2]);
  out << buf << endl;

  /*---extra backbone--*/
  /*N-N*/
  el_col_sprint(buf, _N_,     _N_,     min_NN);
  out << buf << endl;
  el_col_sprint(buf, _N_,     _N_HB_,  min_NN);
  out << buf << endl;
  el_col_sprint(buf, _N_,     _PRO_N_, min_NN);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_,  _N_HB_,  min_NN);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_,  _PRO_N_, min_NN);
  out << buf << endl;
  el_col_sprint(buf, _PRO_N_, _PRO_N_, min_NN);
  out << buf << endl;
  /*N-C*/
  el_col_sprint(buf, _PRO_N_, _C_,     min_NC);
  out << buf << endl;
  /*N-O*/
  el_col_sprint(buf, _PRO_N_, _O_,     min_NO);
  out << buf << endl;
  el_col_sprint(buf, _PRO_N_, _O_HB_,  min_NO);
  out << buf << endl;
  /*N-CA*/
  el_col_sprint(buf, _N_,     _CA_,    N_HDR+CA_HDR);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_,  _CA_,    N_HDR+CA_HDR);
  out << buf << endl;
  el_col_sprint(buf, _PRO_N_, _CA_,    N_HDR+CA_HDR);
  out << buf << endl;
  
  /*C-C*/
  el_col_sprint(buf, _C_,    _C_,      C_HDR+C_HDR);
  out << buf << endl;
  /*C-Ca*/
  el_col_sprint(buf, _C_,    _CA_,     C_HDR+CA_HDR);
  out << buf << endl;

  /*Ca-Ca*/
  el_col_sprint(buf, _CA_,   _CA_,     CA_HDR+CA_HDR);
  out << buf << endl;
  
  /*O-O*/
  el_col_sprint(buf, _O_,    _O_,      O_HDR+O_HDR);
  out << buf << endl;
  el_col_sprint(buf, _O_,    _O_HB_,   O_HDR+O_HDR);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _O_HB_,   O_HDR+O_HDR);
  out << buf << endl;
  
  /*indexed side chain related*/
  for(int i=0; i<p.getLength(); i++){
    aa_t type=static_cast<aa_t>(p.getResidue(i)->getID());
    int index=cbIndex[i]+_O_HB_+1;
    atom* CB=NULL;
    atom* CG=NULL;
    atom* CG2=NULL;
    double cb_hdr, cb_ir;
    double cg_hdr, cg_ir;
    double cg2_hdr, cg2_ir;
    if(type!=GLY){
      CB=p.getResidue(i)->getCB();
      cb_hdr=CB_HDR;
      cb_ir=CB_IR;
      /*CB-Backbone*/
      el_col_sprint(buf, static_cast<dmd_atom_t>(index), _N_,     min_NCb);
      out << buf << endl;
      el_col_sprint(buf, static_cast<dmd_atom_t>(index), _N_HB_,  min_NCb);
      out << buf << endl;
      el_col_sprint(buf, static_cast<dmd_atom_t>(index), _PRO_N_, min_NCb);
      out << buf << endl;
      
      el_col_sprint(buf, static_cast<dmd_atom_t>(index), _CA_,    cb_hdr+CA_HDR);
      out << buf << endl;
      
      el_col_sprint(buf, static_cast<dmd_atom_t>(index), _C_,     min_CCb);
      out << buf << endl;
      
      el_col_sprint(buf, static_cast<dmd_atom_t>(index), _O_,     min_OCb);
      out << buf << endl;
      el_col_sprint(buf, static_cast<dmd_atom_t>(index), _O_HB_,  min_OCb);
      out << buf << endl;
      
      if(type!=ALA){
	CG=p.getResidue(i)->getG();
	cg_hdr = G1_CONST[type][G1_HDR];
	cg_ir=G1_CONST[type][G1_IR];
	/*CG-Backbone*/
	el_col_sprint(buf, static_cast<dmd_atom_t>(index+1), _N_,     cg_hdr+N_HDR);
	out << buf << endl;
	el_col_sprint(buf, static_cast<dmd_atom_t>(index+1), _N_HB_,  cg_hdr+N_HDR);
	out << buf << endl;
	el_col_sprint(buf, static_cast<dmd_atom_t>(index+1), _PRO_N_, cg_hdr+N_HDR);
	out << buf << endl;
	
	el_col_sprint(buf, static_cast<dmd_atom_t>(index+1), _CA_,    cg_hdr+CA_HDR);
	out << buf << endl;
	
	el_col_sprint(buf, static_cast<dmd_atom_t>(index+1), _C_,     cg_hdr+C_HDR);
	out << buf << endl;
	
	el_col_sprint(buf, static_cast<dmd_atom_t>(index+1), _O_,     cg_hdr+O_HDR);
	out << buf << endl;
	el_col_sprint(buf, static_cast<dmd_atom_t>(index+1), _O_HB_,  cg_hdr+O_HDR);
	out << buf << endl;
	
	if(type==ILE||type==THR||type==VAL){
	  CG2=p.getResidue(i)->getG2();
	  if(type==ILE){
	    cg2_hdr=ILE_G2[G2_HDR]; cg2_ir =ILE_G2[G2_IR];
	  }
	  else if(type==THR){
	    cg2_hdr=THR_G2[G2_HDR]; cg2_ir =THR_G2[G2_IR];
	  }
	  else{
	    cg2_hdr=VAL_G2[G2_HDR]; cg2_ir =VAL_G2[G2_IR];
	  }
	  /*CG2-Backbone*/
	  el_col_sprint(buf, static_cast<dmd_atom_t>(index+2), _N_,     cg2_hdr+N_HDR);
	  out << buf << endl;
	  el_col_sprint(buf, static_cast<dmd_atom_t>(index+2), _N_HB_,  cg2_hdr+N_HDR);
	  out << buf << endl;
	  el_col_sprint(buf, static_cast<dmd_atom_t>(index+2), _PRO_N_, cg2_hdr+N_HDR);
	  out << buf << endl;
	  
	  el_col_sprint(buf, static_cast<dmd_atom_t>(index+2), _CA_,    cg2_hdr+CA_HDR);
	  out << buf << endl;
	  
	  el_col_sprint(buf, static_cast<dmd_atom_t>(index+2), _C_,     cg2_hdr+C_HDR);
	  out << buf << endl;
	  
	  el_col_sprint(buf, static_cast<dmd_atom_t>(index+2), _O_,     cg2_hdr+O_HDR);
	  out << buf << endl;
	  el_col_sprint(buf, static_cast<dmd_atom_t>(index+2), _O_HB_,  cg2_hdr+O_HDR);
	  out << buf << endl;
	}
      }
    }
    
  }
}

void printBOND_TYPE(ostream& out, pseudoPep& p, int* cbIndex){
  out << txtKeyWords[LINK_PAIRS] << endl;
  char buf[1024];
  /*N-Ca, Ca_N1*/
  double n_len;
  double n_dev = getTriDev_angle(CA_C, C_N, BOND_DEV, CA_C_N*rpi, n_len);
  double_bond_sprint(buf, N_CA, BOND_DEV, n_len, n_dev);
  out << _N_ << " " << _CA_ << " " << buf << endl;
  /*Ca-C, C-Ca1*/
  n_dev = getTriDev_angle(C_N, N_CA, BOND_DEV, C_N_CA*rpi, n_len);
  double_bond_sprint(buf, CA_C, BOND_DEV, n_len, n_dev);
  out << _C_ << " " << _CA_ << " " << buf << endl;
  /*C-N1, N-C*/
  n_dev = getTriDev_angle(N_CA, CA_C, BOND_DEV, N_CA_C*rpi, n_len);
  double_bond_sprint(buf, C_N, BOND_DEV, n_len, n_dev);
  out << _N_ << " " << _C_ << " " << buf << endl;
  /*C-C, used to alow certain phi region for PRO*/
  n_dev=(PRO_C_C1_MAX-PRO_C_C1_MIN)/(PRO_C_C1_MAX+PRO_C_C1_MIN);
  single_bond_sprint(buf, (PRO_C_C1_MAX+PRO_C_C1_MIN)/2.0, n_dev);
  //sprintf(buf, "%ld %ld %lf %lf", _C_, _C_, PRO_C_C1_MIN, PRO_C_C1_MAX);
  out << _C_ << " " << _C_ << " " << buf << endl;
  /*C-O*/
  single_bond_sprint(buf, C_O, BOND_DEV_O);
  out << _C_ << " " << _O_ << " " << buf << endl;
  /*CA-O*/
  n_dev = getTriDev_angle(CA_C, C_O, BOND_DEV_O, CA_C_O*rpi, n_len);
  single_bond_sprint(buf, n_len, n_dev);
  out << _O_ << " " << _CA_ << " " << buf << endl;
  /*O-N1*/
  n_dev = getTriDev_angle(C_O, C_N, BOND_DEV_O, (360.0-CA_C_N-CA_C_O)*rpi, n_len);
  single_bond_sprint(buf, n_len, n_dev);
  out << _O_ << " " << _N_ << " " << buf << endl;
  /*CA-CA*/
  double c_ca = sqrt(C_N*C_N + N_CA*N_CA - 2.0*C_N*N_CA*cos(C_N_CA*rpi));
  n_dev = getTriDev_angle(CA_C, c_ca, BOND_DEV, CA_C_CA*rpi, n_len);
  single_bond_sprint(buf, n_len, n_dev);
  out << _CA_ << " " << _CA_ << " " << buf << endl;
  /*sidechain-backbone bonding*/
  for(int i=0; i<p.getLength(); i++){
    aa_t type=static_cast<aa_t>(p.getResidue(i)->getID());
    int index=cbIndex[i]+_O_HB_+1;
    if(type!=GLY){
      /*CA-CB*/
      single_bond_sprint(buf, CA_CB, BOND_DEV);
      out << _CA_ << " " << index << " " << buf << endl;
      /*N-CB*/
      n_dev = getTriDev_angle(N_CA, CA_CB, BOND_DEV, N_CA_CB*rpi, n_len);
      if(type==PRO)
	n_dev = getTriDev_angle(N_CA, CA_CB, BOND_DEV, PRO_N_CA_CB*rpi, n_len);
      single_bond_sprint(buf, n_len, n_dev);
      out << _N_ << " " << index << " " << buf << endl;
      /*C-CB*/
      n_dev = getTriDev_angle(CA_C, CA_CB, BOND_DEV, C_CA_CB*rpi, n_len);
      single_bond_sprint(buf, n_len, n_dev);
      if(type==PRO){
	n_dev = getTriDev_angle(CA_C, CA_CB, BOND_DEV, PRO_C_CA_CB*rpi, n_len);
	double_bond_sprint(buf, n_len, n_dev, PRO_C_CB, PRO_C_CB_D/PRO_C_CB);
      }
      out << _C_ << " " << index << " " << buf << endl;
      
      if(type!=ALA){
	const double* cg_para = G1_CONST[type];
	if(cg_para[G1_CB_D]==0){
	  single_bond_sprint(buf, cg_para[G1_CB], BOND_DEV_G);
	}
	else{
	  single_bond_sprint_len(buf, cg_para[G1_CB], cg_para[G1_CB_D]);
	}
	out << index << " " << index+1 << " " << buf << endl;
	
	if(cg_para[G1_CA_D]==0){
	  n_dev = getTriDev_length(cg_para[G1_CB], CA_CB, BOND_DEV_G, cg_para[G1_CA]);
	  single_bond_sprint(buf, cg_para[G1_CA], n_dev);
	}
	else{
	  single_bond_sprint_len(buf, cg_para[G1_CA], cg_para[G1_CA_D]);
	}
	out << _CA_ << " " << index+1 << " " << buf << endl;
	
	/*N-CG, C-CG*/
	const double* dihedral = Dihedral_C_CG[type];
	if(dihedral[0]!=INF){
	  //c-cg
	  sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf %lf %lf",
		  _C_, index+1, 
		  Dihedral_MIN, 
		  dihedral[0], Dihedral_E,
		  dihedral[1], -Dihedral_E, 
		  dihedral[2], Dihedral_E, 
		  Dihedral_MAX);
	  out << buf << endl;
	  //n-cg
	  sprintf(buf, "%ld %ld %lf %lf",
		  _N_, index+1, Dihedral_MIN, Dihedral_MAX);
	  out << buf << endl;
	}
	if(type==PRO){//different constraint for PROLINE
	  single_bond_sprint_len(buf, PRO_N_CG, PRO_N_CG_D);
	  out << _N_ << " " << index+1 << " " << buf << endl;
	  
	  /*this constraint is actually used to
	    avoid the small distance of c-PRO_CB*/
	  single_bond_sprint_len(buf, PRO_C_CG, PRO_C_CG_D);
	  out << _C_ << " " << index+1 << " " << buf << endl;
	}
	if(type==ILE){
	  /*ILE, G2*/
	  //CG2-CB
	  single_bond_sprint(buf, ILE_G2[G2_CB], BOND_DEV_G);
	  out << index << " " << index+2 << " " << buf << endl;
	  //CG2-CA
	  n_dev = getTriDev_length(CA_CB, ILE_G2[G2_CB], BOND_DEV_G, ILE_G2[G2_CA]);
	  single_bond_sprint(buf, ILE_G2[G2_CA], n_dev);
	  out << _CA_ << " " << index+2 << " " << buf << endl;
	  //CG2-CG1
	  n_dev = getTriDev_length(ILE_G2[G2_CB], G1_CONST[ILE][G1_CB], 
				   BOND_DEV_G, ILE_G2[G1_G2]);
	  single_bond_sprint(buf, ILE_G2[G1_G2], n_dev);
	  out << index+1 << " " << index+2 << " " << buf << endl;
	  
	  //N-CG2
	  sprintf(buf, "%ld %ld %lf %lf",
		  _N_, index+2, Dihedral_MIN, Dihedral_MAX);
	  out << buf << endl;
	  //C-CG2
	  const double* dihedral = Dihedral_C_CG[ILE];
	  sprintf(buf, "%ld %ld %lf %lf",
		  _C_, index+2, Dihedral_MIN, Dihedral_MAX);
	  out << buf << endl;
	  
	}
	else if(type==THR){
	  /*THR, G2*/
	  //CG2-CB
	  single_bond_sprint(buf, THR_G2[G2_CB], BOND_DEV_G);
	  out << index << " " << index+2 << " " << buf << endl;
	  //CG2-CA
	  n_dev = getTriDev_length(CA_CB, THR_G2[G2_CB], BOND_DEV_G, THR_G2[G2_CA]);
	  single_bond_sprint(buf, THR_G2[G2_CA], n_dev);
	  out << _CA_ << " " << index+2 << " " << buf << endl;
	  //CG2-CG1
	  n_dev = getTriDev_length(THR_G2[G2_CB], G1_CONST[THR][G1_CB], 
				   BOND_DEV_G, THR_G2[G1_G2]);
	  single_bond_sprint(buf, THR_G2[G1_G2], n_dev);
	  out << index+1 << " " << index+2 << " " << buf << endl;
	  
	  //N-CG2
	  sprintf(buf, "%ld %ld %lf %lf",
		  _N_, index+2, Dihedral_MIN, Dihedral_MAX);
	  out << buf << endl;
	  //C-CG2
	  sprintf(buf, "%ld %ld %lf %lf",
		  _C_, index+2, Dihedral_MIN, Dihedral_MAX);
	  out << buf << endl;
	}
	else if(type==VAL){
	  /*VAL, G2*/
	  //CG2-CB
	  single_bond_sprint(buf, VAL_G2[G2_CB], BOND_DEV_G);
	  out << index << " " << index+2 << " " << buf << endl;
	  //CG2-CA
	  n_dev = getTriDev_length(CA_CB, VAL_G2[G2_CB], BOND_DEV_G, VAL_G2[G2_CA]);
	  single_bond_sprint(buf, VAL_G2[G2_CA], n_dev);
	  out << _CA_ << " " << index+2 << " " << buf << endl;
	  //CG2-CG1
	  n_dev = getTriDev_length(VAL_G2[G2_CB], G1_CONST[VAL][G1_CB], 
				   BOND_DEV_G, VAL_G2[G1_G2]);
	  single_bond_sprint(buf, VAL_G2[G1_G2], n_dev);
	  out << index+1 << " " << index+2 << " " << buf << endl;
	  
	  //N-CG2
	  sprintf(buf, "%ld %ld %lf %lf",
		  _N_, index+2, Dihedral_MIN, Dihedral_MAX);
	  out << buf << endl;
	  //C-CG2
	  sprintf(buf, "%ld %ld %lf %lf",
		  _C_, index+2, Dihedral_MIN, Dihedral_MAX);
	  out << buf << endl;
	}
	
      }
    }
  }
  /*sidechain-sidechain*/
  for(int i=0; i<p.getLength(); i++){
    aa_t type=static_cast<aa_t>(p.getResidue(i)->getID());
    int index=cbIndex[i]+_O_HB_+1;
    atom* CB=NULL;
    atom* CG=NULL;
    atom* CG2=NULL;
    double cb_hdr, cb_ir;
    double cg_hdr, cg_ir;
    double cg2_hdr, cg2_ir;
    if(type!=GLY){
      CB=p.getResidue(i)->getCB();
      cb_hdr=CB_HDR;
      cb_ir=CB_IR;
      
      if(type!=ALA){
	CG=p.getResidue(i)->getG();
	cg_hdr = G1_CONST[type][G1_HDR];
	cg_ir=G1_CONST[type][G1_IR];

	if(type==ILE||type==THR||type==VAL){
	  CG2=p.getResidue(i)->getG2();
	  if(type==ILE){
	    cg2_hdr=ILE_G2[G2_HDR];
	    cg2_ir =ILE_G2[G2_IR];
	  }
	  else if(type==THR){
	    cg2_hdr=THR_G2[G2_HDR];
	    cg2_ir =THR_G2[G2_IR];
	  }
	  else{
	    cg2_hdr=VAL_G2[G2_HDR];
	    cg2_ir =VAL_G2[G2_IR];
	  }
	}
      }
    }
    
    for(int j=i+1; j<p.getLength(); j++){
      aa_t typep=static_cast<aa_t>(p.getResidue(j)->getID());
      int indexp=cbIndex[j]+_O_HB_+1;
      atom* CBp=NULL;
      atom* CGp=NULL;
      atom* CG2p=NULL;
      double hdr, ir;
      if(typep!=GLY){
	CBp=p.getResidue(j)->getCB();
	hdr=CB_HDR;
	ir =CB_IR;
	if(CB){
	  if(side_sprint(buf, CB, index, CBp, indexp,  hdr+cb_hdr, ir+cb_ir))
	    out << buf << endl;
	}
	if(CG){
	  if(side_sprint(buf, CG, index+1, CBp, indexp,  hdr+cg_hdr, ir+cg_ir))
	    out << buf << endl;
	}
	if(CG2){
	  if(side_sprint(buf, CG2,index+2, CBp, indexp,  hdr+cg2_hdr, ir+cg2_ir))
	    out << buf << endl;
	}
	if(typep!=ALA){
	  CGp=p.getResidue(j)->getG();
	  hdr=G1_CONST[typep][G1_HDR];
	  ir=G1_CONST[typep][G1_IR];
	  if(CB){
	    if(side_sprint(buf, CB, index,  CGp, indexp+1,  hdr+cb_hdr, ir+cb_ir))
	      out << buf << endl;
	  }
	  if(CG){
	    if((type==PHE || type==TRP || type==TYR) && (typep==PHE || typep==TRP || typep==TYR)){
	      double thdr=0;
	      if(type==PHE)     thdr+=PHE_ARO_HDR;
	      else if(type==TRP)thdr+=TRP_ARO_HDR;
	      else              thdr+=TYR_ARO_HDR;
	      
	      if(typep==PHE)     thdr+=PHE_ARO_HDR;
	      else if(typep==TRP)thdr+=TRP_ARO_HDR;
	      else               thdr+=TYR_ARO_HDR;
	      if(side_sprint(buf, CG, index+1, CGp, indexp+1, thdr, ir+cg_ir))
		out << buf << endl;
	    }
	    else{
	      if(side_sprint(buf, CG, index+1, CGp, indexp+1,  hdr+cg_hdr, ir+cg_ir))
		out << buf << endl;
	    }
	  }
	  if(CG2){
	    if(side_sprint(buf, CG2,index+2, CGp,indexp+1,  hdr+cg2_hdr, ir+cg2_ir))
	      out << buf << endl;
	  }
	  
	  
	  if(typep==ILE||typep==THR||typep==VAL){
	    CG2p=p.getResidue(j)->getG2();
	    if(typep==ILE){
	      hdr=ILE_G2[G2_HDR];
	      ir =ILE_G2[G2_IR];
	    }
	    else if(typep==THR){
	      hdr=THR_G2[G2_HDR];
	      ir =THR_G2[G2_IR];
	    }
	    else{
	      hdr=VAL_G2[G2_HDR];
	      ir =VAL_G2[G2_IR];
	    }
	    if(CB){
	      if(side_sprint(buf, CB, index,   CG2p, indexp+2, hdr+cb_hdr, ir+cb_ir))
		out << buf << endl;
	    }
	    if(CG){
	      if(side_sprint(buf, CG, index+1, CG2p, indexp+2, hdr+cg_hdr, ir+cg_ir))
		out << buf << endl;
	    }
	    if(CG2){
	      if(side_sprint(buf, CG2,index+2, CG2p, indexp+2, hdr+cg2_hdr, ir+cg2_ir))
		out << buf << endl;
	    }
	  }
	}
	
      }
    }

  }
  
  /*mainchain hydrogen bond related*/
  /*NP-OP*/
  sprintf(buf,"%ld %ld %lf %lf %lf", _N_HB_, _O_HB_, 
	  HB_N_O[0], HB_N_O[1]+eps, -fabs(EHB_M_M));
  out << buf << endl;
  /*NP-C*/
  sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf 0.00000", 
  	  _N_HB_, _C_,
  	  HB_N_C[0], 
  	  HB_N_C[1], fabs(EHB_M_M)/4.0, 
	  HB_N_C[2], fabs(EHB_M_M)/2.0,
  	  HB_N_C[3]+eps);
  out << buf << endl;
  //sprintf(buf, "%ld %ld %lf %lf", 
  //_N_HB_, _C_,
  //HB_N_C[1],
  //HB_N_C[2]), 
  //out << buf << endl;
  /*OP-C*/
  sprintf(buf, "%ld %ld %lf %lf %lf %lf 0.00000", 
  	  _O_HB_, _C_,
  	  HB_O_C[0], 
  	  HB_O_C[1], fabs(EHB_M_M)*3.0/4.0, 
  	  HB_O_C[2]+eps); 
  out << buf << endl;
  //sprintf(buf, "%ld %ld %lf %lf", 
  //	  _O_HB_, _C_,
  //	  HB_O_C[1],
  //	  HB_O_C[2]);
  //out << buf << endl;
  /*OP-CA*/
  sprintf(buf, "%ld %ld %lf %lf %lf %lf 0.0000", 
  	  _O_HB_, _CA_,
  	  HB_O_CA[0], 
  	  HB_O_CA[1], fabs(EHB_M_M)*3.0/4.0, 
  	  HB_O_CA[2]+eps);
  out << buf << endl;
  //sprintf(buf, "%ld %ld %lf %lf", 
  //	  _O_HB_, _CA_,
  //	  HB_O_CA[1],
  //	  HB_O_CA[2]);
  //out << buf << endl;

}

void printREACT(ostream& out){
  out << txtKeyWords[REACT] << endl;
  char buf[1024];
  /*mainchain hydrogen bonding*/
  sprintf(buf, "%ld %ld %ld %ld 1",
  	  _N_, _O_, _N_HB_, _O_HB_);
  out << buf << endl;
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _N_, _C_, _N_HB_, _C_);
  out << buf << endl;
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _O_, _C_, _O_HB_, _C_);
  out << buf << endl;
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _O_, _CA_, _O_HB_, _CA_);
  out << buf << endl;
}

void printATOM_LIST(ostream& out, pseudoPep& p, randomGenerator& r, int* cbIndex){
  out << txtKeyWords[LIST_ATOMS] << endl;
  char buf[1024];
  vec v;
  int index=1;
  atom* CB;
  atom* CG;
  atom* CG2;
  dmd_atom_t type;
  for(int i=0; i<p.getLength(); i++){
    pseudoAA* theAA = p.getResidue(i);
    /*N*/
    v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
    if(theAA->getID()==PRO)
      list_atom_sprint(buf, index++, _PRO_N_, *theAA->getN()->getR(), v);
    else{
      if(i==0)/*first N atom will be hydrogen-active!*/
	list_atom_sprint(buf, index++, _N_HB_, *theAA->getN()->getR(), v);
      else
	list_atom_sprint(buf, index++, _N_, *theAA->getN()->getR(), v);
    }
    out << buf << endl;
      
    /*C*/
    v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
    list_atom_sprint(buf, index++, _C_, *theAA->getC()->getR(), v);
    out << buf << endl;
    /*O*/
    v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
    if(i<p.getLength()-1){
      list_atom_sprint(buf, index++, _O_, *theAA->getO()->getR(), v);
    }
    else{/*last O atom will be hydrogen-active!*/
      list_atom_sprint(buf, index++, _O_HB_, *theAA->getO()->getR(), v);
    }
    out << buf << endl;
    /*CA*/
    v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
    list_atom_sprint(buf, index++, _CA_, *theAA->getCA()->getR(), v);
    out << buf << endl;
    /*CB*/
    if(CB = theAA->getCB()){
      type = static_cast<dmd_atom_t>(_O_HB_+1+cbIndex[i]);
      v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
      list_atom_sprint(buf, index++, type, *CB->getR(), v);
      out << buf << endl;
      /*CG*/
      if(CG = theAA->getG()){
	type = static_cast<dmd_atom_t>(_O_HB_+2+cbIndex[i]);
	double mass = G1_CONST[theAA->getID()][G1_M];
	v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss())/sqrt(mass);
	list_atom_sprint(buf, index++, type, *CG->getR(), v);
	out << buf << endl;
	if(CG2 = theAA->getG2()){
	  if(theAA->getID()==ILE) 
	    type = static_cast<dmd_atom_t>(_O_HB_+3+cbIndex[i]);
	  else if(theAA->getID()==THR) 
	    type = static_cast<dmd_atom_t>(_O_HB_+3+cbIndex[i]);
	  else if(theAA->getID()==VAL) 
	    type = static_cast<dmd_atom_t>(_O_HB_+3+cbIndex[i]);
	  else{
	    cerr << "error in type" << endl;
	    exit(1);
	  }
	  v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
	  list_atom_sprint(buf, index++, type, *CG2->getR(), v);
	  out << buf << endl;
	}
      }
    }
  }
}

void printBOND_LIST(ostream& out, pseudoPep& p, int* cbIndex){
  out << txtKeyWords[LIST_BONDS] << endl;
  char buf[1024];
  for(int i=0; i<p.getLength(); i++){
    /*intra amino acids*/
    pseudoAA* theAA = p.getResidue(i);
    /*N-CA*/
    out << theAA->getN()->getIndex() << " " << theAA->getCA()->getIndex() << endl;
    /*CA-C*/
    out << theAA->getCA()->getIndex() << " " << theAA->getC()->getIndex() << endl;
    /*N-C*/
    out << theAA->getN()->getIndex() << " " << theAA->getC()->getIndex() << endl;
    /*C-O*/
    out << theAA->getC()->getIndex() << " " << theAA->getO()->getIndex() << endl;    
    /*CA-O*/
    out << theAA->getCA()->getIndex() << " " << theAA->getO()->getIndex() << endl;        
    if(theAA->getCB()){/*not GLY*/
      /*CA-CB*/
      out << theAA->getCA()->getIndex() << " " << theAA->getCB()->getIndex() << endl;
      /*N-CB*/
      out << theAA->getN()->getIndex() << " " << theAA->getCB()->getIndex() << endl;
      /*C-CB*/
      out << theAA->getC()->getIndex() << " " << theAA->getCB()->getIndex() << endl;
      if(theAA->getG()){/*not ALA*/
	/*CB-CG*/
	out << theAA->getCB()->getIndex() << " " << theAA->getG()->getIndex() << endl;
	/*CA-CG*/
	out << theAA->getCA()->getIndex() << " " << theAA->getG()->getIndex() << endl;
	/*N-CG, C-CG*/
	int id = theAA->getID();
	if(Dihedral_C_CG[id][0]!=INF){
	  out << theAA->getC()->getIndex() << " " << theAA->getG()->getIndex() << endl;
	  out << theAA->getN()->getIndex() << " " << theAA->getG()->getIndex() << endl;
	}	
	if(theAA->getG2()){/*beta-branched*/
	  /*CG2-CB*/
	  out << theAA->getCB()->getIndex() << " " << theAA->getG2()->getIndex() << endl;
	  /*CG2-CA*/
	  out << theAA->getCA()->getIndex() << " " << theAA->getG2()->getIndex() << endl;
	  /*CG2-CG*/
	  out << theAA->getG()->getIndex() << " " << theAA->getG2()->getIndex() << endl;
	  /*C-CG2*/
	  out << theAA->getC()->getIndex() << " " << theAA->getG2()->getIndex() << endl;
	  /*N-CG2*/
	  out << theAA->getN()->getIndex() << " " << theAA->getG2()->getIndex() << endl;
	}
      }
    }
    /*extra constraints for PRO*/
    if(theAA->getID()==PRO){
      out << theAA->getN()->getIndex() << " " << theAA->getG()->getIndex() << endl;
    }

    /*inter amino acids*/
    if(i>0){
      pseudoAA* preAA = p.getResidue(i-1);
      /*Ca0 -- N1*/
      out << preAA->getCA()->getIndex() << " " << theAA->getN()->getIndex() << endl;
      /*C0 -- Ca1*/
      out << preAA->getC()->getIndex() << " " << theAA->getCA()->getIndex() << endl;
      /*C0 -- N1*/
      out << preAA->getC()->getIndex() << " " << theAA->getN()->getIndex() << endl;
      /*O0--N1*/
      out << preAA->getO()->getIndex() << " " << theAA->getN()->getIndex() << endl;
      /*Ca0 -- Ca1*/
      out << preAA->getCA()->getIndex() << " " << theAA->getCA()->getIndex() << endl;
      if(theAA->getID()==PRO){
	out << preAA->getC()->getIndex() << " " << theAA->getCB()->getIndex() << endl;
	out << preAA->getC()->getIndex() << " " << theAA->getG()->getIndex() << endl;
	out << preAA->getC()->getIndex() << " " << theAA->getC()->getIndex() << endl;
      }
    }
  }
  
  for(int i=0; i<p.getLength(); i++){
    aa_t type=static_cast<aa_t>(p.getResidue(i)->getID());
    int index=cbIndex[i]+_O_HB_+1;
    atom* CB=NULL;
    atom* CG=NULL;
    atom* CG2=NULL;
    double cb_hdr, cb_ir;
    double cg_hdr, cg_ir;
    double cg2_hdr, cg2_ir;
    if(type!=GLY){
      CB=p.getResidue(i)->getCB();
      cb_hdr=CB_HDR;
      cb_ir=CB_IR;
      
      if(type!=ALA){
	CG=p.getResidue(i)->getG();
	cg_hdr = G1_CONST[type][G1_HDR];
	cg_ir=G1_CONST[type][G1_IR];
	
	if(type==ILE||type==THR||type==VAL){
	  CG2=p.getResidue(i)->getG2();
	  if(type==ILE){
	    cg2_hdr=ILE_G2[G2_HDR];
	    cg2_ir =ILE_G2[G2_IR];
	  }
	  else if(type==THR){
	    cg2_hdr=THR_G2[G2_HDR];
	    cg2_ir =THR_G2[G2_IR];
	  }
	  else{
	    cg2_hdr=VAL_G2[G2_HDR];
	    cg2_ir =VAL_G2[G2_IR];
	  }
	}
      }
    }
    for(int j=i+1; j<p.getLength(); j++){
      aa_t typep=static_cast<aa_t>(p.getResidue(j)->getID());
      int indexp=cbIndex[j]+_O_HB_+1;
      atom* CBp=NULL;
      atom* CGp=NULL;
      atom* CG2p=NULL;
      double hdr, ir;
      if(typep!=GLY){
	CBp=p.getResidue(j)->getCB();
	hdr=CB_HDR;
	ir =CB_IR;
	if(CB){
	  if(CB->getDist(*CBp)<(ir+cb_ir)) 
	    out << CB->getIndex() << " " << CBp->getIndex() << endl;
	}
	if(CG){
	  if(CG->getDist(*CBp)<(ir+cg_ir))
	    out << CG->getIndex() << " " << CBp->getIndex() << endl;
	}
	if(CG2){
	  if(CG2->getDist(*CBp)<(ir+cg2_ir))
	    out << CG2->getIndex() << " " << CBp->getIndex() << endl;
	}
	if(typep!=ALA){
	  CGp=p.getResidue(j)->getG();
	  hdr=G1_CONST[typep][G1_HDR];
	  ir=G1_CONST[typep][G1_IR];
	  if(CB){
	    if(CB->getDist(*CGp)<(ir+cb_ir)) 
	      out << CB->getIndex() << " " << CGp->getIndex() << endl;
	  }
	  if(CG){
	    if(CG->getDist(*CGp)<(ir+cg_ir))
	      out << CG->getIndex() << " " << CGp->getIndex() << endl;
	  }
	  if(CG2){
	    if(CG2->getDist(*CGp)<(ir+cg2_ir))
	      out << CG2->getIndex() << " " << CGp->getIndex() << endl;
	  }
	  
	  if(typep==ILE||typep==THR||typep==VAL){
	    CG2p=p.getResidue(j)->getG2();
	    if(typep==ILE){
	      hdr=ILE_G2[G2_HDR];
	      ir =ILE_G2[G2_IR];
	    }
	    else if(typep==THR){
	      hdr=THR_G2[G2_HDR];
	      ir =THR_G2[G2_IR];
	    }
	    else{
	      hdr=VAL_G2[G2_HDR];
	      ir =VAL_G2[G2_IR];
	    }
	    if(CB){
	      if(CB->getDist(*CG2p)<(ir+cb_ir)) 
		out << CB->getIndex() << " " << CG2p->getIndex() << endl;
	    }
	    if(CG){
	      if(CG->getDist(*CG2p)<(ir+cg_ir))
		out << CG->getIndex() << " " << CG2p->getIndex() << endl;
	    }
	    if(CG2){
	      if(CG2->getDist(*CG2p)<(ir+cg2_ir))
		out << CG2->getIndex() << " " << CG2p->getIndex() << endl;
	    }
	  }
	}
      }
    }
  }
}

void printPERM_BOND_LIST(ostream& out, pseudoPep& p, int* cbIndex){
  out << txtKeyWords[LIST_PERM_BONDS] << endl;
  char buf[1024];
  for(int i=0; i<p.getLength(); i++){
    /*intra amino acids*/
    pseudoAA* theAA = p.getResidue(i);
    /*N-CA*/
    out << theAA->getN()->getIndex() << " " << theAA->getCA()->getIndex() << " "
	<< _N_ << " " << _CA_ << endl;
    /*CA-C*/
    out << theAA->getCA()->getIndex() << " " << theAA->getC()->getIndex() << " "
	<< _CA_ << " " << _C_ << endl;
    /*N-C*/
    out << theAA->getN()->getIndex() << " " << theAA->getC()->getIndex() << " "
	<< _N_ << " " << _C_ << endl;    
    /*C-O*/
    out << theAA->getC()->getIndex() << " " << theAA->getO()->getIndex() << " " 
	<< _C_ << " " << _O_ << endl;    
    /*CA-O*/
    out << theAA->getCA()->getIndex() << " " << theAA->getO()->getIndex() << " " 
	<< _CA_ << " " << _O_ << endl;        
    if(theAA->getCB()){/*not GLY*/
      dmd_atom_t  cbt = static_cast<dmd_atom_t>(cbIndex[i]+1+_O_HB_);
      if(theAA->getID()==PRO) cbt = _PRO_CB_;
      /*CA-CB*/
      out << theAA->getCA()->getIndex() << " " << theAA->getCB()->getIndex() << " "
	  << _CA_ << " " << cbt << endl;
      /*N-CB*/
      out << theAA->getN()->getIndex() << " " << theAA->getCB()->getIndex() << " "
	  << _N_ << " " << cbt << endl;
      /*C-CB*/
      out << theAA->getC()->getIndex() << " " << theAA->getCB()->getIndex() << " "
	  << _C_ << " " << cbt  << endl;
      if(theAA->getG()){/*not ALA*/
	dmd_atom_t gt = static_cast<dmd_atom_t>(cbIndex[i]+2+_O_HB_);
	/*CB-CG*/
	out << theAA->getCB()->getIndex() << " " << theAA->getG()->getIndex() << " "
	    << cbt << " " << gt << endl;
	/*CA-CG*/
	out << theAA->getCA()->getIndex() << " " << theAA->getG()->getIndex() << " "
	    << _CA_ << " " << gt << endl;
	if(Dihedral_C_CG[theAA->getID()][0]!=INF){
	  /*N-CG*/
	  out << theAA->getN()->getIndex() << " " << theAA->getG()->getIndex() << " "
	      << _N_ << " " << gt << endl;
	  /*C-CG*/
	  out << theAA->getC()->getIndex() << " " << theAA->getG()->getIndex() << " "
	      << _C_ << " " << gt << endl;
	}
	if(theAA->getG2()){/*beta-branched*/
	  dmd_atom_t g2t=static_cast<dmd_atom_t>(cbIndex[i]+3+_O_HB_);
	  /*CB-CG2*/
	  out << theAA->getCB()->getIndex() << " " << theAA->getG2()->getIndex() << " "
	      << cbt << " " << g2t << endl;
	  /*CA-CG2*/
	  out << theAA->getCA()->getIndex() << " " << theAA->getG2()->getIndex() << " "
	      << _CA_ << " " << g2t << endl;
	  /*CG-CG2*/
	  out << theAA->getG()->getIndex() << " " << theAA->getG2()->getIndex() << " "
	      << gt << " " << g2t << endl;
	  /*N-CG2*/
	  out << theAA->getN()->getIndex() << " " << theAA->getG2()->getIndex() << " "
	      << _N_ << " " << g2t << endl;
	  /*C-CG2*/
	  out << theAA->getC()->getIndex() << " " << theAA->getG2()->getIndex() << " "
	      << _C_ << " " << g2t << endl;
	}
	/*extra constraints for PRO: N-CG*/
	if(theAA->getID()==PRO){
	  out << theAA->getN()->getIndex() << " " << theAA->getG()->getIndex() << " "
	      << _N_ << " " << gt << endl;
	}
      }
    }

    /*inter amino acids*/
    if(i>0){
      pseudoAA* preAA = p.getResidue(i-1);
      /*Ca0 -- N1*/
      out << preAA->getCA()->getIndex() << " " << theAA->getN()->getIndex() << " "
	  << _CA_ << " " << _N_ << endl;
      /*C0 -- Ca1*/
      out << preAA->getC()->getIndex() << " " << theAA->getCA()->getIndex() << " "
	  << _C_ << " " << _CA_ << endl;
      /*C0 -- N1*/
      out << preAA->getC()->getIndex() << " " << theAA->getN()->getIndex() << " "
	  << _C_ << " " << _N_ << endl;
      /*O0--N1*/
      out << preAA->getO()->getIndex() << " " << theAA->getN()->getIndex() << " " 
	  << _O_ << " " << _N_ << endl;
      /*Ca0 -- Ca1*/
      out << preAA->getCA()->getIndex() << " " << theAA->getCA()->getIndex() << " "
	  << _CA_ << " " << _CA_ << endl;
      if(theAA->getID()==PRO){
	dmd_atom_t cb_t = static_cast<dmd_atom_t>(cbIndex[i]+1+_O_HB_);
	dmd_atom_t cg_t = static_cast<dmd_atom_t>(cbIndex[i]+2+_O_HB_);
	out << preAA->getC()->getIndex() << " " << theAA->getCB()->getIndex() << " "
	    << _C_ << " " << cb_t << endl;
	out << preAA->getC()->getIndex() << " " << theAA->getG()->getIndex() << " "
	    << _C_ << " " << cg_t << endl;
	out << preAA->getC()->getIndex() << " " << theAA->getC()->getIndex() << " "
	    << _C_ << " " << _C_ << endl;
      }
    }
  }
  
  for(int i=0; i<p.getLength(); i++){
    aa_t type=static_cast<aa_t>(p.getResidue(i)->getID());
    int index=cbIndex[i]+_O_HB_+1;
    atom* CB=NULL;
    atom* CG=NULL;
    atom* CG2=NULL;
    double cb_hdr, cb_ir;
    double cg_hdr, cg_ir;
    double cg2_hdr, cg2_ir;
    if(type!=GLY){
      CB=p.getResidue(i)->getCB();
      cb_hdr=CB_HDR;
      cb_ir=CB_IR;
      
      if(type!=ALA){
	CG=p.getResidue(i)->getG();
	cg_hdr = G1_CONST[type][G1_HDR];
	cg_ir=G1_CONST[type][G1_IR];
	
	if(type==ILE||type==THR||type==VAL){
	  CG2=p.getResidue(i)->getG2();
	  if(type==ILE){
	    cg2_hdr=ILE_G2[G2_HDR];
	    cg2_ir =ILE_G2[G2_IR];
	  }
	  else if(type==THR){
	    cg2_hdr=THR_G2[G2_HDR];
	    cg2_ir =THR_G2[G2_IR];
	  }
	  else{
	    cg2_hdr=VAL_G2[G2_HDR];
	    cg2_ir =VAL_G2[G2_IR];
	  }
	}
      }
    }
    for(int j=i+1; j<p.getLength(); j++){
      aa_t typep=static_cast<aa_t>(p.getResidue(j)->getID());
      int indexp=cbIndex[j]+_O_HB_+1;
      atom* CBp=NULL;
      atom* CGp=NULL;
      atom* CG2p=NULL;
      double hdr, ir;
      if(typep!=GLY){
	CBp=p.getResidue(j)->getCB();
	hdr=CB_HDR;
	ir =CB_IR;
	if(CB){
	  if(CB->getDist(*CBp)<(ir+cb_ir)) {
	    out << CB->getIndex() << " " << CBp->getIndex() << " "
		<< index << " " << indexp << endl;
	    //cout << "CB-CB" << endl;
	  }
	}
	if(CG){
	  if(CG->getDist(*CBp)<(ir+cg_ir)){
	    out << CG->getIndex() << " " << CBp->getIndex() << " "
		<< index+1 << " " << indexp << endl;
	    //cout << "CG-CB" << endl;
	  }
	}
	if(CG2){
	  if(CG2->getDist(*CBp)<(ir+cg2_ir)){
	    out << CG2->getIndex() << " " << CBp->getIndex() << " " 
		<< index+2 << " " << indexp << endl;
	    //cout << "CG2-CB" << endl;
	  }
	}
	if(typep!=ALA){
	  CGp=p.getResidue(j)->getG();
	  hdr=G1_CONST[typep][G1_HDR];
	  ir=G1_CONST[typep][G1_IR];
	  if(CB){
	    if(CB->getDist(*CGp)<(ir+cb_ir)) {
	      out << CB->getIndex() << " " << CGp->getIndex() << " " 
		  << index << " " << indexp+1 << endl;
	      //cout << "CB-CG" << endl;
	    }
	  }
	  if(CG){
	    if(CG->getDist(*CGp)<(ir+cg_ir))
	      out << CG->getIndex() << " " << CGp->getIndex() << " " 
		  << index+1 << " " << indexp+1 << endl;
	  }
	  if(CG2){
	    if(CG2->getDist(*CGp)<(ir+cg2_ir))
	      out << CG2->getIndex() << " " << CGp->getIndex() << " "
		  << index+2 << " " << indexp+1 << endl;
	  }
	  
	  if(typep==ILE||typep==THR||typep==VAL){
	    CG2p=p.getResidue(j)->getG2();
	    if(typep==ILE){
	      hdr=ILE_G2[G2_HDR];
	      ir =ILE_G2[G2_IR];
	    }
	    else if(typep==THR){
	      hdr=THR_G2[G2_HDR];
	      ir =THR_G2[G2_IR];
	    }
	    else{
	      hdr=VAL_G2[G2_HDR];
	      ir =VAL_G2[G2_IR];
	    }
	    if(CB){
	      if(CB->getDist(*CG2p)<(ir+cb_ir)){
		out << CB->getIndex() << " " << CG2p->getIndex() << " " 
		    << index << " " << indexp+2 <<  endl;
		//cout << "CB-CG2" << endl;
	      }
	    }
	    if(CG){
	      if(CG->getDist(*CG2p)<(ir+cg_ir))
		out << CG->getIndex() << " " << CG2p->getIndex() << " "
		    << index+1 << " " << indexp+2 << endl;
	    }
	    if(CG2){
	      if(CG2->getDist(*CG2p)<(ir+cg2_ir))
		out << CG2->getIndex() << " " << CG2p->getIndex() << " " 
		    << index+2 << " " << indexp+2 << endl;
	    }
	  }
	}
      }
    }
  }
}

void printHBA_LIST(ostream& out, pseudoPep& p){
  out << txtKeyWords[HB_LIST] << endl;
  char buf[1024];
  for(int i=0; i<p.getLength(); i++){
    /*N*/
    if(i>0){
      sprintf(buf, "%ld %ld %ld %ld",
	      p.getResidue(i)->getN()->getIndex(), i, 
	      p.getResidue(i-1)->getC()->getIndex(), 
	      p.getResidue(i)->getCA()->getIndex());
      out << buf << endl;
    }
    /*O*/
    if(i<p.getLength()-1){
      sprintf(buf, "%ld %ld %ld",
	      p.getResidue(i)->getO()->getIndex(), i, 
	      p.getResidue(i)->getC()->getIndex());
      out << buf << endl;
    }
  }
}

/*Backbone can form hydrogen bonds, while the sidechain keep the
  contactmap. Therefore, while been feed in the PDB file, the program will
  produce a dmd-compatible txt format to use dmd to relax the protein
  conformation. */

int main(int argc, char* argv[]){
  if(argc<3){
    cout << "usage: preRelax.linux PDB boxSize" << endl;
    exit(1);
  }
  
  pseudoPep p(argv[1], static_cast<input_t>(2));
  pseudoAA* lastAA=p.getResidue(p.getLength()-1);
  lastAA->getO()->getR()->Rotate(*lastAA->getCA()->getR(),*lastAA->getC()->getR(), -PI/3.0);
  p.makeIndex();
  
  randomGenerator ran(_RAN2_, -100);

  double size = atof(argv[2]);
  vec cm = getCM(p);
  vec shift = vec(size/2.0, size/2.0, size/2.0) - cm;
  p.shift(shift);
  
  int* cbSideIndex =  new int[p.getLength()];
  int dummy=0;
  for(int i=0; i<p.getLength();i++){
    cbSideIndex[i]=dummy;
    aa_t type=static_cast<aa_t>(p.getResidue(i)->getID());
    if(type==GLY)cbSideIndex[i]=-1;
    else if(type==ALA)dummy++;
    else if(type==ILE||type==VAL||type==THR) dummy+=3;
    else dummy+=2;
  }
 
  /*print the txt*/
  printSYS_SIZE(cout, size);
  printNUM_ATOMS(cout, p);
  printATOM_TYPE(cout, p);
  printNONEL_COL(cout, p, cbSideIndex);
  printBOND_TYPE(cout, p, cbSideIndex);
  printREACT(cout);
  printATOM_LIST(cout, p, ran, cbSideIndex);
  printBOND_LIST(cout, p, cbSideIndex);
  printPERM_BOND_LIST(cout, p, cbSideIndex);
  printHBA_LIST(cout, p);
}
