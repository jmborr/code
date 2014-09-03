#include "PDBLib.h"
#include "pseudoAA.h"
using namespace std;

/*------------------------------------------------------------------
  SER, THR, ASN, ASP, GLN, GLU contains polar Oxygen, which can form
  hydrogen bond with backbone Nitrogen.
  ----------------*/

const static double ndn   =  1.4872 ;
const static double ndc   = -0.7784 ;
const static double ndca  = -0.7088 ;
const static double N_H_d = 1.02;
const static double cdc   =  1.3379 ;
const static double cdca  = -0.6144 ;
const static double cdn   = -0.7235 ;
const static double C_O_d = 1.24;

int isHBonded(aa& iaa, aa& pre_iaa, aa& jaa, aa& pre_jaa,
	      double& N_O, double& N_C, double& O_CA, double& O_C, int& sep){
  /*Determine whethe the two peptide planes are hbonded*/
  vec* Ni = iaa.getN()->getR();
  vec* Ci = pre_iaa.getC()->getR();
  vec* Cai= iaa.getAtomAt(0)->getR();
  vec NH = ndn*(*Ni) + ndc*(*Ci) + ndca*(*Cai);
  NH *= N_H_d;
  vec Hi = *Ni + NH;
  
  vec* Cj = pre_jaa.getC()->getR();
  vec* Caj= pre_jaa.getAtomAt(0)->getR();
  vec* Nj = jaa.getN()->getR();
  atom* O_atom = pre_jaa.getO();
  vec* Oj; 
  vec CO;
  if(!O_atom){
    CO = cdc*(*Cj) + cdn*(*Nj) + cdca*(*Caj);
    CO *= C_O_d;
    Oj = new vec(*Cj + CO);
  }
  else{
    Oj  = O_atom->getR();
    CO = (*Oj) - (*Cj);
  }
  
  
  int isHBonded=0;
  if(Oj->getDist(Hi) < 2.5){
    if((*Oj - Hi)*NH>0 && (Hi - *Oj)*CO>0)isHBonded=1;
  } 

  if(isHBonded){
    N_O = Oj->getDist(*Ni);
    N_C = Cj->getDist(*Ni);
    O_CA= Oj->getDist(*Cai);
    O_C = Oj->getDist(*Ci);
    sep = +1;
    return 1;
  }
  
  Ni = jaa.getN()->getR();
  Ci = pre_jaa.getC()->getR();
  Cai= jaa.getAtomAt(0)->getR();
  NH = ndn*(*Ni) + ndc*(*Ci) + ndca*(*Cai);
  NH *= N_H_d;
  Hi = *Ni + NH;
  
  Cj = pre_iaa.getC()->getR();
  Caj= pre_iaa.getAtomAt(0)->getR();
  Nj = iaa.getN()->getR();
  O_atom  = pre_iaa.getO();
  if(!O_atom){
    CO = cdc*(*Cj) + cdn*(*Nj) + cdca*(*Caj);
    CO *= C_O_d;
    Oj = new vec(*Cj + CO);
  }
  else{
    Oj  = O_atom->getR();
    CO = (*Oj) - (*Cj);
  }

  if(Oj->getDist(Hi) < 2.5){
    if((*Oj - Hi)*NH>0 && (Hi - *Oj)*CO>0)isHBonded=1;
  }

  if(isHBonded){
    N_O = Oj->getDist(*Ni);
    N_C = Cj->getDist(*Ni);
    O_CA= Oj->getDist(*Cai);
    O_C = Oj->getDist(*Ci);
    sep = -1;
    return 1;
  }
  
  return 0;
  
}


int isHBonded_NCAP(aa& acceptor, aa& donor, aa& prev, 
	      double& N_PO, double& C_PO, double& CA_PO){
  aa_t id=aa_Name2Type(cstr(acceptor.getName()).getStr());
  vec* N  = donor.getN()->getR();
  vec* Ca = donor.getAtomAt(0)->getR();
  vec* C  = prev.getC()->getR();
  vec NH = ndn*(*N) + ndc*(*C) + ndca*(*Ca);
  NH *= N_H_d;
  vec H = *N + NH;
  atom* O1 = NULL;
  atom* O2 = NULL;
  vec PO;
  pseudo_aa a(acceptor);
  switch(id){
  case SER: //OG
    O1 = acceptor.getAtomAt(2);
    PO = *O1->getR();
    break;
  case THR: //OG1
    O1 = acceptor.getAtomAt(2);
    if(O1->getName()!='O')//OG2?
      O1 = acceptor.getAtomAt(3);
    PO = *O1->getR();
    break;
  case ASN: //OD1
    O1 = acceptor.getAtomAt(3);
    if(O1->getName()!='O')//OD2?
      O1 = acceptor.getAtomAt(4);
    PO = *a.getSide();
    break;
  case ASP: //OD1 OD2
    O1 = acceptor.getAtomAt(3);
    O2 = acceptor.getAtomAt(4);
    PO = *a.getSide();
    break;
  case GLN: //OE1
    O1 = acceptor.getAtomAt(4);
    if(O1->getName()!='O')//OE2?
      O1 = acceptor.getAtomAt(5);
    PO = *acceptor.getAtomAt(3)->getR();
    break;
  case GLU: //OE1 OE2
    O1 = acceptor.getAtomAt(4);
    O2 = acceptor.getAtomAt(5);
    PO = *acceptor.getAtomAt(3)->getR();
    break;
  case TYR: //OH
    O1 = acceptor.getAtomAt(8);
    PO = *a.getSide();
    break;
  case HIS://ND1
    O1 = acceptor.getAtomAt(6);
    if(O1->getName()!='N')//ND2
      O1 = acceptor.getAtomAt(5);
    PO = *a.getSide();
  default:
    break;
  }
  int isBonded=0;
  if(O1 && O1->getR()->getDist(H) < 2.5){
    if((*O1->getR()-H)*NH > 0) isBonded=1;
  }
  else if(O2 && O2->getR()->getDist(H) < 2.5){
    if((*O2->getR()-H)*NH > 0) isBonded=1;
  }
  if(isBonded){
    N_PO = PO.getDist(*N);
    C_PO = PO.getDist(*C);
    CA_PO= PO.getDist(*Ca);
  }
  return isBonded;
}

int isHBonded_CCAP(aa& donor, aa& acceptor, aa& prev, 
		   double& S_O, double& S_C){
  aa_t id=aa_Name2Type(cstr(donor.getName()).getStr());
  vec* N  = acceptor.getN()->getR();
  vec* Ca = prev.getAtomAt(0)->getR();
  vec* C  = prev.getC()->getR();
  vec CO = cdn*(*N) + cdc*(*C) + cdca*(*Ca);
  CO *= C_O_d;
  vec O = *C + CO;
  
  atom* D1 = NULL;
  //atom* D2 = NULL;
  vec S;
  pseudo_aa a(donor);
  switch(id){
  case SER: //OG
    D1 = donor.getAtomAt(2);
    S = *D1->getR();
    break;
  case THR: //OG1
    D1 = donor.getAtomAt(2);
    if(D1->getName()!='O')//OG2?
      D1 = donor.getAtomAt(3);
    S = *D1->getR();
    break;
  case ASN: //ND1
    D1 = donor.getAtomAt(4);
    if(D1->getName()!='N')//ND2?
      D1 = donor.getAtomAt(3);
    S = *a.getSide();
    break;
  case GLN: //NE1
    D1 = donor.getAtomAt(5);
    if(D1->getName()!='N')//NE2?
      D1 = donor.getAtomAt(4);
    S = *donor.getAtomAt(3)->getR();
    break;
  case TYR: //OH
    D1 = donor.getAtomAt(8);
    S = *a.getSide();
    break;
  case TRP://NE1
    D1 = donor.getAtomAt(5);
    S = *a.getSide();
    break;
  case HIS://NE2
    D1 = donor.getAtomAt(3);
    if(D1->getName()!='N')//NE1?
      D1 = donor.getAtomAt(4);
    S = *a.getSide();
    break;
  default:
    break;
  }

  int isBonded=0;
  if(D1 && D1->getR()->getDist(O) < 3.5){
    if((*D1->getR()-O)*CO > 0) 
      isBonded=1;
  }
  if(isBonded){
    S_O = S.getDist(O);
    S_C = S.getDist(*C);
  }
  return isBonded;
}

int main(int argc, char* argv[]){
  if(argc<2){
    cout << "com opdb" << endl;
    exit(1);
  }
  protein p(argv[1]);
  double n_dist, c_dist, ca_dist;
  double n_o, n_c, o_ca, o_c;
  int sep;
  double S_O, S_C;
  for(int i=1; i<p.getLength(); i++){
    aa* curr = p.getElementAt(i);
    aa* prev = p.getElementAt(i-1);
    for(int j=1; j<i-1; j++){
      aa* theAA = p.getElementAt(j);
      aa* preAA = p.getElementAt(j-1);
      if(isHBonded_NCAP(*curr, *theAA, *preAA, n_dist, c_dist, ca_dist)){
	cout <<"N-Cap: " <<  curr->getName() << i+1 <<  " " << j+1 << " " 
	     << n_dist << " " <<  c_dist << " " << ca_dist << endl;
      }
      if(isHBonded_NCAP(*theAA, *curr, *prev, n_dist, c_dist, ca_dist)){
	cout <<"N-Cap: " <<  theAA->getName() << j+1 <<  " " << i+1 << " " 
	     << n_dist << " " <<  c_dist << " " << ca_dist << endl;
      }
      if(isHBonded_CCAP(*curr, *theAA, *preAA, S_O, S_C))
	cout << "C-Cap: " << curr->getName() << i+1 << " " << j << " " 
	     << S_O << " " << S_C << endl;
      if(isHBonded_CCAP(*theAA, *curr, *prev, S_O, S_C))
	cout << "C-Cap: " << theAA->getName() << j+1 << " " << i << " "
	     << S_O << " " << S_C << endl;
      if(isHBonded(*curr, *prev, *theAA, *preAA, n_o, n_c, o_ca, o_c, sep)){
	if(sep==1){
	  printf("%4ld %4ld %4ld %10.6lf %10.6lf %10.6lf %10.6lf\n", i+1, j, i-j+sep, n_o, n_c, o_ca, o_c);
	}
	else{
	  printf("%4ld %4ld %4ld %10.6lf %10.6lf %10.6lf %10.6lf\n", i, j+1, i-j+sep, n_o, n_c, o_ca, o_c);
	}
      }
    }
  }
}
