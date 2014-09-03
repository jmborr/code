/*****************************************
The peptide is composed of amino acids.
*****************************************/

#ifndef _PSEUDOPEP_H_
#define _PSEUDOPEP_H_

#include "pseudoAA.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

/*parameter to define vector CO, NH*/
const static double ndn   =  1.4872 ;
const static double ndc   = -0.7784 ;
const static double ndca  = -0.7088 ;
const static double cdc   =  1.3379 ;
const static double cdca  = -0.6144 ;
const static double cdn   = -0.7235 ;

//comments line: start with "//" or "#"
int isComments(string line){
  if(!line.length()) return -1; //blank line
  int non_whitespace=0;
  for(int i=0; i<line.length(); i++){
    if(line.substr(i,1)!=" " && line.substr(i,1)!="\t") break;
    non_whitespace++;
  }
  if(line.length()==non_whitespace) return -1; //blank line also
  
  if(line.substr(non_whitespace,1)=="#") return 1;
  else if(line.substr(non_whitespace,2)=="//") return 1;
  return 0;
}

inline void str2chars(const string line, char* buf){
  line.copy(buf, line.length());
  buf[line.length()]='\0';
}

typedef enum{_SEQ_=0, _PDB_SEQ_, _PDB_BB_, _CRYSTAL_}input_t;

class pseudoPep{
 protected:
  vector<pseudoAA*> chain;
  int n_dmd_atoms;
 public:
  pseudoPep(const char* inpt_f, input_t t){
    pseudoAA* newAA;
    pseudoAA* prev=NULL;
    if(t==_SEQ_){
      ifstream in(inpt_f);
      string line;
      char name[80];
      while(getline(in,line)){
	if(!isComments(line)){
	  str2chars(line, name);
	  newAA = new pseudoAA(name);
	  if(prev){
	    if(newAA->getID()==PRO)
	      newAA->joinTo(*prev, 120, -60);
	    else{
	      newAA->joinTo(*prev, 120, -120);
	      //newAA->joinTo(*prev, -45, -60); HELIX
	    }
	  }
	  chain.push_back(newAA);
	  prev = newAA;
	}
      }
      in.close();
    }
    else if(t==_PDB_SEQ_){
      protein prot(inpt_f);
      char name[80];
      for(int i=0; i<prot.getLength(); i++){
	str2chars(prot.getElementAt(i)->getName(), name);
	newAA = new pseudoAA(name);
	if(prev){
	  if(newAA->getID()==PRO)
	    newAA->joinTo(*prev, 120, -60);
	  else
	    newAA->joinTo(*prev, 120, -120);
	}
	chain.push_back(newAA);
	prev = newAA;
      }
    }
    else if(t==_PDB_BB_){
      pseudoProtein prot(inpt_f);
      double phi, psi, prev_psi, chi;
      for(int i=0; i<prot.getLength(); i++){
	aa* theAA = prot.getElementAt(i);
	newAA = new pseudoAA(*theAA);
	prot.getPhiPsi(i, phi, psi);
	if(prot.getPseudoChi(i, chi) && theAA->getName()!="PRO"){
	  if(newAA->getNSideChains()==3){
	    chi -= 60.0*rpi;
	  }
	  newAA->rotateChi(chi);
	}
	if(prev){
	  phi = (phi-PI)/rpi;
	  prev_psi = (prev_psi-PI)/rpi;
	  prev->setPsi(prev_psi);
	  newAA->setPhi(phi);
	}
	chain.push_back(newAA);
	prev = newAA;
	prev_psi = psi;
      }
    }
    else if(t==_CRYSTAL_){
      pseudoProtein prot(inpt_f);
      double phi, psi, prev_psi;
      for(int i=0; i<prot.getLength(); i++){
	aa* theAA = prot.getElementAt(i);
	pseudo_aa* thePseudoAA = prot.getPseudoElementAt(i);
	newAA = new pseudoAA(*theAA, *thePseudoAA);
	prot.getPhiPsi(i, phi, psi);
	if(prev){
	  phi = (phi-PI)/rpi;
	  prev_psi = (prev_psi-PI)/rpi;
	  prev->setPsi(prev_psi);
	  newAA->setPhi(phi);
	}
	chain.push_back(newAA);
	prev = newAA;
	prev_psi = psi;
      }
    }
  }

  ~pseudoPep(){
    for(int i=0; i<chain.size(); i++){
      delete chain[i];
    }
    chain.clear();
  }

  pseudoAA* getResidue(int i){
    if(i>=0 && i<chain.size())
      return chain[i];
    return NULL;
  }

  int getLength(){
    return chain.size();
  }

  void constructOxygen(){
    for(int i=0; i<chain.size()-1; i++){
      vec c_o = (*chain[i]->getCA()->getR())*cdca 
	+ (*chain[i]->getC()->getR())*cdc 
	+ (*chain[i+1]->getN()->getR())*cdn;
      c_o *= C_O;
      c_o += *chain[i]->getC()->getR();
      *chain[i]->getO()->getR() = c_o;
    }
    /*last CO is put in psi=-180*/
    int last = chain.size();
    vec ca_n = (*chain[last-1]->getN()->getR()) - (*chain[last-1]->getCA()->getR());
    ca_n.Normalize();
    vec ca_c = (*chain[last-1]->getC()->getR()) - (*chain[last-1]->getCA()->getR());
    ca_c.Normalize();
    vec c_o = ca_n*sin(CA_C_O*rpi)/sin(N_CA_C*rpi)
      + ca_c*sin((CA_C_O+N_CA_C-180)*rpi)/sin(N_CA_C*rpi);
    c_o *= C_O;
    c_o += *chain[last-1]->getC()->getR();
    *chain[last-1]->getO()->getR() = c_o;
  }

  void constructHydrogen(){
    for(int i=1; i<chain.size(); i++){
      if(chain[i]->getID()!=PRO){
	vec n_h = (*chain[i]->getCA()->getR())*ndca
	  + (*chain[i-1]->getC()->getR())*ndc
	  + (*chain[i]->getN()->getR())*ndn;
	n_h *= N_H;
	n_h += *chain[i]->getN()->getR();
	if(chain[i]->getH()){
	  *chain[i]->getH()->getR() = n_h;
	}
	else{
	  chain[i]->createH(n_h);
	}
      }
    }
    /*first NH3 is not included*/
  }

  void printPDB(const char* pdb_f){
    ofstream out(pdb_f);
    int atomIndex=1;
    for(int i=0; i<chain.size(); i++)
      chain[i]->printPDB(out, i+1, atomIndex);
    out.close();
  }
  
  void printPDB(ostream& out){
    int atomIndex=1;
    for(int i=0; i<chain.size(); i++)
      chain[i]->printPDB(out, i+1, atomIndex);
  }
  
  int makeIndex(){
    int index=1;
    for(int i=0; i<chain.size(); i++){
      chain[i]->getN()->setIndex(index++);
      chain[i]->getC()->setIndex(index++);
      chain[i]->getO()->setIndex(index++);
      chain[i]->getCA()->setIndex(index++);
      for(int j=0; j<chain[i]->getNSideChains(); j++){
	chain[i]->getSideChains()[j].setIndex(index++);
      }
    }
    n_dmd_atoms = index-1;
    return index-1;
  }
  
  void shift(vec& r){
    for(int i=0; i<chain.size(); i++){
      chain[i]->shift(r);
    }
  }

  void printCMap(ofstream& out){
    int contact=0;
    for(int i=0; i<chain.size(); i++){
      pseudoAA* theAA = chain[i];
      atom* CBi = theAA->getCB();
      double CBir = 2.0;
      atom* CGi = theAA->getG();
      double CGir = 0;
      if(CGi)CGir = G1_CONST[theAA->getID()][G1_IR];
      atom* CG2i = theAA->getG2();
      double CG2ir = 0;
      if(CG2i){
	if(theAA->getID()==ILE) CG2ir = ILE_G2[G2_IR];
	else if(theAA->getID()==THR) CG2ir = THR_G2[G2_IR];
	else if(theAA->getID()==VAL) CG2ir = VAL_G2[G2_IR];
      }
      for(int j=i+1; j<chain.size(); j++){
	pseudoAA* nxtAA = chain[j];
	atom* CBj = nxtAA->getCB();
	double CBjr = 2.0;
	atom* CGj = nxtAA->getG();
	double CGjr = 0;
	if(CGj) CGjr = G1_CONST[nxtAA->getID()][G1_IR];
	atom* CG2j = nxtAA->getG2();
	double CG2jr = 0;
	if(CG2j){
	  if(nxtAA->getID()==ILE) CG2jr = ILE_G2[G2_IR];
	  else if(nxtAA->getID()==THR) CG2jr = THR_G2[G2_IR];
	  else if(nxtAA->getID()==VAL) CG2jr = VAL_G2[G2_IR];
	}
	if(!contact && CBi && CGj && CBi->getDist(*CGj) < (CBir+CGjr+1.0)) contact=1;
	if(!contact && CBj && CGi && CBj->getDist(*CGi) < (CBjr+CGir+1.0)) contact=1;
	if(!contact && CBi && CG2j && CBi->getDist(*CG2j) < (CBir+CG2jr+1.0)) contact=1;
	if(!contact && CBj && CG2i && CBj->getDist(*CG2i) < (CBjr+CG2ir+1.0)) contact=1;
	if(!contact && CGi && CGj && CGi->getDist(*CGj) < (CGir+CGjr+1.0)) contact=1;
	if(!contact && CGj && CGi && CGj->getDist(*CGi) < (CGjr+CGir+1.0)) contact=1;
	if(!contact && CGi && CG2j && CGi->getDist(*CG2j) < (CGir+CG2jr+1.0)) contact=1;
	if(!contact && CGj && CG2i && CGj->getDist(*CG2i) < (CGjr+CG2ir+1.0)) contact=1;
	if(!contact && CG2i && CG2j && CG2i->getDist(*CG2j) < (CG2ir+CG2jr+1.0)) contact=1;
	if(!contact && CG2j && CG2i && CG2j->getDist(*CG2i) < (CG2jr+CG2ir+1.0)) contact=1;
	if(contact){
	  out << i << " " << j << endl;
	  contact=0;
	}
      }
    }
  }
};

#endif/*_PSEUDO_PEP_H_*/
