#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include "pseudoPep.h"
#include "pseudoTXT.h"
#include "PDBLib.h"
#include "eMatrix.h"
#include "random.h"
using namespace std;

typedef struct {
  dmd_atom_t t;
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
} pseudo_t;

typedef struct {
  int p;
  int q;
} bond_t;

int isKeywords(string line){
  //cout << sizeof(txtKeyWords)/sizeof(string) << endl;
  for(int i=0; i<sizeof(txtKeyWords)/sizeof(string); i++){
    if(txtKeyWords[i]==line) return i+1;
  }
  return 0;
}

void processTXT(char* txt, vector<pseudo_t *>& atoms, vector<bond_t *>& hbonds){
  ifstream in(txt);
  string line;
  int theKey=0;
  int tmpKey=0;
  int nTotal=0;
  char buf[1024];
  int iatom, itype;
  double x,y,z,vx,vy,vz,boundary;
  int p,q;
  vector<bond_t*> all_bonds;
  vector<bond_t*> perm_bonds;
  pseudo_t* a;
  bond_t* b;
  while(getline(in, line)){
    tmpKey=isKeywords(line);
    if(theKey){
      if(!tmpKey){
        //not a key
        if(!isComments(line)){
          //not a comment --> data entry line
          switch(theKey){
          case 1://size
	    str2ch(line, buf);
            sscanf(buf, "%lf", &boundary);
            break;
          case 2://NUMBER OF ATOMS
            break;
          case 3://type of atoms
            break;
          case 4: //non-elastic collisions
          case 5: //elastic collisions
          case 6: //bond
          case 7: //reaction
            break;
          case 8://list of atoms
	    str2ch(line, buf);
	    sscanf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf",
		   &iatom, &itype, &x, &y, &z, 
		   &vx, &vy, &vz);
	    a = (pseudo_t*)malloc(sizeof(pseudo_t));
	    a->t=static_cast<dmd_atom_t>(itype);
	    a->x=x;a->y=y;a->z=z;
	    a->vx=vx;a->vy=vy;a->vz=vz;
	    atoms.push_back(a);
            break;
          case 9:
	    str2ch(line,buf);
	    sscanf(buf, "%ld %ld", &p, &q);
	    b = (bond_t*)malloc(sizeof(bond_t));
	    b->p=p; b->q=q;
	    all_bonds.push_back(b);
            break;
	  case 10:
	    str2ch(line,buf);
	    sscanf(buf, "%ld %ld", &p, &q);
	    b = (bond_t*)malloc(sizeof(bond_t));
	    b->p=p;b->q=q;
	    perm_bonds.push_back(b);
	    break;
          default:
	    break;
	  }
	}
      }
      else if(theKey!=tmpKey){
	theKey=tmpKey;
      }
    }
    else if(tmpKey){
      theKey=tmpKey;
    }
  }
  in.close();

  /*considering PBC*/
  for(int i=1; i<atoms.size(); i++){
    if(atoms[i]->x-atoms[i-1]->x > boundary/2.0f ) atoms[i]->x-=boundary;
    if(atoms[i]->x-atoms[i-1]->x < -boundary/2.0f) atoms[i]->x+=boundary;
    if(atoms[i]->y-atoms[i-1]->y > boundary/2.0f ) atoms[i]->y-=boundary;
    if(atoms[i]->y-atoms[i-1]->y < -boundary/2.0f) atoms[i]->y+=boundary;
    if(atoms[i]->z-atoms[i-1]->z > boundary/2.0f ) atoms[i]->z-=boundary;
    if(atoms[i]->z-atoms[i-1]->z < -boundary/2.0f) atoms[i]->z+=boundary;
  }

  /*filter the permbond from allbond*/
  int p1,q1;
  int found;
  for(int i=0; i<all_bonds.size(); i++){
    p=all_bonds[i]->p;
    q=all_bonds[i]->q;
    found=0;
    for(int j=0; j<perm_bonds.size(); j++){
      p1=perm_bonds[j]->p;
      q1=perm_bonds[j]->q;
      if((p==p1&&q==q1)||(p==q1&&q==p1)){
	found=1;
	break;
      }
    }
    if(!found) hbonds.push_back(all_bonds[i]);
    else delete all_bonds[i];
  }
  all_bonds.clear();
  for(int i=0; i<perm_bonds.size(); i++) delete perm_bonds[i];
  perm_bonds.clear();
}

int isContact(atom* cb1, atom* cb2, double hdr, double ir){
  double dist = cb1->getDist(*cb2);
  if(dist<ir)
    return 1;
  else 
    return 0;
}

int main(int argc, char* argv[]){
  if(argc<4){
    cout << "usage: readTXT.linux txt inputMethod inputFile" << endl;
    cout << " inputMethod : 0/1/2/3" << endl;
    cout << "    0 --- input seq file only" << endl;
    cout << "    1 --- input PDB file but using SEQ only" << endl;
    cout << "       these two case, the generated peptide is streched" << endl;
    cout << "    2 --- input PDB file and fix BACKBONE and use approximate SideChain INFOR" << endl;
    cout << "       the sidechain geometry is satisfied by the DESIGN" << endl;
    cout << "    3 --- input PDB file and fix BACKBONE and SIDECHAIN from PDB" << endl;
    cout << "       the sidechain geometry might be away from DESIGN, need relax" << endl;
    cout << " inputFile : the corresponding file(SEQ or PDB) according to inputMethod" << endl;
    exit(1);
  }
  vector<pseudo_t *> atoms;
  vector<bond_t *> hbonds;  
  /*read TXT*/
  processTXT(argv[1], atoms, hbonds);
  /*create the protein*/
  pseudoPep p(argv[3], (input_t)atoi(argv[2]));
  pseudoAA* lastAA=p.getResidue(p.getLength()-1);
  lastAA->getO()->getR()->Rotate(*lastAA->getCA()->getR(),*lastAA->getC()->getR(), -PI/3.0);
  int nAtom = p.makeIndex();
  if(nAtom!=atoms.size()){
    cout << "The txt and pdb files are not matching" << endl;
    exit(1);
  }
  /*reset the protein*/
  int index=0;
  int* array= new int[nAtom](-1);;
  for(int iAA=0; iAA<p.getLength(); iAA++){
    pseudoAA* theAA = p.getResidue(iAA);
    *theAA->getN()->getR()  = vec(atoms[index]->x, atoms[index]->y, atoms[index]->z);
    array[index]=iAA+1;
    index++;
    *theAA->getC()->getR()  = vec(atoms[index]->x, atoms[index]->y, atoms[index]->z);
    index++;
    *theAA->getO()->getR()  = vec(atoms[index]->x, atoms[index]->y, atoms[index]->z);
    array[index]=iAA+1;
    index++;
    *theAA->getCA()->getR() = vec(atoms[index]->x, atoms[index]->y, atoms[index]->z);
    index++;
    if(theAA->getCB()){
      *theAA->getCB()->getR() = vec(atoms[index]->x, atoms[index]->y, atoms[index]->z);
      index++;
    }
    if(theAA->getG()){
      *theAA->getG()->getR()  = vec(atoms[index]->x, atoms[index]->y, atoms[index]->z);
      index++;
    }
    if(theAA->getG2()){
      *theAA->getG2()->getR() = vec(atoms[index]->x, atoms[index]->y, atoms[index]->z);
      index++;
    }
  }
  /*print out the hydrogen bonds*/
  int ihbonds=0;
  for(int i=0; i<hbonds.size(); i++){
    int p=hbonds[i]->p;
    int q=hbonds[i]->q;
    if((atoms[p-1]->t==_O_HB_ && atoms[q-1]->t==_N_HB_)||
       (atoms[p-1]->t==_N_HB_ && atoms[q-1]->t==_O_HB_)){
      cout << array[p-1] << " " << array[q-1] << endl;
      ihbonds++;
    }
  }
  cout << "&" << endl;
  /*print out the contact maps*/
  for(int i=0; i<p.getLength(); i++){
    aa_t type=static_cast<aa_t>(p.getResidue(i)->getID());
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
    
    for(int j=i+2; j<p.getLength(); j++){
      aa_t typep=static_cast<aa_t>(p.getResidue(j)->getID());
      atom* CBp=NULL;
      atom* CGp=NULL;
      atom* CG2p=NULL;
      double hdr, ir;
      int contact=0;
      if(typep!=GLY){
	CBp=p.getResidue(j)->getCB();
	hdr=CB_HDR;
	ir =CB_IR;
	if(CB){
	  if(isContact(CB, CBp, hdr+cb_hdr, ir+cb_ir))
	    if(!contact)contact=1;
	}
	if(CG){
	  if(isContact(CG, CBp, hdr+cg_hdr, ir+cg_ir))
	    if(!contact)contact=1;
	}
	if(CG2){
	  if(isContact(CG2, CBp, hdr+cg2_hdr, ir+cg2_ir))
	    if(!contact)contact=1;
	}
	if(typep!=ALA){
	  CGp=p.getResidue(j)->getG();
	  hdr=G1_CONST[typep][G1_HDR];
	  ir=G1_CONST[typep][G1_IR];
	  if(CB){
	    if(isContact(CB, CGp, hdr+cb_hdr, ir+cb_ir))
	      if(!contact)contact=1;
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
	      if(isContact(CG, CGp, thdr, ir+cg_ir))
		if(!contact)contact=1;
	    }
	    else{
	      if(isContact(CG, CGp, hdr+cg_hdr, ir+cg_ir))
		if(!contact)contact=1;
	    }
	  }
	  if(CG2){
	    if(isContact(CG2, CGp, hdr+cg2_hdr, ir+cg2_ir))
	      if(!contact)contact=1;
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
	      if(isContact(CB, CG2p, hdr+cb_hdr, ir+cb_ir))
		if(!contact)contact=1;
	    }
	    if(CG){
	      if(isContact(CG, CG2p, hdr+cg_hdr, ir+cg_ir))
		if(!contact)contact=1;
	    }
	    if(CG2){
	      if(isContact(CG2, CG2p, hdr+cg2_hdr, ir+cg2_ir))
		if(!contact)contact=1;
	    }
	  }
	}
	
      }
      if(contact) cout << i+1 << " " << j+1 << endl;
    }
  }
}

