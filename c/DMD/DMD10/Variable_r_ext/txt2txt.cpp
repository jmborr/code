#include "pseudoAA.h"
#include "pseudoPep.h"
#include "PDBLib.h"
#include "random.h"
#include "pseudoTXT.h"
#include <cstring>
using namespace std;

int isKeywords(string line){
  //cout << sizeof(txtKeyWords)/sizeof(string) << endl;                                                      
  for(int i=0; i<sizeof(txtKeyWords)/sizeof(string); i++){
    if(txtKeyWords[i]==line) return i+1;
  }
  return 0;
}

inline void str2ch(const string line, char* buf){
  line.copy(buf, line.length());
  buf[line.length()]='\0';
}

void process_txt(const char* txt, double& size, vector<double>& xyz){
  ifstream in(txt);
  string line;
  int theKey=0;
  int tmpKey=0;
  char buf[1024];
  int iatom, itype;
  double x,y,z;
  double vx, vy, vz;

  while(getline(in, line)){
    tmpKey=isKeywords(line);
    if(theKey){
      if(!tmpKey){
        //not a key                                                                                          
        if(!isComments(line)){
          //not a comment --> data entry line                                                                
          switch(theKey){
          case 1://SYSTEM SIZE                                                                               
            str2ch(line, buf);
            sscanf(buf, "%lf", &size);//read size                                                            
            break;
          case 8://list of atoms                                                                             
            str2ch(line, buf);
            sscanf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf",
                   &iatom, &itype, &x, &y, &z, &vx, &vy, &vz);
	    xyz.push_back(x);
	    xyz.push_back(y);
	    xyz.push_back(z);
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
}


void single_bond_sprint(char* buf, double ave, double dev){
  sprintf(buf, "%lf %lf", ave*(1.0-dev), ave*(1.0+dev));
}

void single_bond_sprint_len(char* buf, double ave, double len_dev){
  sprintf(buf, "%lf %lf", ave-len_dev, ave+len_dev);
}

void double_bond_sprint(char* buf, double l_ave, double l_dev, 
			double r_ave, double r_dev){
  if(l_ave*(1.0+l_dev) > r_ave*(1.0-r_dev)){
    cerr << "Error: overlaping bonds" << endl;
    exit(1);
  }
  sprintf(buf, "%lf %lf %lf %lf %lf %lf", 
	  l_ave*(1.0-l_dev), l_ave*(1.0+l_dev), -max_height, 
	  r_ave*(1.0-r_dev),  max_height, r_ave*(1.0+r_dev));
}

void type_sprint(char* buf, dmd_atom_t type, double mass, double hdr, double ir){
  if(ir==0) ir = hdr + eps;
  sprintf(buf, "%ld %lf %lf %lf", type, mass, hdr, ir);
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

void el_col_sprint(char* buf, dmd_atom_t a, dmd_atom_t b, double min_d){
  sprintf(buf, "%ld %ld %lf", a, b, min_d);
}

void sel_col_sprint(char* buf, dmd_atom_t a, dmd_atom_t b, double hdr, double min_d, 
		    int nstep=1, double e=Dihedral_E){
  double delta = min_d - hdr;
  double step = delta/nstep;
  int pt = sprintf(buf, "%ld %ld %lf", a, b, hdr);
  double dist = hdr;
  for(int i=0; i<nstep; i++){
    dist += step;
    pt += sprintf(&buf[pt], " %lf %lf", dist, e);
  }
}

void hp_sprint(char* buf, dmd_atom_t a, dmd_atom_t b,
	       double hd, double id, double e){
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf",
	    a, b, hd, id, -fabs(e), id+ir_ext, -fabs(e)/2.0);
  }
  else{ sprintf(buf, "%ld %ld %lf %lf %lf", a, b, hd, id, -fabs(e));}
}

void sb_sprint(char* buf, dmd_atom_t a, dmd_atom_t b,
	       double hd, double id, double e){
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf %lf",
	    a, b, hd, id-ir_ext, -fabs(e), id, -fabs(e)/4.0, id+ir_ext,
	    -fabs(e)/4.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", a, b, hd, id, -fabs(e));
  }
}

void rsb_sprint(char* buf, dmd_atom_t a, dmd_atom_t b,
		double hd, double id, double e){
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf %lf",
	    a, b, hd, id-ir_ext, fabs(e), id, fabs(e)/4.0,
	    id+ir_ext, fabs(e)/4.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", a, b, hd, id, fabs(e));
  }
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
  else if(lastAtom=lastAA->getD());
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

void printATOM_TYPE(ostream& out){
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

  /*CB atoms*/
  int type_shift = _CYS_CB_;
  for(int i=0; i<20; i++){
    dmd_atom_t cb_t = static_cast<dmd_atom_t>(type_shift + i);
    if(i!=GLY){
      type_sprint(buf, cb_t, CB_M, CB_HDR, 0.0);
    }
    else{/*gly CB is a never-used type*/
      type_sprint(buf, cb_t, 1.0, eps, eps);
    }
    out << buf << endl;
  }
  
  /*CG atoms*/
  type_shift = _CYS_CG_;
  for(int i=0; i<20; i++){
    dmd_atom_t cg_t = static_cast<dmd_atom_t>(type_shift + i);
    if(i!=ALA && i!=GLY){
      const double* cg_para = G1_CONST[i];
      type_sprint(buf, cg_t, cg_para[G1_M], cg_para[G1_HDR], cg_para[G1_IR]);
    }
    else{/*ala, gly CG is a never-used type*/
      type_sprint(buf, cg_t, 1.0, eps, eps);
    }
    out << buf << endl;
  }

  /*CG2 atoms*/
  /*ILE*/
  type_sprint(buf, _ILE_CG2_, ILE_G2[G2_M], ILE_G2[G2_HDR], ILE_G2[G2_IR]);
  out << buf << endl;
  /*THR*/
  type_sprint(buf, _THR_CG2_, THR_G2[G2_M], THR_G2[G2_HDR], THR_G2[G2_IR]);
  out << buf << endl;
  /*VAL*/
  type_sprint(buf, _VAL_CG2_, VAL_G2[G2_M], VAL_G2[G2_HDR], VAL_G2[G2_IR]);
  out << buf << endl;

  /*CD atoms*/
  /*ARG*/
  type_sprint(buf, _ARG_CD_, ARG_D[D_M], ARG_D[D_HDR], ARG_D[D_IR]);
  out << buf << endl;
  /*LYS*/
  type_sprint(buf, _LYS_CD_, LYS_D[D_M], LYS_D[D_HDR], LYS_D[D_IR]);
  out << buf << endl;
  /*TRP*/
  type_sprint(buf, _TRP_CD_, TRP_D[D_M], TRP_D[D_HDR], TRP_D[D_IR]);
  out << buf << endl;
  
  /*DERIVED HBA types*/
  type_sprint(buf, _ASP_HBA_, 
	      G1_CONST[ASP][G1_M], G1_CONST[ASP][G1_HDR], G1_CONST[ASP][G1_IR]);
  out << buf << endl;
  type_sprint(buf, _ASN_HBA_, 
	      G1_CONST[ASN][G1_M], G1_CONST[ASN][G1_HDR], G1_CONST[ASN][G1_IR]);
  out << buf << endl;
  type_sprint(buf, _GLN_HBA_, 
	      G1_CONST[GLN][G1_M], G1_CONST[GLN][G1_HDR], G1_CONST[GLN][G1_IR]);
  out << buf << endl;
  type_sprint(buf, _GLU_HBA_, 
	      G1_CONST[GLU][G1_M], G1_CONST[GLU][G1_HDR], G1_CONST[GLU][G1_IR]);
  out << buf << endl;
  type_sprint(buf, _SER_HBA_, 
	      G1_CONST[SER][G1_M], G1_CONST[SER][G1_HDR], G1_CONST[SER][G1_IR]);
  out << buf << endl;
  type_sprint(buf, _THR_HBA_, 
	      THR_G2[G2_M], THR_G2[G2_HDR], THR_G2[G2_IR]);
  out << buf << endl;
  /*derived HBD types*/
  type_sprint(buf, _ASN_HBD_,
	      G1_CONST[ASN][G1_M], G1_CONST[ASN][G1_HDR], G1_CONST[ASN][G1_IR]);
  out << buf << endl;
  type_sprint(buf, _GLN_HBD_,
	      G1_CONST[GLN][G1_M], G1_CONST[GLN][G1_HDR], G1_CONST[GLN][G1_IR]);
  out << buf << endl;
  type_sprint(buf, _SER_HBD_,
	      G1_CONST[SER][G1_M], G1_CONST[SER][G1_HDR], G1_CONST[SER][G1_IR]);
  out << buf << endl;
  type_sprint(buf, _THR_HBD_,
	      THR_G2[G2_M], THR_G2[G2_HDR], THR_G2[G2_IR]);
  out << buf << endl;
}

void printNONEL_COL(ostream& out){
  out << txtKeyWords[NONEL_COL] << endl;
  char buf[1024];
  /*mainchain hydogen related*/
  /*O-N*/
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _O_, min_NO, HB_N_O[2]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _O_HB_, min_NO, HB_N_O[2]);
  out << buf << endl; 
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _O_HB_, min_NO, HB_N_O[2]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _O_, min_NO, HB_N_O[2]);
  out << buf << endl;
  /*N-C*/
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _C_, min_NC, HB_N_C[4]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _C_, min_NC, HB_N_C[4]);
  out << buf << endl;
  /*O-C*/
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_, _C_, min_CO, HB_O_C[3]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_HB_, _C_, min_CO, HB_O_C[3]);
  out << buf << endl;  
  /*O-Ca*/
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_, _CA_, min_OCa, HB_O_CA[3]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_HB_, _CA_, min_OCa, HB_O_CA[3]);
  out << buf << endl;

  /*N-N*/
  sel_col_sprint(buf, _N_, _N_, hdr_NN, min_NN, 2, Dihedral_E/2);
  out << buf << endl;
  sel_col_sprint(buf, _N_, _N_HB_, hdr_NN, min_NN, 2, Dihedral_E/2);
  out << buf << endl;
  sel_col_sprint(buf, _N_HB_, _N_HB_, hdr_NN, min_NN, 2, Dihedral_E/2);
  out << buf << endl;
  sel_col_sprint(buf, _N_, _PRO_N_, hdr_NN, min_PRO_NN, 1, Dihedral_E/2.0);
  out << buf << endl;
  sel_col_sprint(buf, _N_HB_, _PRO_N_, hdr_NN, min_PRO_NN, 1, Dihedral_E/2.0);
  out << buf << endl;
  sel_col_sprint(buf, _PRO_N_, _PRO_N_, hdr_NN, min_PRO_NN, 1, Dihedral_E/2.0);
  out << buf << endl;
  int type_shift = _CYS_CB_;
  /*O-Cb, OP-Cb*/
  sel_col_sprint(buf, _O_, _CB_, hdr_OCb, min_OCb, 2, Dihedral_E/2);
  out << buf << endl;
  sel_col_sprint(buf, _O_HB_, _CB_, hdr_OCb, min_OCb, 2, Dihedral_E/2);
  out << buf << endl;
  for(int i=0; i<20; i++){
    dmd_atom_t cb_t = static_cast<dmd_atom_t>(type_shift + i);
    if(cb_t!=_PRO_CB_){
      sel_col_sprint(buf, _O_, cb_t, hdr_OCb, min_OCb, 2, Dihedral_E/2);
      out << buf << endl;
      sel_col_sprint(buf, _O_HB_, cb_t, hdr_OCb, min_OCb, 2, Dihedral_E/2);
      out << buf << endl;
    }
    else{
      sel_col_sprint(buf, _O_, cb_t, hdr_PRO_OCb, min_PRO_OCb, 2, Dihedral_E/2.0);
      out << buf << endl;
      sel_col_sprint(buf, _O_HB_, cb_t, hdr_PRO_OCb, min_PRO_OCb, 2, Dihedral_E/2.0);
      out << buf << endl;
    }
  }
  /*PRO_CB-N*/
  sel_col_sprint(buf, _N_, _PRO_CB_, hdr_NCb, min_PRO_NCb, 5, Dihedral_E);
  out << buf << endl;
  sel_col_sprint(buf, _N_HB_, _PRO_CB_, hdr_NCb, min_PRO_NCb, 5, Dihedral_E);
  out << buf << endl;
  sel_col_sprint(buf, _PRO_N_, _PRO_CB_, hdr_NCb, min_PRO_NCb, 5, Dihedral_E);
  out << buf << endl;
      
  /*sidechain main related*/
  /*HBA*/
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _ASP_CG_, 
	  N_HDR+G1_CONST[ASP][G1_HDR], HB_ASP_N[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _ASN_CG_, 
	  N_HDR+G1_CONST[ASN][G1_HDR], HB_ASN_N[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _GLN_CG_, 
	  N_HDR+G1_CONST[GLN][G1_HDR], HB_GLN_N[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _GLU_CG_, 
	  N_HDR+G1_CONST[GLU][G1_HDR], HB_GLU_N[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _SER_CG_, min_NO, HB_SER_N[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_, _THR_CG2_, min_NO, HB_THR_N[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _ASP_CG_, 
	  N_HDR+G1_CONST[ASP][G1_HDR], HB_ASP_N[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _ASN_CG_, 
	  N_HDR+G1_CONST[ASN][G1_HDR], HB_ASN_N[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _GLN_CG_, 
	  N_HDR+G1_CONST[GLN][G1_HDR], HB_GLN_N[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _GLU_CG_, 
	  N_HDR+G1_CONST[GLU][G1_HDR], HB_GLU_N[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _SER_CG_, min_NO, HB_SER_N[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _N_HB_, _THR_CG2_, min_NO, HB_THR_N[1]);
  out << buf << endl;
  /*HBD*/
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_, _ASN_CG_, 
	  O_HDR+G1_CONST[ASN][G1_HDR], HBD_ASN_O[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_, _GLN_CG_,
	  O_HDR+G1_CONST[GLN][G1_HDR], HBD_GLN_O[1]);
  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_, _SER_CG_,
	  HBD_SER_O[0], HBD_SER_O[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_, _THR_CG2_,
	  HBD_THR_O[0], HBD_THR_O[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_HB_, _ASN_CG_, 
	  O_HDR+G1_CONST[ASN][G1_HDR], HBD_ASN_O[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_HB_, _GLN_CG_,
	  O_HDR+G1_CONST[GLN][G1_HDR], HBD_GLN_O[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_HB_, _SER_CG_,
	  HBD_SER_O[0], HBD_SER_O[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf 0.0000", _O_HB_, _THR_CG2_,
	  HBD_THR_O[0], HBD_THR_O[1]);
  out << buf << endl;
  
  /*define hydrophobic interaction*/
  dmd_atom_t type_i, type_j;
  double hdr_i, hdr_j;
  double ir_i, ir_j;
  for(int i=0; i<sizeof(dmd_atom_hp)/sizeof(hp_t); i++){
    if(dmd_atom_hp[i]!=IRRELEVANT){
      /*type determination of ith atom*/
      type_i = static_cast<dmd_atom_t>(_CYS_CB_+i);
      if(type_i>=_CYS_CB_ && type_i<=_PRO_CB_){
	hdr_i = CB_HDR;
	ir_i  = CB_IR;
      }
      else if(type_i>=_CYS_CG_ && type_i<=_PRO_CG_){
	hdr_i = G1_CONST[type_i-_CYS_CG_][G1_HDR];
	ir_i  = G1_CONST[type_i-_CYS_CG_][G1_IR];
      }
      else if(type_i==_ILE_CG2_){
	hdr_i = ILE_G2[G2_HDR];
	ir_i  = ILE_G2[G2_IR];
      }
      else if(type_i==_THR_CG2_){
	hdr_i = THR_G2[G2_HDR];
	ir_i  = THR_G2[G2_IR];
      }
      else if(type_i==_VAL_CG2_) {
	hdr_i = VAL_G2[G2_HDR];
	ir_i  = VAL_G2[G2_IR];
      }      
      else if(type_i==_ARG_CD_){
	hdr_i = ARG_D[D_HDR];
	ir_i  = ARG_D[D_IR];
      }
      else if(type_i==_LYS_CD_){
	hdr_i = LYS_D[D_HDR];
	ir_i  = LYS_D[D_IR];
      }
      else if(type_i==_TRP_CD_){
	hdr_i = TRP_D[D_HDR];
	ir_i  = TRP_D[D_IR];
      }
      for(int j=i; j<sizeof(dmd_atom_hp)/sizeof(hp_t); j++){
	/*type determination of ith atom*/
	if(dmd_atom_hp[j]!=IRRELEVANT){
	  type_j = static_cast<dmd_atom_t>(_CYS_CB_+j);
	  if(type_j>=_CYS_CB_ && type_j<=_PRO_CB_){
	    hdr_j = CB_HDR;
	    ir_j  = CB_IR;
	  }
	  else if(type_j>=_CYS_CG_ && type_j<=_PRO_CG_){
	    hdr_j = G1_CONST[type_j-_CYS_CG_][G1_HDR];
	    ir_j  = G1_CONST[type_j-_CYS_CG_][G1_IR];
	  }
	  else if(type_j==_ILE_CG2_){
	    hdr_j = ILE_G2[G2_HDR];
	    ir_j  = ILE_G2[G2_IR];
	  }
	  else if(type_j==_THR_CG2_){
	    hdr_j = THR_G2[G2_HDR];
	    ir_j  = THR_G2[G2_IR];
	  }
	  else if(type_j==_VAL_CG2_){
	    hdr_j = VAL_G2[G2_HDR];
	    ir_j  = VAL_G2[G2_IR];
	  }  
	  else if(type_j==_ARG_CD_){
	    hdr_j = ARG_D[D_HDR];
	    ir_j  = ARG_D[D_IR];
	  }
	  else if(type_j==_LYS_CD_){
	    hdr_j = LYS_D[D_HDR];
	    ir_j  = LYS_D[D_IR];
	  }
	  else if(type_j==_TRP_CD_){
	    hdr_j = TRP_D[D_HDR];
	    ir_j  = TRP_D[D_IR];
	  }
	  /*determine the interaction potential*/
	  if(dmd_atom_hp[i]==APOLOR && dmd_atom_hp[j]==APOLOR){/*APOLOR-APOLOR*/
	    if((type_i==_PHE_CG_ || type_i==_TRP_CD_ || type_i==_TYR_CG_ || type_i==_PRO_CG_) && 
	       (type_j==_PHE_CG_ || type_j==_TRP_CD_ || type_j==_TYR_CG_ || type_j==_PRO_CG_)){
	      /*AROMATIC interaction is more than HYDROPHOBIC;
		The interaction potential will be assigned later
		{PHE,TYR,TRP,PRO}--{PHE,TYR,TRP,PRO}*/
	    }
	    else{
	      hp_sprint(buf, type_i, type_j, hdr_i+hdr_j, ir_i+ir_j,
			HP_E);
	      out << buf << endl;
	    }
	  }
	  else if((dmd_atom_hp[i]==APOLOR && dmd_atom_hp[j]==AMPHI) ||
		  (dmd_atom_hp[i]==AMPHI && dmd_atom_hp[j]==APOLOR)){
	    if((type_i==_PHE_CG_ || type_i==_TRP_CD_ || type_i==_TYR_CG_ || type_i==_PRO_CG_) && 
	       (type_j==_PHE_CG_ || type_j==_TRP_CD_ || type_j==_TYR_CG_ || type_j==_PRO_CG_)){
	      /*PRO as an amphipatic amino acid with aromatic ring, it has a strong
		AROMATIC interaction with hydrophobic aromatic residues;
		The interaction is assgined later
		{PHE,TYR,TRP,PRO}--{PHE,TYR,TRP,PRO}*/
	    }
	    else{
	      hp_sprint(buf, type_i, type_j, hdr_i+hdr_j, ir_i+ir_j,
			HA_E);
	      out << buf << endl;
	    }
	  }
	}
      }
    }
  }
  /*end the definition of HYDROPHOBIC interactions*/
  
  /*SB related:
    We treat the Electrostatic interaction as short-range attraction*/
  double hcd = 0;
  double itd = 0;
  /*ARG(D)-ASP(G)*/
  hcd = ARG_D[D_HDR]+G1_CONST[ASP][G1_HDR];
  itd = ARG_D[D_IR] +G1_CONST[ASP][G1_IR];
  sb_sprint(buf, _ARG_CD_, _ASP_CG_, hcd, itd, SB_E);
  out << buf << endl;
  sb_sprint(buf, _ARG_CD_, _ASP_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  /*ARG(D)-GLU(G)*/
  hcd = ARG_D[D_HDR]+G1_CONST[GLU][G1_HDR];
  itd = ARG_D[D_IR] +G1_CONST[GLU][G1_IR];
  sb_sprint(buf, _ARG_CD_, _GLU_CG_, hcd, itd, SB_E);
  out << buf << endl;
  sb_sprint(buf, _ARG_CD_, _GLU_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  /*ARG(D)-LYS(D)*/
  hcd = ARG_D[D_HDR]+LYS_D[D_HDR];
  itd = ARG_D[D_IR] +LYS_D[D_IR];
  rsb_sprint(buf, _ARG_CD_, _LYS_CD_, hcd, itd, SB_E);
  out << buf << endl;
  /*ARG-ARG*/
  hcd = ARG_D[D_HDR]+ARG_D[D_HDR];
  itd = ARG_D[D_IR] +ARG_D[D_IR];
  rsb_sprint(buf, _ARG_CD_, _ARG_CD_, hcd, itd, SB_E);
  out << buf << endl;  
  /*LYS-ASP*/
  hcd = LYS_D[D_HDR]+G1_CONST[ASP][G1_HDR];
  itd = LYS_D[D_IR] +G1_CONST[ASP][G1_IR];
  sb_sprint(buf, _LYS_CD_, _ASP_CG_, hcd, itd, SB_E);
  out << buf << endl;
  sb_sprint(buf, _LYS_CD_, _ASP_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  /*LYS-GLU*/
  hcd = LYS_D[D_HDR]+G1_CONST[GLU][G1_HDR];
  itd = LYS_D[D_IR] +G1_CONST[GLU][G1_IR];
  sb_sprint(buf, _LYS_CD_, _GLU_CG_, hcd, itd, SB_E);
  out << buf << endl;
  sb_sprint(buf, _LYS_CD_, _GLU_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  /*LYS-LYS*/
  hcd = LYS_D[D_HDR]+LYS_D[D_HDR];
  itd = LYS_D[D_IR] +LYS_D[D_IR];
  rsb_sprint(buf, _LYS_CD_, _LYS_CD_, hcd, itd, SB_E);
  out << buf << endl;
  /*GLU-ASP*/
  hcd = G1_CONST[GLU][G1_HDR]+G1_CONST[ASP][G1_HDR];
  itd = G1_CONST[GLU][G1_IR] +G1_CONST[ASP][G1_IR];
  rsb_sprint(buf, _GLU_CG_,  _ASP_CG_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _GLU_HBA_, _ASP_CG_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _GLU_CG_,  _ASP_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _GLU_HBA_, _ASP_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  /*GLU-GLU*/
  hcd = G1_CONST[GLU][G1_HDR]+G1_CONST[GLU][G1_HDR];
  itd = G1_CONST[GLU][G1_IR] +G1_CONST[GLU][G1_IR];
  rsb_sprint(buf, _GLU_CG_,  _GLU_CG_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _GLU_HBA_, _GLU_CG_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _GLU_CG_,  _GLU_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _GLU_HBA_, _GLU_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  /*ASP-ASP*/
  hcd = G1_CONST[ASP][G1_HDR]+G1_CONST[ASP][G1_HDR];
  itd = G1_CONST[ASP][G1_IR] +G1_CONST[ASP][G1_IR];
  rsb_sprint(buf, _ASP_CG_,  _ASP_CG_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _ASP_HBA_, _ASP_CG_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _ASP_CG_,  _ASP_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  rsb_sprint(buf, _ASP_HBA_, _ASP_HBA_, hcd, itd, SB_E);
  out << buf << endl;
  /*end SB*/
  
  /*define the AROMATIC interactions*/
  /*ARO --- ARO; the hardcore diameter is increased due to ringpacking*/
  /*PHE-PHE*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _PHE_CG_, _PHE_CG_, 
	    ARO_R_PHE+ARO_R_PHE, G1_CONST[PHE][G1_IR] + G1_CONST[PHE][G1_IR],
	    -fabs(AROMATIC_E), 
	    G1_CONST[PHE][G1_IR] + G1_CONST[PHE][G1_IR]+ir_ext,
	  -fabs(AROMATIC_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _PHE_CG_, _PHE_CG_,
	    ARO_R_PHE+ARO_R_PHE, G1_CONST[PHE][G1_IR] + G1_CONST[PHE][G1_IR],
	    -fabs(AROMATIC_E) ) ;
  }
  out << buf << endl;
  /*PHE-TRP*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _PHE_CG_, _TRP_CD_, 
	    ARO_R_PHE+ARO_R_TRP, 
	    G1_CONST[PHE][G1_IR] + TRP_D[D_IR], -fabs(AROMATIC_E),
	    G1_CONST[PHE][G1_IR] + TRP_D[D_IR]+ir_ext, -fabs(AROMATIC_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _PHE_CG_, _TRP_CD_, 
	    ARO_R_PHE+ARO_R_TRP, 
	    G1_CONST[PHE][G1_IR] + TRP_D[D_IR], -fabs(AROMATIC_E) );
  }
  out << buf << endl;
  /*TRP-TRP*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _TRP_CD_, _TRP_CD_, 
	    ARO_R_TRP+ARO_R_TRP, 
	    TRP_D[D_IR] +TRP_D[D_IR], -fabs(AROMATIC_E),
	    TRP_D[D_IR] +TRP_D[D_IR]+ir_ext, -fabs(AROMATIC_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _TRP_CD_, _TRP_CD_, 
	    ARO_R_TRP+ARO_R_TRP, 
	    TRP_D[D_IR] +TRP_D[D_IR], -fabs(AROMATIC_E) );
  }
  out << buf << endl;
  
  /*TRP-TYR*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _TRP_CD_, _TYR_CG_, 
	    ARO_R_TRP + ARO_R_TYR, 
	    TRP_D[D_IR] +G1_CONST[TYR][G1_IR], -fabs(AROMATIC_E),
	    TRP_D[D_IR] +G1_CONST[TYR][G1_IR]+ir_ext, -fabs(AROMATIC_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _TRP_CD_, _TYR_CG_, 
	    ARO_R_TRP + ARO_R_TYR, 
	    TRP_D[D_IR] +G1_CONST[TYR][G1_IR], -fabs(AROMATIC_E) );
  }
  out << buf << endl;
  /*PHE-TYR*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _PHE_CG_, _TYR_CG_, 
	    ARO_R_PHE + ARO_R_TYR, 
	    G1_CONST[PHE][G1_IR] +G1_CONST[TYR][G1_IR], -fabs(AROMATIC_E),
	    G1_CONST[PHE][G1_IR] +G1_CONST[TYR][G1_IR]+ir_ext,
	    -fabs(AROMATIC_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _PHE_CG_, _TYR_CG_, 
	    ARO_R_PHE + ARO_R_TYR, 
	    G1_CONST[PHE][G1_IR] +G1_CONST[TYR][G1_IR], -fabs(AROMATIC_E) );
  }
  out << buf << endl;
  /*TYR-TYR*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _TYR_CG_, _TYR_CG_, 
	    ARO_R_TYR + ARO_R_TYR, 
	    G1_CONST[TYR][G1_IR] +G1_CONST[TYR][G1_IR], -fabs(AROMATIC_E),
	    G1_CONST[TYR][G1_IR] +G1_CONST[TYR][G1_IR]+ir_ext,
	    -fabs(AROMATIC_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _TYR_CG_, _TYR_CG_, 
	    ARO_R_TYR + ARO_R_TYR, 
	    G1_CONST[TYR][G1_IR] +G1_CONST[TYR][G1_IR], -fabs(AROMATIC_E) );
  }
  out << buf << endl;
  
  /*ARO---PRO*/
  /*TRP-PRO*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _TRP_CD_, _PRO_CG_, 
	    TRP_D[D_HDR]+G1_CONST[PRO][G1_HDR], 
	    TRP_D[D_IR] +G1_CONST[PRO][G1_IR], -fabs(ARO_PRO_E),
	    TRP_D[D_IR] +G1_CONST[PRO][G1_IR]+ir_ext, -fabs(ARO_PRO_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _TRP_CD_, _PRO_CG_, 
	    TRP_D[D_HDR]+G1_CONST[PRO][G1_HDR], 
	    TRP_D[D_IR] +G1_CONST[PRO][G1_IR], -fabs(ARO_PRO_E) );
  }
  out << buf << endl;
  /*TYR-PRO*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _TYR_CG_, _PRO_CG_, 
	    G1_CONST[TYR][G1_HDR]+G1_CONST[PRO][G1_HDR], 
	    G1_CONST[TYR][G1_IR] +G1_CONST[PRO][G1_IR], -fabs(ARO_PRO_E),
	    G1_CONST[TYR][G1_IR] +G1_CONST[PRO][G1_IR]+ir_ext,
	    -fabs(ARO_PRO_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _TYR_CG_, _PRO_CG_, 
	    G1_CONST[TYR][G1_HDR]+G1_CONST[PRO][G1_HDR], 
	    G1_CONST[TYR][G1_IR] +G1_CONST[PRO][G1_IR], -fabs(ARO_PRO_E) );
  }
  out << buf << endl;
  /*PHE-PRO*/
  if( ir_ext>0 ){
    sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf", _PHE_CG_, _PRO_CG_, 
	    G1_CONST[PHE][G1_HDR]+G1_CONST[PRO][G1_HDR], 
	    G1_CONST[PHE][G1_IR] +G1_CONST[PRO][G1_IR], -fabs(ARO_PRO_E),
	    G1_CONST[PHE][G1_IR] +G1_CONST[PRO][G1_IR]+ir_ext,
	    -fabs(ARO_PRO_E)/2.0);
  }
  else{
    sprintf(buf, "%ld %ld %lf %lf %lf", _PHE_CG_, _PRO_CG_, 
	    G1_CONST[PHE][G1_HDR]+G1_CONST[PRO][G1_HDR], 
	    G1_CONST[PHE][G1_IR] +G1_CONST[PRO][G1_IR], -fabs(ARO_PRO_E) );
  }
  out << buf << endl;
  
  /*HIS?*/
  /*end the definition of AROMATIC interactions*/
}

void printEL_COL(ostream& out){
  out << txtKeyWords[EL_COL] << endl;
  char buf[1024];
  /*N-C*/
  /*
  el_col_sprint(buf, _N_, _C_, min_NC);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _C_, min_NC);
  out << buf << endl;
  */
  el_col_sprint(buf, _PRO_N_, _C_, min_NC);
  out << buf << endl;
  /*N-Cb*/
  el_col_sprint(buf, _N_, _CB_, min_NCb);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _CB_, min_NCb);
  out << buf << endl;
  el_col_sprint(buf, _PRO_N_, _CB_, min_NCb);
  out << buf << endl;
  int type_shift = _CYS_CB_;
  for(int i=0; i<20; i++){
    dmd_atom_t cb_t = static_cast<dmd_atom_t>(type_shift + i);
    if(cb_t!=_PRO_CB_){
      el_col_sprint(buf, _N_, cb_t, min_NCb);
      out << buf << endl;
      el_col_sprint(buf, _N_HB_, cb_t, min_NCb);
      out << buf << endl;
      el_col_sprint(buf, _PRO_N_, cb_t, min_NCb);
      out << buf << endl;
    }
  }
  /*C-CB*/
  el_col_sprint(buf, _C_, _CB_, min_CCb);
  out << buf << endl;
  for(int i=0; i<20; i++){
    dmd_atom_t cb_t = static_cast<dmd_atom_t>(type_shift + i);
    el_col_sprint(buf, _C_, cb_t, min_CCb);
    out << buf << endl;
  }
  /*O-Ca, OP-Ca*/
  /*
  el_col_sprint(buf, _O_, _CA_, min_OCa);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _CA_, min_OCa);
  out << buf << endl;
  */

  /*main chain hydrogen bonding related*/
  /*NP-OP*/
  /*
  el_col_sprint(buf, _N_HB_, _O_HB_, min_NO);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _O_, min_NO);
  out << buf << endl;
  el_col_sprint(buf, _N_, _O_HB_, min_NO);
  out << buf << endl; 
  */
  /*side-main hydrogen bonding related*/
  el_col_sprint(buf, _N_, _ASP_HBA_, N_HDR+G1_CONST[ASP][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_, _ASN_HBA_, N_HDR+G1_CONST[ASN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_, _GLN_HBA_, N_HDR+G1_CONST[GLN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_, _GLU_HBA_, N_HDR+G1_CONST[GLU][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_, _SER_HBA_, min_NO);
  out << buf << endl;
  el_col_sprint(buf, _N_, _THR_HBA_, min_NO);
  out << buf << endl;  
  el_col_sprint(buf, _N_HB_, _ASP_HBA_, N_HDR+G1_CONST[ASP][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _ASN_HBA_, N_HDR+G1_CONST[ASN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _GLN_HBA_, N_HDR+G1_CONST[GLN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _GLU_HBA_, N_HDR+G1_CONST[GLU][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _SER_HBA_, min_NO);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _THR_HBA_, min_NO);
  out << buf << endl;  
  
  el_col_sprint(buf, _O_, _ASN_HBD_, O_HDR+G1_CONST[ASN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _O_, _GLN_HBD_, O_HDR+G1_CONST[GLN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _O_, _SER_HBD_, HBD_SER_O[0]);
  out << buf << endl;
  el_col_sprint(buf, _O_, _THR_HBD_, HBD_THR_O[0]);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _ASN_HBD_, O_HDR+G1_CONST[ASN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _GLN_HBD_, O_HDR+G1_CONST[GLN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _SER_HBD_, HBD_SER_O[0]);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _THR_HBD_, HBD_THR_O[0]);
  out << buf << endl;
  
  el_col_sprint(buf, _O_, _ASN_HBA_, O_HDR+G1_CONST[ASN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _ASN_HBA_, O_HDR+G1_CONST[ASN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_, _ASN_HBD_, N_HDR+G1_CONST[ASN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _ASN_HBD_, N_HDR+G1_CONST[ASN][G1_HDR]);
  out << buf << endl;
  
  el_col_sprint(buf, _O_, _GLN_HBA_, O_HDR+G1_CONST[GLN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _GLN_HBA_, O_HDR+G1_CONST[GLN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_, _GLN_HBD_, N_HDR+G1_CONST[GLN][G1_HDR]);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _GLN_HBD_, N_HDR+G1_CONST[GLN][G1_HDR]);
  out << buf << endl;
 
  el_col_sprint(buf, _O_, _SER_HBA_, HBD_SER_O[0]);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _SER_HBA_, HBD_SER_O[0]);
  out << buf << endl;
  el_col_sprint(buf, _N_, _SER_HBD_, min_NO);
  out << buf << endl;
  el_col_sprint(buf, _N_HB_, _SER_HBD_, min_NO);
  out << buf << endl;
  
  el_col_sprint(buf, _O_, _THR_HBA_, HBD_THR_O[0]);
  out << buf << endl;
  el_col_sprint(buf, _O_HB_, _THR_HBA_, HBD_THR_O[0]);
  out << buf << endl;
  el_col_sprint(buf, _N_, _THR_HBD_, min_NO);
  out << buf << endl;  
  el_col_sprint(buf, _N_HB_, _THR_HBD_, min_NO);
  out << buf << endl;  
  
}

void printBOND_TYPE(ostream& out){
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

  /*CA-CB*/
  single_bond_sprint(buf, CA_CB, BOND_DEV);
  out << _CA_ << " " << _CB_ << " " << buf << endl;
  /*CA-PRO_CB*/  
  out << _CA_ << " " << _PRO_CB_ << " " << buf << endl;
  
  /*N-CB*/
  n_dev = getTriDev_angle(N_CA, CA_CB, BOND_DEV, N_CA_CB*rpi, n_len);
  single_bond_sprint(buf, n_len, n_dev);
  out << _N_ << " " << _CB_ << " " << buf << endl;

  /*N-PRO_CB*/
  n_dev = getTriDev_angle(N_CA, CA_CB, BOND_DEV, PRO_N_CA_CB*rpi, n_len);
  single_bond_sprint(buf, n_len, n_dev);
  out << _N_ << " " << _PRO_CB_ << " " << buf << endl;

  /*C-CB*/
  n_dev = getTriDev_angle(CA_C, CA_CB, BOND_DEV, C_CA_CB*rpi, n_len);
  single_bond_sprint(buf, n_len, n_dev);
  out << _C_ << " " << _CB_ << " " << buf << endl;
  /*C-PRO_CB, 
    C-PRO_CB1, used to fix the phi angle of PRO*/
  n_dev = getTriDev_angle(CA_C, CA_CB, BOND_DEV, PRO_C_CA_CB*rpi, n_len);
  double_bond_sprint(buf, n_len, n_dev, PRO_C_CB, PRO_C_CB_D/PRO_C_CB);
  out << _C_ << " " << _PRO_CB_ << " " << buf << endl;

  /*C-C, used to alow certain phi region for PRO*/
  sprintf(buf, "%ld %ld %lf %lf", _C_, _C_, PRO_C_C1_MIN, PRO_C_C1_MAX);
  out << buf << endl;

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

  /*CG-CB, CG-CA and C-CG, N-CG*/
  int type_shift = _CYS_CG_;
  for(int i=0; i<20; i++){
    dmd_atom_t cg_t = static_cast<dmd_atom_t>(type_shift + i);
    if(i!=ALA && i!=GLY){
      const double* cg_para = G1_CONST[i];
      if(cg_para[G1_CB_D]==0){
	single_bond_sprint(buf, cg_para[G1_CB], BOND_DEV_G);
      }
      else{
	single_bond_sprint_len(buf, cg_para[G1_CB], cg_para[G1_CB_D]);
      }
      if(i!=PRO)
	out << _CB_ << " " << cg_t << " " << buf << endl;
      else
	out << _PRO_CB_ << " " << cg_t << " " << buf << endl;

      if(cg_para[G1_CA_D]==0){
	n_dev = getTriDev_length(cg_para[G1_CB], CA_CB, BOND_DEV_G, cg_para[G1_CA]);
	single_bond_sprint(buf, cg_para[G1_CA], n_dev);
      }
      else{
	single_bond_sprint_len(buf, cg_para[G1_CA], cg_para[G1_CA_D]);
      }
      out << _CA_ << " " << cg_t << " " << buf << endl;

      /*N-CG, C-CG*/
      const double* dihedral = Dihedral_C_CG[i];
      if(dihedral[0]!=INF){
	//c-cg
	sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf %lf %lf",
		_C_, cg_t, 
		Dihedral_MIN, 
		dihedral[0], Dihedral_E,
		dihedral[1], -Dihedral_E, 
		dihedral[2], Dihedral_E, 
		Dihedral_MAX);
	out << buf << endl;
	//n-cg
	sprintf(buf, "%ld %ld %lf %lf",
		_N_, cg_t, Dihedral_MIN, Dihedral_MAX);
	out << buf << endl;
      }
      
      if(i==PRO){//different constraint for PROLINE
	single_bond_sprint_len(buf, PRO_N_CG, PRO_N_CG_D);
	out << _N_ << " " << cg_t << " " << buf << endl;
	
	/*this constraint is actually used to
	  avoid the small distance of c-PRO_CB*/
	single_bond_sprint_len(buf, PRO_C_CG, PRO_C_CG_D);
	out << _C_ << " " << cg_t << " " << buf << endl;
      }
      
    }
  }

  /*ILE, G2*/
  //CG2-CB
  single_bond_sprint(buf, ILE_G2[G2_CB], BOND_DEV_G);
  out << _CB_ << " " << _ILE_CG2_ << " " << buf << endl;
  //CG2-CA
  n_dev = getTriDev_length(CA_CB, ILE_G2[G2_CB], BOND_DEV_G, ILE_G2[G2_CA]);
  single_bond_sprint(buf, ILE_G2[G2_CA], n_dev);
  out << _CA_ << " " << _ILE_CG2_ << " " << buf << endl;
  //CG2-CG1
  //n_dev = getTriDev_length(ILE_G2[G2_CB], G1_CONST[ILE][G1_CB], 
  //BOND_DEV_G, ILE_G2[G1_G2]);
  n_dev = ILE_G2[CG_G2_D]/ILE_G2[G1_G2];
  single_bond_sprint(buf, ILE_G2[G1_G2], n_dev);
  out << _ILE_CG_ << " " << _ILE_CG2_ << " " << buf << endl;

  //N-CG2
  sprintf(buf, "%ld %ld %lf %lf",
	  _N_, _ILE_CG2_, Dihedral_MIN, Dihedral_MAX);
  out << buf << endl;
  //C-CG2
  const double* dihedral = Dihedral_C_CG[ILE];
  sprintf(buf, "%ld %ld %lf %lf",
	  _C_, _ILE_CG2_, Dihedral_MIN, Dihedral_MAX);
  out << buf << endl;
  
  /*THR, G2*/
  //CG2-CB
  single_bond_sprint(buf, THR_G2[G2_CB], BOND_DEV_G);
  out << _CB_ << " " << _THR_CG2_ << " " << buf << endl;
  //CG2-CA
  n_dev = getTriDev_length(CA_CB, THR_G2[G2_CB], BOND_DEV_G, THR_G2[G2_CA]);
  single_bond_sprint(buf, THR_G2[G2_CA], n_dev);
  out << _CA_ << " " << _THR_CG2_ << " " << buf << endl;
  //CG2-CG1
  n_dev = getTriDev_length(THR_G2[G2_CB], G1_CONST[THR][G1_CB], 
			   BOND_DEV_G, THR_G2[G1_G2]);
  single_bond_sprint(buf, THR_G2[G1_G2], n_dev);
  out << _THR_CG_ << " " << _THR_CG2_ << " " << buf << endl;
 
  //N-CG2
  sprintf(buf, "%ld %ld %lf %lf",
	  _N_, _THR_CG2_, Dihedral_MIN, Dihedral_MAX);
  out << buf << endl;
  //C-CG2
  sprintf(buf, "%ld %ld %lf %lf",
	  _C_, _THR_CG2_, Dihedral_MIN, Dihedral_MAX);
  out << buf << endl;

  /*VAL, G2*/
  //CG2-CB
  single_bond_sprint(buf, VAL_G2[G2_CB], BOND_DEV_G);
  out << _CB_ << " " << _VAL_CG2_ << " " << buf << endl;
  //CG2-CA
  n_dev = getTriDev_length(CA_CB, VAL_G2[G2_CB], BOND_DEV_G, VAL_G2[G2_CA]);
  single_bond_sprint(buf, VAL_G2[G2_CA], n_dev);
  out << _CA_ << " " << _VAL_CG2_ << " " << buf << endl;
  //CG2-CG1
  n_dev = getTriDev_length(VAL_G2[G2_CB], G1_CONST[VAL][G1_CB], 
			   BOND_DEV_G, VAL_G2[G1_G2]);
  single_bond_sprint(buf, VAL_G2[G1_G2], n_dev);
  out << _VAL_CG_ << " " << _VAL_CG2_ << " " << buf << endl;
  
  //N-CG2
  sprintf(buf, "%ld %ld %lf %lf",
	  _N_, _VAL_CG2_, Dihedral_MIN, Dihedral_MAX);
  out << buf << endl;
  //C-CG2
  sprintf(buf, "%ld %ld %lf %lf",
	  _C_, _VAL_CG2_, Dihedral_MIN, Dihedral_MAX);
  out << buf << endl;

  /*TRP, CD*/
  /*CD-CG*/
  single_bond_sprint(buf, TRP_D[D_CG], BOND_DEV_G);
  out << _TRP_CG_ << " " << _TRP_CD_ << " " << buf << endl;
  /*CD-CB*/
  n_dev = getTriDev_length(G1_CONST[TRP][G1_CB], TRP_D[D_CG], BOND_DEV_G, TRP_D[D_CB]);
  single_bond_sprint(buf, TRP_D[D_CB], n_dev);
  out << _TRP_CD_ << " " << _CB_ << " " << buf << endl;
  /*CD-CA; Dihedral*/
  sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf %lf %lf",
	  _CA_, _TRP_CD_, 
	  Dihedral_MIN, 
	  Dihedral_TRP_D[0], Dihedral_E,
	  Dihedral_TRP_D[1], -Dihedral_E, 
	  Dihedral_TRP_D[2], Dihedral_E, 
	  Dihedral_MAX);
  out << buf << endl;

  /*ARG, CD*/
  /*CD-CG*/
  single_bond_sprint(buf, ARG_D[D_CG], BOND_DEV_G);
  out << _ARG_CG_ << " " << _ARG_CD_ << " " << buf << endl;
  /*CD-CB*/
  n_dev = ARG_D[D_CB_D]/ARG_D[D_CB];
  single_bond_sprint(buf, ARG_D[D_CB], n_dev);
  out << _ARG_CD_ << " " << _CB_ << " " << buf << endl;

  /*LYS, CD*/
  /*CD-CG*/
  single_bond_sprint(buf, LYS_D[D_CG], BOND_DEV_G);
  out << _LYS_CG_ << " " << _LYS_CD_ << " " << buf << endl;
  /*CD-CB*/
  n_dev = LYS_D[D_CB_D]/LYS_D[D_CB];
  single_bond_sprint(buf, LYS_D[D_CB], n_dev);
  out << _LYS_CD_ << " " << _CB_ << " " << buf << endl;
  
  /*mainchain hydrogen bond related*/
  /*NP-OP*/
  sprintf(buf,"%ld %ld %lf %lf %lf %lf 0.00000", _N_HB_, _O_HB_, 
	  HB_N_O[0], HB_N_O[1], -fabs(EHB_M_M)/4.0, HB_N_O[2]+eps);
  out << buf << endl;
  /*NP-C*/
  sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf %lf %lf 0.00000", 
  	  _N_HB_, _C_,
  	  HB_N_C[0], 
  	  HB_N_C[1],  fabs(EHB_M_M)/4.0, 
	  HB_N_C[2],  fabs(EHB_M_M)/4.0,
	  HB_N_C[3], -fabs(EHB_M_M)/4.0,
  	  HB_N_C[4]+eps);
  out << buf << endl;
  //sprintf(buf, "%ld %ld %lf %lf", 
  //_N_HB_, _C_,
  //HB_N_C[1],
  //HB_N_C[2]), 
  //out << buf << endl;
  /*OP-C*/
  sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf 0.00000", 
  	  _O_HB_, _C_,
  	  HB_O_C[0], 
  	  HB_O_C[1],  fabs(EHB_M_M)/2.0, 
  	  HB_O_C[2], -fabs(EHB_M_M)/4.0, 
  	  HB_O_C[3]+eps); 
  out << buf << endl;
  //sprintf(buf, "%ld %ld %lf %lf", 
  //	  _O_HB_, _C_,
  //	  HB_O_C[1],
  //	  HB_O_C[2]);
  //out << buf << endl;
  /*OP-CA*/
  sprintf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf 0.0000", 
  	  _O_HB_, _CA_,
  	  HB_O_CA[0], 
  	  HB_O_CA[1],  fabs(EHB_M_M)/2.0, 
  	  HB_O_CA[2], -fabs(EHB_M_M)/4.0, 
  	  HB_O_CA[3]+eps);
  out << buf << endl;
  //sprintf(buf, "%ld %ld %lf %lf", 
  //	  _O_HB_, _CA_,
  //	  HB_O_CA[1],
  //	  HB_O_CA[2]);
  //out << buf << endl;


  /*side-main hydrogen bond related*/
  //HBA
  /*ASP*/
  sprintf(buf, "%ld %ld %lf %lf %lf", _N_HB_, _ASP_HBA_, 
	  HB_ASP_N[0], HB_ASP_N[1]+eps, -fabs(EHB_S_M));
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _ASP_HBA_, _C_,
	  HB_ASP_C[0], HB_ASP_C[1]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _ASP_HBA_, _CA_,
	  HB_ASP_CA[0], HB_ASP_CA[1]);
  out << buf << endl;
  /*ASN*/
  sprintf(buf, "%ld %ld %lf %lf %lf", _N_HB_, _ASN_HBA_, 
	  HB_ASN_N[0], HB_ASN_N[1]+eps, -fabs(EHB_S_M));
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _ASN_HBA_, _C_,
	  HB_ASN_C[0], HB_ASN_C[1]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _ASN_HBA_, _CA_,
	  HB_ASN_CA[0], HB_ASN_CA[1]);
  out << buf << endl;
//  /*GLU*/
//  sprintf(buf, "%ld %ld %lf %lf %lf", _N_HB_, _GLU_HBA_, 
//	  HB_GLU_N[0], HB_GLU_N[1]+eps, -fabs(EHB_S_M));
//  out << buf << endl;
//  sprintf(buf, "%ld %ld %lf %lf", _GLU_HBA_, _C_,
//	  HB_GLU_C[0], HB_GLU_C[1]);
//  out << buf << endl;
//  sprintf(buf, "%ld %ld %lf %lf", _GLU_HBA_, _CA_,
//	  HB_GLU_CA[0], HB_GLU_CA[1]);
//  out << buf << endl;
//  /*GLN*/
//  sprintf(buf, "%ld %ld %lf %lf %lf", _N_HB_, _GLN_HBA_, 
//	  HB_GLN_N[0], HB_GLN_N[1]+eps, -fabs(EHB_S_M));
//  out << buf << endl;
//  sprintf(buf, "%ld %ld %lf %lf", _GLN_HBA_, _C_,
//	  HB_GLN_C[0], HB_GLN_C[1]);
//  out << buf << endl;
//  sprintf(buf, "%ld %ld %lf %lf", _GLN_HBA_, _CA_,
//	  HB_GLN_CA[0], HB_GLN_CA[1]);
//  out << buf << endl;
  /*SER*/
  sprintf(buf, "%ld %ld %lf %lf %lf", _N_HB_, _SER_HBA_, 
	  HB_SER_N[0], HB_SER_N[1]+eps, -fabs(EHB_S_M));
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _SER_HBA_, _C_,
	  HB_SER_C[0], HB_SER_C[1]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _SER_HBA_, _CA_,
	  HB_SER_CA[0], HB_SER_CA[1]);
  out << buf << endl;
  /*THR*/
  sprintf(buf, "%ld %ld %lf %lf %lf", _N_HB_, _THR_HBA_, 
	  HB_THR_N[0], HB_THR_N[1]+eps, -fabs(EHB_S_M));
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _THR_HBA_, _C_,
	  HB_THR_C[0], HB_THR_C[1]);
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _THR_HBA_, _CA_,
	  HB_THR_CA[0], HB_THR_CA[1]);
  out << buf << endl;
  //HBD
//  sprintf(buf, "%ld %ld %lf %lf %lf", _O_HB_, _ASN_HBD_,
//	  HBD_ASN_O[0], HBD_ASN_O[1]+eps, -fabs(EHB_S_M));
//  out << buf << endl;
//  sprintf(buf, "%ld %ld %lf %lf", _C_, _ASN_HBD_,
//	  HBD_ASN_C[0], HBD_ASN_C[1]);
//  out << buf << endl;

//  sprintf(buf, "%ld %ld %lf %lf %lf", _O_HB_, _GLN_HBD_,
//	  HBD_GLN_O[0], HBD_GLN_O[1]+eps, -fabs(EHB_S_M));
//  out << buf << endl;
//  sprintf(buf, "%ld %ld %lf %lf", _C_, _GLN_HBD_,
//	  HBD_GLN_C[0], HBD_GLN_C[1]);
//  out << buf << endl;

  sprintf(buf, "%ld %ld %lf %lf %lf", _O_HB_, _SER_HBD_,
	  HBD_SER_O[0], HBD_SER_O[1]+eps, -fabs(EHB_S_M));
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _C_, _SER_HBD_,
	  HBD_SER_C[0], HBD_SER_C[1]);
  out << buf << endl;
  
  sprintf(buf, "%ld %ld %lf %lf %lf", _O_HB_, _THR_HBD_,
	  HBD_THR_O[0], HBD_THR_O[1]+eps, -fabs(EHB_S_M));
  out << buf << endl;
  sprintf(buf, "%ld %ld %lf %lf", _C_, _THR_HBD_,
	  HBD_THR_C[0], HBD_THR_C[1]);
  out << buf << endl;

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
  /*side-main hydrogen bonding*/
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _N_, _ASP_CG_, _N_HB_, _ASP_HBA_);
  out << buf << endl;
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _N_, _ASN_CG_, _N_HB_, _ASN_HBA_);
  out << buf << endl;
//  sprintf(buf, "%ld %ld %ld %ld 1",
//	  _N_, _GLU_CG_, _N_HB_, _GLU_HBA_);
//  out << buf << endl;
//  sprintf(buf, "%ld %ld %ld %ld 1",
//	  _N_, _GLN_CG_, _N_HB_, _GLN_HBA_);
//  out << buf << endl;
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _N_, _SER_CG_, _N_HB_, _SER_HBA_);
  out << buf << endl;
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _N_, _THR_CG2_, _N_HB_, _THR_HBA_);
  out << buf << endl;

//  sprintf(buf, "%ld %ld %ld %ld 1",
//	  _C_, _ASN_CG_, _O_HB_, _ASN_HBD_);
//  out << buf << endl;
//  sprintf(buf, "%ld %ld %ld %ld 1",
//	  _C_, _GLN_CG_, _O_HB_, _GLN_HBD_);
//  out << buf << endl;
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _O_, _SER_CG_, _O_HB_, _SER_HBD_);
  out << buf << endl;
  sprintf(buf, "%ld %ld %ld %ld 1",
	  _O_, _THR_CG_, _O_HB_, _THR_HBD_);
  out << buf << endl;
}

void printATOM_LIST(ostream& out, pseudoPep& p, vector<double>xyz, randomGenerator& r){
  out << txtKeyWords[LIST_ATOMS] << endl;
  char buf[1024];
  vec v;
  vec pos;
  int index=1;
  atom* CB;
  atom* CG;
  atom* CG2;
  atom* CD;
  dmd_atom_t type;
  for(int i=0; i<p.getLength(); i++){
    pseudoAA* theAA = p.getResidue(i);
    /*N*/
    v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
    pos = vec(xyz[3*index-3],xyz[3*index-2],xyz[3*index-1]);
    if(theAA->getID()==PRO){
      list_atom_sprint(buf, index++, _PRO_N_, pos, v);
    }
    else{
      
      if(i==0)/*first N atom will be hydrogen-active!*/
	list_atom_sprint(buf, index++, _N_HB_, pos, v);
      else
	list_atom_sprint(buf, index++, _N_, pos, v);
    }
    out << buf << endl;
      

    /*C*/
    v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
    pos = vec(xyz[3*index-3],xyz[3*index-2],xyz[3*index-1]);
    list_atom_sprint(buf, index++, _C_, pos, v);
    out << buf << endl;
    /*O*/
    v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
    pos = vec(xyz[3*index-3],xyz[3*index-2],xyz[3*index-1]);
    if(i<p.getLength()-1){
      list_atom_sprint(buf, index++, _O_, pos, v);
    }
    else{
      list_atom_sprint(buf, index++, _O_HB_, pos, v);
    }
    out << buf << endl;
    /*CA*/
    pos = vec(xyz[3*index-3],xyz[3*index-2],xyz[3*index-1]);
    v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
    list_atom_sprint(buf, index++, _CA_, pos, v);
    out << buf << endl;
    /*CB*/
    if(CB = theAA->getCB()){
      type = static_cast<dmd_atom_t>(theAA->getID()+_CYS_CB_);
      v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
      pos = vec(xyz[3*index-3],xyz[3*index-2],xyz[3*index-1]);
      list_atom_sprint(buf, index++, type, pos, v);
      out << buf << endl;
      /*CG*/
      if(CG = theAA->getG()){
	type = static_cast<dmd_atom_t>(theAA->getID()+_CYS_CG_);
	double mass = G1_CONST[theAA->getID()][G1_M];
	v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss())/sqrt(mass);
	pos = vec(xyz[3*index-3],xyz[3*index-2],xyz[3*index-1]);
	list_atom_sprint(buf, index++, type, pos, v);
	out << buf << endl;
	/*CG2*/
	if(CG2 = theAA->getG2()){
	  if(theAA->getID()==ILE) type = _ILE_CG2_;
	  else if(theAA->getID()==THR) type = _THR_CG2_;
	  else if(theAA->getID()==VAL) type = _VAL_CG2_;
	  else{
	    cerr << "error in type" << endl;
	    exit(1);
	  }
	  v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
	  pos = vec(xyz[3*index-3],xyz[3*index-2],xyz[3*index-1]);
	  list_atom_sprint(buf, index++, type, pos, v);
	  out << buf << endl;
	}
	/*CD*/
	if(CD = theAA->getD()){
	  if(theAA->getID()==TRP) type = _TRP_CD_;
	  else if(theAA->getID()==ARG) type = _ARG_CD_;
	  else if(theAA->getID()==LYS) type = _LYS_CD_;
	  else{
	    cerr << "error in type" << endl;
	    exit(1);
	  }
	  v = vec(r.nextGauss(), r.nextGauss(), r.nextGauss());
	  pos = vec(xyz[3*index-3],xyz[3*index-2],xyz[3*index-1]);
	  list_atom_sprint(buf, index++, type, pos, v);
	  out << buf << endl;
	}
      }
    }
  }
}

void printBOND_LIST(ostream& out, pseudoPep& p){
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
	if(theAA->getD()){
	  dmd_atom_t dt;
	  if(theAA->getID()==ARG) dt = _ARG_CD_;
	  else if(theAA->getID()==LYS) dt = _LYS_CD_;
	  else if(theAA->getID()==TRP) dt = _TRP_CD_;
	  /*CG-CD*/
	  out << theAA->getG()->getIndex() << " " << theAA->getD()->getIndex() << endl;
	  /*CB-CD*/
	  out << theAA->getCB()->getIndex() << " " << theAA->getD()->getIndex() << endl;
	  /*CA-CD---Dihedral of TRP*/
	  if(theAA->getID()==TRP){
	    out << theAA->getCA()->getIndex() << " " << theAA->getD()->getIndex() << endl;
	  }
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
}

void printPERM_BOND_LIST(ostream& out, pseudoPep& p){
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
      dmd_atom_t  cbt = _CB_;
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
	dmd_atom_t gt = static_cast<dmd_atom_t>(_CYS_CG_ + theAA->getID());
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
	  dmd_atom_t g2t;
	  if(theAA->getID() == ILE) g2t = _ILE_CG2_;
	  else if(theAA->getID() == THR) g2t = _THR_CG2_;
	  else if(theAA->getID() == VAL) g2t = _VAL_CG2_;
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
	if(theAA->getD()){/*bulky residues with Cd*/
	  dmd_atom_t dt;
	  if(theAA->getID()==ARG) dt = _ARG_CD_;
	  else if(theAA->getID()==LYS) dt = _LYS_CD_;
	  else if(theAA->getID()==TRP) dt = _TRP_CD_;
	  /*CG-CD*/
	  out << theAA->getG()->getIndex() << " " << theAA->getD()->getIndex() << " "
	      << gt << " " << dt << endl;
	  /*CB-CD*/
	  out << theAA->getCB()->getIndex() << " " << theAA->getD()->getIndex() << " "
	      << cbt << " " << dt << endl;
	  /*CA-CD, TRP--Dihedral*/
	  if(theAA->getID()==TRP){
	    out << theAA->getCA()->getIndex() << " " << theAA->getD()->getIndex() << " " 
		<< _CA_ << " " << _TRP_CD_ << endl;
	  }
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
	out << preAA->getC()->getIndex() << " " << theAA->getCB()->getIndex() << " "
	    << _C_ << " " << _PRO_CB_ << endl;
	out << preAA->getC()->getIndex() << " " << theAA->getG()->getIndex() << " "
	    << _C_ << " " << _PRO_CG_ << endl;
	out << preAA->getC()->getIndex() << " " << theAA->getC()->getIndex() << " "
	    << _C_ << " " << _C_ << endl;
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
    /*HBA*/
    aa_t type = static_cast<aa_t>(p.getResidue(i)->getID());
    /*do not consider GLN/GLU*/
    if(type==ASP || type==ASN || type==SER){
      sprintf(buf, "%ld %ld", p.getResidue(i)->getG()->getIndex(), i);
      out << buf << endl;
    }
    else if(type==THR){
      sprintf(buf, "%ld %ld", p.getResidue(i)->getG2()->getIndex(), i);
      out << buf << endl;
    }
  }
}

void print_GID(ostream& out, pseudoPep& p){
  out << txtKeyWords[GID_LIST] << endl;
  char buf[1024];
  for(int i=0; i<p.getLength(); i++){
    aa_t type = static_cast<aa_t>(p.getResidue(i)->getID());
    if(type!=GLY){
      /*CB*/
      sprintf(buf, "%ld %ld", p.getResidue(i)->getCB()->getIndex(), i);
      out << buf<< endl;
      /*CG*/
      if(type!=ALA){
	sprintf(buf, "%ld %ld", p.getResidue(i)->getG()->getIndex(), i);
	out << buf << endl;
	/*CG2*/
	if(type==ILE || type==VAL || type==THR){
	  sprintf(buf, "%ld %ld", p.getResidue(i)->getG2()->getIndex(), i);
	  out << buf << endl;
	}
	/*CD*/
	if(type==ARG || type==LYS || type==TRP){
	  sprintf(buf, "%ld %ld", p.getResidue(i)->getD()->getIndex(), i);
	  out << buf << endl;
	}
      }
    }
  }
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
	if(theAA->getD()){
	  tmp += *theAA->getD()->getR();
	  natom++;
	}
      }
    }
  }
  tmp /= static_cast<double>(natom);
  return tmp;
}

typedef enum {
  NON_ENER_KEY=-1,E_HB_MM=0, E_HB_SM, E_HH, E_HA, E_SB, E_AROMATIC, E_ARO_PRO, E_ROTAMER
} ener_t;

const static char* ener_key[] = {
  "E_HB_MM",
  "E_HB_SM",
  "E_HH",
  "E_HA",
  "E_SB",
  "E_AROMATIC",
  "E_ARO_PRO",
  "E_ROTAMER"
};

void process_para(const char* para_f){
  FILE* in = fopen(para_f,"r");
  char keyword[100];
  char val[100];
  ener_t key;
  if(!in) cerr << "Can not open the parameter file!\n" << endl;
  while(!feof(in)){
    fscanf(in, "%s%s", keyword, val);
    key = NON_ENER_KEY;
    for(int i=0; i<sizeof(ener_key)/sizeof(char*); i++){
      if(strcmp(ener_key[i], keyword)==0){
	key = static_cast<ener_t>(i);
	break;
      }
    }
    switch(key){
    case E_HB_MM:
      setMMHB_E(atof(val));
      break;
    case E_HB_SM:
      setSMHB_E(atof(val));
      break;
    case E_HH:
      setHP_E(atof(val));
      break;
    case E_HA:
      setHA_E(atof(val));
      break;
    case E_SB:
      setSB_E(atof(val));
      break;
    case E_AROMATIC:
      setAROMATIC_E(atof(val));
      break;
    case E_ARO_PRO:
      setARO_PRO_E(atof(val));
      break;
    case E_ROTAMER:
      setDI_E(atof(val));
      break;
    default:
      break;
    }
  }
  fclose(in);
}

int main(int argc, char* argv[]){
  if(argc<5){
    cout << "usage: command seq relaxedTXT E_Para_F  ir_ext [is_ExN]" << endl;
    cout << "ir_ext: external radius. If set to zero, then we do not add external radius";
    cout << " is_ExN -- OPT parameter to exclude neightbouring interaction: 0/1" << endl;
    exit(1);
  }
  
  process_para(argv[3]);
  ir_ext = atof( argv[4] );
  int is_en=0;
  if(argc==6){
    is_en = atoi(argv[5]);
  }

  pseudoPep p(argv[1], (input_t)0);
  pseudoAA* lastAA=p.getResidue(p.getLength()-1);
  lastAA->getO()->getR()->Rotate(*lastAA->getCA()->getR(),*lastAA->getC()->getR(), -PI/3.0);
  p.makeIndex();
  
  randomGenerator ran(_RAN2_, -100);
  
  double size;
  vec cm = getCM(p);
  vec shift = vec(size/2.0, size/2.0, size/2.0) - cm;
  p.shift(shift);
  
  vector<double>xyz;
  process_txt(argv[2], size, xyz);
  
  /*create the TXT*/
  printSYS_SIZE(cout, size);
  printNUM_ATOMS(cout, p);
  printATOM_TYPE(cout);
  printNONEL_COL(cout);
  printEL_COL(cout);
  printBOND_TYPE(cout);
  printREACT(cout);
  
  printATOM_LIST(cout, p, xyz, ran);
  
  printBOND_LIST(cout, p);
  printPERM_BOND_LIST(cout, p);
  printHBA_LIST(cout, p);
  if(is_en) print_GID(cout, p);
  /*END of the txt output*/
}
