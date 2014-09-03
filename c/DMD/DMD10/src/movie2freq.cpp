#include "pseudoPep.h"
#include "pseudoAA.h"
#include "movie.h"
#include "pseudoTXT.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

typedef enum {
  NON_ENER_KEY=-1,E_HB_MM=0, E_HB_SM, E_HP, E_SB, E_AROMATIC, E_ROTAMER
} ener_t;

const static char* ener_key[] = {
  "E_HB_MM",
  "E_HB_SM",
  "E_HP",
  "E_SB",
  "E_AROMATIC",
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
    case E_HP:
      setHP_E(atof(val));
      break;
    case E_SB:
      setSB_E(atof(val));
      break;
    case E_AROMATIC:
      setAROMATIC_E(atof(val));
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

void pepSetAtomType(pseudoPep& p){
  for(int i=0; i<p.getLength(); i++){
    pseudoAA* ia = p.getResidue(i);
    ia->getN()->setType(_N_);
    if(ia->getH())ia->getH()->setType(_NON_ATOM_);
    ia->getC()->setType(_C_);
    ia->getO()->setType(_O_);
    ia->getCA()->setType(_CA_);
    if(ia->getCB())      ia->getCB()->setType(_CYS_CB_+ia->getID());
    if(ia->getG())       ia->getG()->setType( _CYS_CG_+ia->getID());
    if(ia->getID()==THR) ia->getG2()->setType(_THR_CG2_);
    if(ia->getID()==ILE) ia->getG2()->setType(_ILE_CG2_);
    if(ia->getID()==VAL) ia->getG2()->setType(_VAL_CG2_);
  }
}


double* createIR(){
  double* table = new double[_VAL_CG2_ - _CYS_CB_ +1];
  for(int i=_CYS_CB_; i<=_PRO_CB_; i++){
    if(dmd_atom_hp[i-_CYS_CB_]==APOLOR ||
       dmd_atom_hp[i-_CYS_CB_]==AMPHI){
      table[i-_CYS_CB_]=CB_IR;
    }
    else{
      table[i-_CYS_CB_]=0;
    }
  }
  for(int i=_CYS_CG_; i<=_PRO_CG_; i++){
    table[i-_CYS_CB_] = G1_CONST[i-_CYS_CG_][G1_IR];
  }
  table[_ILE_CG2_-_CYS_CB_]=ILE_G2[G2_IR];
  table[_THR_CG2_-_CYS_CB_]=THR_G2[G2_IR];
  table[_VAL_CG2_-_CYS_CB_]=VAL_G2[G2_IR];
  return table;
}

static int n1=0;
static int n2=0;
static int n1r=0;
static int n2r=0;

void pep2CM( pseudoPep& p, const double* IR_table, short** bonds, double** cm_ss, double** cm_mm){


  int natom = p.getLength();
  for(int i=0; i<natom; i++){
    pseudoAA* ia = p.getResidue(i);
    int ia_nside = ia->getNSideChains();
    atom* ia_side= ia->getSideChains();
    for(int j=0; j<i-1; j++){
      pseudoAA* ja = p.getResidue(j);
      int ja_nside = ja->getNSideChains();
      atom* ja_side= ja->getSideChains();
      int isContact = 0;
      for(int is = 0; is<ia_nside; is++){
	for(int js = 0; js<ja_nside; js++){
	  if(ia_side[is].getDist(ja_side[js]) <=
	     (IR_table[ia_side[is].getType()-_CYS_CB_]+
	      IR_table[ja_side[js].getType()-_CYS_CB_])){
	    isContact=1;
	  }
	}
      }
      cm_ss[i][j]+=(double)isContact;
      //mainchain hydrogen bonds
      if(j<i-2){
	//int isHB=0;
	int iN = ia->getN()->getIndex();
	int jO = ja->getO()->getIndex();
	int found =0;
	int ii=0;
	while(bonds[jO-1][ii]>=0){
	  if(bonds[jO-1][ii]==iN-1){
	    found=1;
	    break;
	  }
	  ii++;
	}
	int iO = ia->getO()->getIndex();
	int jN = ja->getN()->getIndex();
	ii=0;
	while(bonds[jN-1][ii]>=0){
	  if(bonds[jN-1][ii]==iO-1){
	    found=1;
	    break;
	  }
	  ii++;
	}
	
	/*
	if(ia->getID()!=PRO && ia->getH()){
	  if(ia->getN()->getDist(*ja->getO())<3.2){
	    vec NH = *(ia->getH()->getR())-*(ia->getN()->getR());
	    vec HO = *(ja->getO()->getR())-*(ia->getH()->getR());
	    vec CO = *(ja->getO()->getR())-*(ja->getC()->getR());
	    NH.Normalize();
	    HO.Normalize();
	    CO.Normalize();
	 
	    if( NH*HO>0.5 && CO*HO<-0.5){
	      isHB=1;
	    }
	  }
	}
	if(ja->getID()!=PRO && ja->getH()){
	  if(ja->getN()->getDist(*ia->getO())<3.2){
	    vec NH = *(ja->getH()->getR())-*(ja->getN()->getR());
	    vec HO = *(ia->getO()->getR())-*(ja->getH()->getR());
	    vec CO = *(ia->getO()->getR())-*(ia->getC()->getR());
	    NH.Normalize();
            HO.Normalize();
            CO.Normalize();
	    
	    if( NH*HO>0.5 && CO*HO<-0.5){
              isHB=1;
            }
	  }
	}
	if(found) n1++;
	if(isHB) n2++;
	if(found && !isHB)n1r++;
	if(isHB && !found)n2r++;
	*/
	cm_mm[i][j]+=(double)found;
      }
    }
  }
}

int main(int argc, char* argv[]){
  if(argc < 7){
    cout << "usage: movie2freq.linux movie(bin) inputMethod inputFile startFrame nFrame output" << endl;
    cout << " inputMethod : 0/1" << endl;
    cout << "    0 --- input seq file only" << endl;
    cout << "    1 --- input PDB file but using SEQ only" << endl;
    exit(1);
  }
  pseudoPep p(argv[3], (input_t)atoi(argv[2]));
  int natoms = p.makeIndex();
  //set the atom type for each atom in the peptide
  pepSetAtomType(p);
  //make a table of interaction radius
  double* IR = createIR();
  
  int startFrame=atoi(argv[4]);
  int nFrame = atoi(argv[5]);

  movie m(argv[1]);
  if(m.getErrorStatus()){
    cout << m.getMsg() << endl;
  }
  double** r = m.getCoords();
  double boundary = m.getDimensions()[0];
  short** bonds = m.getBonds();
  
  
  int natoms_m = m.getNAtoms();
  if((natoms_m%natoms)){
    cout << "mismach of input pdb and movie" << endl;
  }
  int nProts = natoms_m / natoms;

  double*** freq_ss = new double**[nProts];
  for(int i=0; i<nProts; i++){
    freq_ss[i] = new double*[natoms];
    freq_ss[i][0] = new double[(natoms*(natoms-1))/2](0);
    for(int j=1; j<natoms; j++) freq_ss[i][j] = freq_ss[i][j-1]+j-1;
  }

  double*** freq_mm = new double**[nProts];
  for(int i=0; i<nProts; i++){
    freq_mm[i] = new double*[natoms];
    freq_mm[i][0] = new double[(natoms*(natoms-1))/2](0);
    for(int j=1; j<natoms; j++) freq_mm[i][j] = freq_mm[i][j-1]+j-1;
  }

  m.jumpto(startFrame);
  for(int iFrame=0; iFrame<nFrame; iFrame++){
    {
      int& np = nProts;
      int& nres = natoms;
      int& npart = natoms_m;

      double cm[np][3];
      for(int ip=0; ip<np; ip++){
	cm[ip][0]=r[ip*nres][0];
	cm[ip][1]=r[ip*nres][1];
	cm[ip][2]=r[ip*nres][2];
	for(int i=1+ip*nres; i<nres*(ip+1); i++){
	  if(r[i][0]-r[i-1][0] > boundary/2.0f ) r[i][0]-=boundary;
	  if(r[i][0]-r[i-1][0] < -boundary/2.0f) r[i][0]+=boundary;
	  if(r[i][1]-r[i-1][1] > boundary/2.0f ) r[i][1]-=boundary;
	  if(r[i][1]-r[i-1][1] < -boundary/2.0f) r[i][1]+=boundary;
	  if(r[i][2]-r[i-1][2] > boundary/2.0f ) r[i][2]-=boundary;
	  if(r[i][2]-r[i-1][2] < -boundary/2.0f) r[i][2]+=boundary;
	  cm[ip][0]+=r[i][0];
	  cm[ip][1]+=r[i][1];
	  cm[ip][2]+=r[i][2];
	} 
	for(int i=0; i<3; i++)cm[ip][i]/=(double)nres;
      }
      
      int alloc[np];
      for(int ip=0; ip<np; ip++)alloc[ip]=1;
      int next=0;
      alloc[next]=0;
      double cmCurr[3];
      
      for(int ncount=1; ncount<np; ncount++){
	for(int j=0; j<3; j++) cmCurr[j]=0;
	int ialloc=0;
	for(int i=0; i<np; i++)if(!alloc[i]){
	  ialloc++;
	  for(int j=0;j<3;j++)cmCurr[j]+=cm[i][j];
	}
	for(int j=0; j<3; j++) cmCurr[j]/=(double)ialloc;
	
	double dist=0, min=3000*pow(boundary,2);
	int imin;
	for(int i=0; i<np; i++)if(alloc[i]){
	  for(int j=0; j<3;j++){
	    double temp = fabs(cm[i][j]-cmCurr[j]);
	    if(temp>boundary/2.0f)temp-=boundary;
	    dist+=temp*temp;
	  }
	  if(min>dist){
	    min = dist;
	    imin = i;
	  }
	}
	
	if(cm[imin][0]-cmCurr[0]>boundary/2.0f){
	  cm[imin][0]-=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][0]-=boundary;
	}
	if(cm[imin][0]-cmCurr[0]<-boundary/2.0f){
	  cm[imin][0]+=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][0]+=boundary;
	}
	
	if(cm[imin][1]-cmCurr[1]>boundary/2.0f){
	  cm[imin][1]-=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][1]-=boundary;
	}
	if(cm[imin][1]-cmCurr[1]<-boundary/2.0f){
	  cm[imin][1]+=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][1]+=boundary;
	}
	
	if(cm[imin][2]-cmCurr[2]>boundary/2.0f){
	  cm[imin][2]-=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][2]-=boundary;
	}
	if(cm[imin][2]-cmCurr[2]<-boundary/2.0f){
	  cm[imin][2]+=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][2]+=boundary;
	}
	alloc[imin]=0;
      }
      
      for(int i=0; i<3; i++)cmCurr[i]=cm[0][i];
      for(int ip=1; ip<np;ip++)
	for(int i=0;i<3;i++)cmCurr[i]+=cm[ip][i];
      for(int i=0; i<3; i++)cmCurr[i]/=(double)np;
      for(int i=0; i<npart; i++){
	r[i][0]-=(cmCurr[0]-boundary/2.0);
	r[i][1]-=(cmCurr[1]-boundary/2.0);
	r[i][2]-=(cmCurr[2]-boundary/2.0);
      }
    }

    int index=0;
    for(int iProts=0; iProts<nProts; iProts++){
      for(int iAA=0; iAA<p.getLength(); iAA++){
	pseudoAA* theAA = p.getResidue(iAA);
	*theAA->getN()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	index++;
	*theAA->getC()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	index++;
	*theAA->getO()->getR() =  vec(r[index][0], r[index][1], r[index][2]);
	index++;
	*theAA->getCA()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	index++;
	if(theAA->getCB()){
	  *theAA->getCB()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	  index++;
	}
	if(theAA->getG()){
	  *theAA->getG()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	  index++;
	}
	if(theAA->getG2()){
	  *theAA->getG2()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	  index++;
	}
      }
      p.constructHydrogen();
      //p.printPDB(out);
      //if(iProts<nProts-1) out<< "TER" << endl;
      pep2CM(p, IR, bonds, freq_ss[iProts], freq_mm[iProts]);
    }
    //out << "END" << endl;
    m.nextFrame();
  }
  ofstream out(argv[6]);
  for(int i=0; i<p.getLength();i++){
    for(int j=0; j<i-1; j++){
      out << i+1 << " " << j+1 << " " << freq_ss[0][i][j]/(double)nFrame << " "
	   << freq_mm[0][i][j]/(double)nFrame << endl;
    }
  }
  out.close();

  cout << n1 << " " << n1r << endl;
  cout << n2 << " " << n2r << endl;
}
