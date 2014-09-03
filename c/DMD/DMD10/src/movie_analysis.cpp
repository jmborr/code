#include "pseudoPep.h"
#include "pseudoAA.h"
#include "movie.h"
#include "krms.h"
#include "pseudoTXT.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

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

int compare(char* str1, char* str2, int len){
  if(strlen(str1)<len && strlen(str2)<len)return 0;
  for(int i=0;i<len; i++){
    if(str1[i]!=str2[i]) return 0;
  }
  return 1;
}

int readCA(char* file, double**& tmp){
  ifstream in(file);
  char line[100];
  vector<double> tmp_array;
  while(in.getline(line, 100)){
    if(compare(line, "ATOM  ", 6) && compare(&line[12], " CA ", 4)){//CA line
      double x, y, z;
      sscanf(&line[30], "%lf%lf%lf", &x, &y, &z);
      tmp_array.push_back(x);
      tmp_array.push_back(y);
      tmp_array.push_back(z);
    }
  }
  in.close();
  int length = tmp_array.size()/3;
  tmp = new double*[length];
  tmp[0] = new double[length*3];
  for(int i=1; i<length; i++) tmp[i] = tmp[i-1]+3;
  for(int i=0; i<length*3; i++) tmp[0][i] = tmp_array[i];
  return length;
}

void shiftR(double** tmp, int len){
  double* r=new double[3](0);
  for(int i=0;i<3; i++)r[i]=0;
  for(int i=0; i<len; i++)
    for(int j=0; j<3; j++){
      r[j] += tmp[i][j];
    }
  
  for(int i=0; i<3; i++) r[i]/=(double)len;
  
  for(int i=0; i<len; i++)
    for(int j=0; j<3; j++)
      tmp[i][j]-=r[j];
  delete []r;
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
    if(ia->getID()==TRP) ia->getD()->setType(_TRP_CD_);
    if(ia->getID()==ARG) ia->getD()->setType(_ARG_CD_);
    if(ia->getID()==LYS) ia->getD()->setType(_LYS_CD_);
  }
}


double* createIR(){
  double* table = new double[_TRP_CD_ - _CYS_CB_ +1];
  for(int i=_CYS_CB_; i<=_PRO_CB_; i++){
    if(dmd_atom_hp[i-_CYS_CB_]==APOLOR ||
       dmd_atom_hp[i-_CYS_CB_]==AMPHI){
      table[i-_CYS_CB_]=CB_IR;
    }
    else{
      table[i-_CYS_CB_]=-10000;
    }
  }
  
  for(int i=_CYS_CG_; i<=_PRO_CG_; i++){
    table[i-_CYS_CB_] = G1_CONST[i-_CYS_CG_][G1_IR];
  }
  
  table[_ILE_CG2_-_CYS_CB_]=ILE_G2[G2_IR];
  table[_THR_CG2_-_CYS_CB_]=THR_G2[G2_IR];
  table[_VAL_CG2_-_CYS_CB_]=VAL_G2[G2_IR];
  
  table[_ARG_CD_-_CYS_CB_]=ARG_D[D_IR];
  table[_LYS_CD_-_CYS_CB_]=LYS_D[D_IR];
  table[_TRP_CD_-_CYS_CB_]=TRP_D[D_IR];
  return table;
}


double** createEMap(){
  int ndim = _TRP_CD_ - _CYS_CB_ + 1;
  double** table = new double*[ndim];
  table[0] = new double[ndim*ndim](0);
  for(int i=1; i<ndim; i++) table[i] = table[i-1]+ndim;
  //h-h interaction
  for(int i=0; i<ndim; i++){
    for(int j=0; j<ndim; j++){
      if(dmd_atom_hp[i]==APOLOR && dmd_atom_hp[j]==APOLOR){
	table[i][j] = HP_E;
      }
      else if((dmd_atom_hp[i]==APOLOR && dmd_atom_hp[j]==AMPHI)||
	      (dmd_atom_hp[i]==AMPHI && dmd_atom_hp[j]==APOLOR)){
	table[i][j] = HP_E*3.0/4.0;
      }
    }
  }
  //aromatic interaction
  table[_PHE_CG_-_CYS_CB_][_PHE_CG_-_CYS_CB_] = AROMATIC_E;
  table[_PHE_CG_-_CYS_CB_][_TRP_CD_-_CYS_CB_] = AROMATIC_E;
  table[_PHE_CG_-_CYS_CB_][_TYR_CG_-_CYS_CB_] = AROMATIC_E;

  table[_TRP_CD_-_CYS_CB_][_PHE_CG_-_CYS_CB_] = AROMATIC_E;
  table[_TRP_CD_-_CYS_CB_][_TRP_CD_-_CYS_CB_] = AROMATIC_E;
  table[_TRP_CD_-_CYS_CB_][_TYR_CG_-_CYS_CB_] = AROMATIC_E;

  table[_TYR_CG_-_CYS_CB_][_PHE_CG_-_CYS_CB_] = AROMATIC_E;
  table[_TYR_CG_-_CYS_CB_][_TRP_CD_-_CYS_CB_] = AROMATIC_E;
  table[_TYR_CG_-_CYS_CB_][_TYR_CG_-_CYS_CB_] = AROMATIC_E;

  table[_TRP_CD_-_CYS_CB_][_PRO_CG_-_CYS_CB_] = AROMATIC_E;
  table[_PHE_CG_-_CYS_CB_][_PRO_CG_-_CYS_CB_] = AROMATIC_E;
  table[_TYR_CG_-_CYS_CB_][_PRO_CG_-_CYS_CB_] = AROMATIC_E;
  table[_PRO_CG_-_CYS_CB_][_TRP_CD_-_CYS_CB_] = AROMATIC_E;
  table[_PRO_CG_-_CYS_CB_][_PHE_CG_-_CYS_CB_] = AROMATIC_E;
  table[_PRO_CG_-_CYS_CB_][_TYR_CG_-_CYS_CB_] = AROMATIC_E;
  //salt bridge interaction
  table[_ARG_CD_-_CYS_CB_][_ARG_CD_-_CYS_CB_] = -SB_E/4.0;
  table[_ARG_CD_-_CYS_CB_][_LYS_CD_-_CYS_CB_] = -SB_E/4.0;
  table[_ARG_CD_-_CYS_CB_][_ASP_CG_-_CYS_CB_] = SB_E;
  table[_ARG_CD_-_CYS_CB_][_GLU_CG_-_CYS_CB_] = SB_E;

  table[_LYS_CD_-_CYS_CB_][_ARG_CD_-_CYS_CB_] = -SB_E/4.0;
  table[_LYS_CD_-_CYS_CB_][_LYS_CD_-_CYS_CB_] = -SB_E/4.0;
  table[_LYS_CD_-_CYS_CB_][_ASP_CG_-_CYS_CB_] = SB_E;
  table[_LYS_CD_-_CYS_CB_][_GLU_CG_-_CYS_CB_] = SB_E;


  table[_ASP_CG_-_CYS_CB_][_ARG_CD_-_CYS_CB_] = SB_E;
  table[_ASP_CG_-_CYS_CB_][_LYS_CD_-_CYS_CB_] = SB_E;
  table[_ASP_CG_-_CYS_CB_][_ASP_CG_-_CYS_CB_] = -SB_E/4.0;
  table[_ASP_CG_-_CYS_CB_][_GLU_CG_-_CYS_CB_] = -SB_E/4.0;
  
  table[_GLU_CG_-_CYS_CB_][_ARG_CD_-_CYS_CB_] = SB_E;
  table[_GLU_CG_-_CYS_CB_][_LYS_CD_-_CYS_CB_] = SB_E;
  table[_GLU_CG_-_CYS_CB_][_ASP_CG_-_CYS_CB_] = -SB_E/4.0;
  table[_GLU_CG_-_CYS_CB_][_GLU_CG_-_CYS_CB_] = -SB_E/4.0;
  return table;
}


double getE(pseudoPep& p, const double* IR_table, short** bonds, double** E, 
	    int& n_hp, int& n_ha, int& n_ar, int& n_pro, int& n_sb, int& n_hb ){
  int natom = p.getLength();
  int nHB = 0;
  int nContact=0;
  n_hp=0; n_ha=0; n_ar=0; n_pro=0; n_sb=0; n_hb=0;
  double totalE=0;
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
	    nContact++;
	    totalE += E[ia_side[is].getType()-_CYS_CB_][ja_side[js].getType()-_CYS_CB_];
	    /*determine the categorized energies*/
	    int ta = ia_side[is].getType();
	    int tb = ja_side[js].getType();
	    if((dmd_atom_hp[ta - _CYS_CB_]==APOLOR && dmd_atom_hp[tb - _CYS_CB_]==APOLOR)||
	       (dmd_atom_hp[ta - _CYS_CB_]==APOLOR && dmd_atom_hp[tb - _CYS_CB_]==AMPHI) ||
	       (dmd_atom_hp[ta - _CYS_CB_]==AMPHI && dmd_atom_hp[tb - _CYS_CB_]==APOLOR) ){
	      /*aromatic*/
	      if((ta==_PHE_CG_ || ta==_TRP_CD_ || ta==_TYR_CG_) && (tb==_PHE_CG_ || tb==_TRP_CD_ || tb==_TYR_CG_ || tb==_PRO_CG_) ||
		 (tb==_PHE_CG_ || tb==_TRP_CD_ || tb==_TYR_CG_) && (ta==_PHE_CG_ || ta==_TRP_CD_ || ta==_TYR_CG_ || ta==_PRO_CG_)){
		if(ta!=_PRO_CG_ && tb!=_PRO_CG_) 
		  n_ar++;
		else
		  n_pro++;
	      }
	      /*normal HP*/
	      else{
		if((dmd_atom_hp[ta - _CYS_CB_]==APOLOR && dmd_atom_hp[tb - _CYS_CB_]==APOLOR))
		  n_hp++;
		else
		  n_ha++;
	      }
	    }
	    /*salt bridge*/
	    if((ta==_ASP_CG_||ta==_GLU_CG_) && (tb==_ARG_CD_||tb==_LYS_CD_)||
	       (tb==_ASP_CG_||tb==_GLU_CG_) && (ta==_ARG_CD_||ta==_LYS_CD_)){
	      n_sb++;
	    }
	  }
	}
      }
      
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
	nHB+=found;
	if(found){
	  totalE += EHB_M_M;
	}
	int iO = ia->getO()->getIndex();
	int jN = ja->getN()->getIndex();
	ii=0;
	found = 0;
	while(bonds[jN-1][ii]>=0){
	  if(bonds[jN-1][ii]==iO-1){
	    found=1;
	    break;
	  }
	  ii++;
	}
	nHB+=found;
	if(found){
	  totalE += EHB_M_M;
	}
      }
    }
  }
  n_hb=nHB;
  return totalE;
}

int main(int argc, char* argv[]){
  if(argc < 9){
    cout << "usage: movie_analysis.linux movie(bin) E_para_f inputMethod inputFile refPDB startFrame nFrame output" << endl;
    cout << " inputMethod : 0/1" << endl;
    cout << "    0 --- input seq file only" << endl;
    cout << "    1 --- input PDB file but using SEQ only" << endl;
    exit(1);
  }

  process_para(argv[2]);

  pseudoPep p(argv[4], (input_t)atoi(argv[3]));
  int natoms = p.makeIndex();
  //set the atom type for each atom in the peptide
  pepSetAtomType(p);
  //make a table of interaction radius
  double* IR = createIR();
  //make a table of interaction energy
  double** map = createEMap();
  
  //read the CA of the reference PDB
  double ** ref_r = NULL;
  int nca_ref = readCA(argv[5], ref_r);

  if(nca_ref!=p.getLength()){
    cout << "The refPDB and Movie are matching" << endl;
    exit(1);
  }
  shiftR(ref_r, nca_ref);
  
  double** the_r = new double*[p.getLength()];
  the_r[0] = new double[p.getLength()*3](0);
  for(int i=1; i<p.getLength(); i++){
    the_r[i] = the_r[i-1]+3;
  }
  
  int startFrame=atoi(argv[6]);
  int nFrame = atoi(argv[7]);

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
  
  m.jumpto(startFrame);
  ofstream out(argv[8]);
  char ch_dummy[100];
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
	if(theAA->getD()){
	  *theAA->getD()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	  index++;
	}
      }
      p.constructHydrogen();
      for(int i=0; i<p.getLength(); i++){
	the_r[i][0] = p.getResidue(i)->getCA()->getR()->getX();
	the_r[i][1] = p.getResidue(i)->getCA()->getR()->getY();
	the_r[i][2] = p.getResidue(i)->getCA()->getR()->getZ();
      }
      shiftR(the_r, p.getLength());
      double rmsd = get_rms(the_r, ref_r, p.getLength());
      int n_hp, n_ha, n_ar, n_pro, n_sb, n_hb;
      double e_total = getE(p, IR, bonds, map, n_hp, n_ha, n_ar, n_pro, n_sb, n_hb);
      sprintf(ch_dummy, "%8.5lf %8.5lf %4ld %4ld %4ld %4ld %4ld %4ld", sqrt(rmsd/p.getLength()), e_total, 
	      n_hp, n_ha, n_ar, n_pro, n_sb, n_hb);
      out << ch_dummy << endl;
    }
    m.nextFrame();
  }
  out.close();
  
}
