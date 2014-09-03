/********************************************************
The classes used to manipulate PDB files.
-----------Feng Ding-------------------------------------
-----------fding@buphy.bu.edu----------------------------
*********************************************************/
#ifndef _PDBLIB_H_
#define _PDBLIB_H_
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;
const double PI = 4.0*atan(1.0);
const double rpi = PI/180;
const static double INF = 1.0e+38;

void newLine(ostream& out){
  out << endl;
}

string dtostr(int d){
  char s[21];
  sprintf(s,"%ld",d);
  return string(s);
}

class histogram{
  int* his;
  double min;
  double max;
  int nbin;
  double step;
  double sum;
  double sqrSum;
  int counter;
  double theMin;
  double theMax;
 public:
  histogram(double imin, double imax, int inbin){
    min = imin;
    max = imax;
    nbin = inbin;
    step = (max-min)/(double)nbin;
    his = new int[nbin];
    for(int i=0; i<nbin; i++)
      his[i]=0;
    counter=0;
    sqrSum=0;
    sum=0;
    theMax = min;
    theMin = max;
  }
  
  void addToBin(double a){
    if(a > min && a <max){
      int i = (int)((a - min)/step);
      his[i]++;
      sum+=a;
      sqrSum+=a*a;
      counter++;
      if(a > theMax) theMax = a;
      if(a < theMin) theMin = a;
    }
  }
  
  int* getHis(){
    return his;
  }
  
  int getNBin(){
    return nbin;
  }
  
  double getMin(){
    return min;
  }
  
  double getMax(){
    return max;
  }
  
  double getStep(){
    return step;
  }
  
  double getTheMax(){
    return theMax;
  }
  
  double getTheMin(){
    return theMin;
  }
  
  void print(ostream& out){
    for(int i=0; i<nbin; i++){
      out << min+i*step << "\t" << his[i] << endl;
    }
  }
  
  void normalPrint(ostream& out){
    for(int i=0; i<nbin; i++){
      if(counter>0)
	out << min+i*step << "\t" << his[i]/(double)counter << endl;
      else
	out << min+i*step << "\t0" << endl;
    }
  }
  
  double getAve(){
    return sum/(double)counter;
  }
  
  double getStd(){
    double tmp1 = sum/(double)counter;
    double tmp2 = sqrSum/(double)counter;
    return sqrt(tmp2 - tmp1*tmp1);
  }
  
  int getNCount(){
    return counter;
  }
};

class cstr{
  char* str;
 public:
  cstr(const string& a){
    int len = a.length();
    str = (char*)malloc((len+1)*sizeof(char));
    a.copy(str,len);
    str[len]='\0';
  }
  ~cstr(){
    delete []str;
  }
  char* getStr(){
    return str;
  }
};


class vec{
  double x;
  double y;
  double z;
 public:
  vec(double ix=0, double iy=0, double iz=0){
    x = ix;
    y = iy;
    z = iz;
  }
  
  vec(const vec& old){ //clone
    x = old.x;
    y = old.y;
    z = old.z;
  }
  
  ~vec(){
    //cout << "delete vector" << endl;
  }
  
  // operator overlading
  friend const vec 
    operator+(const vec& left,
	      const vec& right);
  
  friend const vec 
    operator-(const vec& left,
	      const vec& right);
  
  friend const vec 
    operator*(const vec& left,
	      const double& right);
  
  friend const vec 
    operator*(const double& left,
	      const vec& right);
  
  friend const double 
    operator*(const vec& left,
	      const vec& right); //dot product
  
  friend const vec 
    operator^(const vec& left,
	      const vec& right); //cross product
  
  friend const vec 
    operator/(const vec& left,
	      const double& right);
  
  vec& 
    operator=(const vec& right){
    if(this == &right) return *this;
    x = right.x;
    y = right.y;
    z = right.z;
    return *this;
  }
  
  vec&
    operator+=(const vec& right){
    x += right.x;
    y += right.y;
    z += right.z;
    return *this;
  }
  
  vec&
    operator-=(const vec& right){
    x -= right.x;
    y -= right.y;
    z -= right.z;
    return *this;
  }
  
  vec&
    operator*=(double right){
    x *= right;
    y *= right;
    z *= right;
    return *this;
  }
  
  vec&
    operator/=(double right){
    x /= right;
    y /= right;
    z /= right;
    return *this;
  }
  
  void print(ostream& out)const{
    out << "(" << x << ", " << y << ", " << z << ")";
  }
  
  double getX(){
    return x;
  }
  
  double getY(){
    return y;
  }
  
  double getZ() const {
    return z;
  }
  
  void setX(double ix){
    x = ix;
  }
  
  void setY(double iy){
    y = iy;
  }
  
  void setZ(double iz){
    z = iz;
  }
  
  void addX(double ix){
    x += ix;
  }
  
  void minusX(double ix){
    x -= ix;
  }
  
  void addY(double iy){
    y += iy;
  }
  
  void minusY(double iy){
    y -= iy;
  }
  
  void addZ(double iz){
    z += iz;
  }
  
  void minusZ(double iz){
    z -= iz;
  }
  
  double getDist(const vec& a){
    double tmp=0;
    tmp+=(x-a.x)*(x-a.x);
    tmp+=(y-a.y)*(y-a.y);
    tmp+=(z-a.z)*(z-a.z);
    return sqrt(tmp);
  }
  
  double getMag(){
    return sqrt(x*x+y*y+z*z);
  }
  
  void Normalize(){
    double tmp = sqrt(x*x + y*y + z*z);
    if( tmp > 0){
      x /= tmp;
      y /= tmp;
      z /= tmp;
    }
  }
  
  vec norm() const{
    double tmp=sqrt(x*x + y*y + z*z);
    if(tmp>0){
      return vec(x/tmp, y/tmp, z/tmp);
    }
    else{
      return vec(INF, INF, INF);
    }
  }
  
  void Rotate(const vec& p1,
	      const vec& p2,
	      double theta){
    /*-------------------------
      this rotation is make by unchange the position of P1 and P2, and 
      use the axis defined by the so defined vector from p1 to p2 TO rotate
      theta angle anti-clockwise.
      -------------------------*/
    vec axis = p2 - p1;
    axis.Normalize();
    
    vec tmp = vec(axis);
    vec tmp2 = (*this) - p2;
    double product = tmp2 * tmp;
    tmp *= product;
    
    vec pivot = tmp + p2;
    tmp = (*this) - pivot;
    tmp2 = axis^tmp;
    
    tmp *= cos(theta);
    tmp2 *= sin(theta);
    
    *this = tmp + tmp2 + pivot;
  }

  void EularRotate(double** m){
    /*---------
      rotate the vector according to the input matrix defined
      --------*/
    double r[3];
    r[0] = m[0][0]*x+m[0][1]*y+m[0][2]*z;
    r[1] = m[1][0]*x+m[1][1]*y+m[1][2]*z;
    r[2] = m[2][0]*x+m[2][1]*y+m[2][2]*z;
    *this = vec(r[0], r[1], r[2]);
  }
};

// operator overlading
inline const vec 
operator+(const vec& left,
	  const vec& right){
  return vec(left.x + right.x,
	     left.y + right.y,
	     left.z + right.z);
}

inline const vec 
operator-(const vec& left,
	  const vec& right){
  return vec(left.x - right.x,
	     left.y - right.y,
	     left.z - right.z);
}

inline const vec 
operator*(const vec& left,
	  const double& right){
  return vec(left.x*right,
	     left.y*right,
	     left.z*right);
}

inline const vec 
operator*(const double& left,
	  const vec& right){
  return vec(left*right.x,
	     left*right.y,
	     left*right.z);
}

inline const double 
operator*(const vec& left,
	  const vec& right){                  //dot product
  return (left.x*right.x + left.y*right.y + left.z*right.z);
}

inline const vec 
operator^(const vec& left,
	  const vec& right){                  //cross product
  return vec(left.y*right.z - left.z*right.y,
	     left.z*right.x - left.x*right.z,
	     left.x*right.y - left.y*right.x);
}

inline const vec 
operator/(const vec& left,
	  const double& right){
  return vec(left.x/right,
	     left.y/right,
	     left.z/right);
}

class atom{
  char name;
  vec* r;
  int index;
 public:
  atom(char iname, double x, double y, double z){
    name = iname;
    r = new vec(x,y,z);
    index=0;
  }
  
  atom(char iname, const vec& rvec){
    name = iname;
    r = new vec(rvec);
    index = 0;
  }

  atom(const atom& old){
    name=old.name;
    r=new vec(*old.r);
    index = old.index;
  }
  
  atom(){
    r = new vec(0,0,0);
    index = 0;
  }

  ~atom(){
    delete r;
    //cout << "end deleteing atom r" << endl;
  }
  
  const char getName(){
    return name;
  }
  
  void setName(char ichar){
    name = ichar;
  }

  vec* getR(){
    return r;
  }
  
  void setIndex(int i){
    index =i;
  }
  
  int getIndex(void){
    return index;
  }
  
  void print(ostream& out){
    out << "ATOM: " <<  name;
    r->print(out);
  }
  
  double getDist(atom& a){
    return r->getDist(*(a.r));
  }

  atom& 
    operator=(const atom& right){
    if(this == &right) return *this;
    name = right.name;
    *r = *(right.r);
    return *this;
  }
};

const static string aa_3c[21]={
  "ALA","CYS","ASP","GLU","PHE",
  "GLY","HIS","ILE","LYS","LEU",
  "MET","ASN","PRO","GLN","ARG",
  "SER","THR","VAL","TRP","TYR",
  "MO2"
};

const static string aa_1c[21]={
  "A"  ,  "C",  "D",  "E",  "F",
  "G"  ,  "H",  "I",  "K",  "L",
  "M"  ,  "N",  "P",  "Q",  "R",
  "S"  ,  "T",  "V",  "W",  "Y",
  "O"
};

const double Roseman[21] = {  
  //Ala:
  0.390,  \
  //Cys:
  0.250, \
  //Asp: -
  -3.810, \
  //Glu: -
  -2.910, \
  //Phe:
  2.270, \
  //Gly:
  0.000, \
  //His: -
  -0.640, \
  //Ile:
  1.820, \
  //Lys: +
  -2.770, \
  //Leu:
  1.820, \
  //Met:
  0.960, \
  //Asn:
  -1.910, \
  //Pro:
  0.990, \
  //Gln:
  -1.300, \
  //Arg: +
  -3.950, \
  //Ser:
  -1.240, \
  //Thr:
  -1.000, \
  //Val:
  1.300,  \
  //Trp:
  2.130, \
  //Tyr:
  1.470, \
  //MO2:
  -1.910/*this is a make-up value, not the real one*/};

const double CG_Charge[21] = {  
  //Ala:
  0.0,  \
  //Cys:
  0.0, \
  //Asp: -
  -1.0, \
  //Glu: -
  -1.0, \
  //Phe:
  0.0, \
  //Gly:
  0.0, \
  //His: +
  0.0, \
  //Ile:
  0.0, \
  //Lys: +
  1.0, \
  //Leu:
  0.0, \
  //Met:
  0.0, \
  //Asn:
  0.0, \
  //Pro:
  0.0, \
  //Gln:
  0.0, \
  //Arg: +
  1.0, \
  //Ser:
  0.0, \
  //Thr:
  0.0, \
  //Val:
  0.0,  \
  //Trp:
  0.0, \
  //Tyr:
  0.0, \
  //MO2: oxydized MET
  -1.0/*a make-up value, not the real one*/};

bool check(string name){
  int len = name.size();
  if(len == 3)
    for(int i=0; i<21; i++){
      if(name == aa_3c[i]) return true;
    }
  else if(len == 1)
    for(int i=0; i<21; i++){
      if(name == aa_1c[i]) return true;
    }
  return false;
}

int Name2Int(string name){
  int len = name.size();
  if(name.size()==3)
    for(int i=0; i<21; i++){
      if(name==aa_3c[i])return i;
    }
  else if(len==1)
    for(int i=0; i<21; i++){
      if(name==aa_1c[i])return i;
    }
  else return -1;
}

string Int2Name(int i){
  return aa_3c[i];
}

class aa{
  string name;
  atom* N;
  atom* C;
  atom* O;
  vector<atom*> heavyAtom;//CA, [CB, [C,O,N,S]]
 public:
  aa(string iname){
    name = iname;
    N = NULL;
    C = NULL;
    O = NULL;
  }
  
  aa(const aa& old){
    name = old.name;
    if(old.N)N=new atom(*old.N);
    if(old.C)C=new atom(*old.C);
    if(old.O)O=new atom(*old.O);
    for(int i=0; i<(old.heavyAtom).size(); i++){
      heavyAtom.push_back(new atom(*old.heavyAtom[i]));
    }
  }
  
  ~aa(){
    if(N!=NULL) delete N;
    if(C!=NULL) delete C;
    if(O!=NULL) delete O;
    for(int i=0; i<heavyAtom.size(); i++) delete heavyAtom[i];
    /*
      vector<atom*>::iterator w;
      for(w=heavyAtom.begin(); w!=heavyAtom.end(); w++){
      delete *w;
      }
    */
    heavyAtom.clear();
    //cout << "delete aa" << endl;
  }
  
  void addAtom(char ch, double x, double y, double z){
    heavyAtom.push_back(new atom(ch, x, y, z));
  }
  
  string getName()const{
    return name;
  }
  
  atom* getN()const{
    return N;
  }
  
  atom* getC()const{
    return C;
  }
  
  atom* getO()const{
    return O;
  }
  
  void setN(double x, double y, double z){
    N = new atom('N', x, y ,z);
  }
  
  void setC(double x, double y, double z){
    C = new atom('C', x, y ,z);
  }
  
  void setO(double x, double y, double z){
    O = new atom('O', x, y ,z);
  }
  
  atom* getAtomAt(int i)const{
    if(i >= 0 && i < heavyAtom.size())return heavyAtom[i];
  }
  
  int getSize()const{
    return heavyAtom.size();
  }
  
  void print(ostream& out){
    out << name << endl;
    out << "Peptide:" << endl;
    if(N!=NULL){N->print(out); newLine(out);}
    if(C!=NULL){C->print(out); newLine(out);}
    if(O!=NULL){O->print(out); newLine(out);}
    out << "CA + side chaine:" << endl;
    vector<atom*>::iterator w;
    for(w=heavyAtom.begin(); w!=heavyAtom.end(); w++){
      (*w)->print(out) ;
      newLine(out);
    }
  }
  bool SizeCheck()const;

  aa& 
    operator=(const aa& old){
    if(this == &old) return *this;
    name = old.name;
    if(old.N)N=new atom(*old.N);
    if(old.C)C=new atom(*old.C);
    if(old.O)O=new atom(*old.O);
    for(int i=0; i<(old.heavyAtom).size(); i++){
      heavyAtom.push_back(new atom(*old.heavyAtom[i]));
    }
    return *this;
  }
};

bool aa::SizeCheck()const{
  int natom = heavyAtom.size();
  switch(Name2Int(name)){
  case 0:
    return (natom == 2);
    break;
  case 1:
    return (natom == 3);
    break;
  case 2:
    return (natom == 5);
    break;
  case 3:
    return (natom == 6);
    break;
  case 4:
    return (natom == 8);
    break;
  case 5:
    return (natom == 1);
    break;
  case 6:
    return (natom == 7);
    break;
  case 7:
    return (natom == 5);
    break;
  case 8:
    return (natom == 6);
    break;
  case 9:
    return (natom == 5);
    break;
  case 10:
    return (natom == 5);
    break;
  case 11:
    return (natom == 5);
    break;
  case 12:
    return (natom == 4);
    break;
  case 13:
    return (natom == 6);
    break;
  case 14:
    return (natom == 8);
    break;
  case 15:
    return (natom == 3);
    break;
  case 16:
    return (natom == 4);
    break;
  case 17:
    return (natom == 4);
    break;
  case 18:
    return (natom == 11);
    break;
  case 19:
    return (natom == 9);
    break;
  case 20:/*MO2, oxydized MET*/
    return (natom == 7);
    break;
  default:
    return false;
  }
}

class pseudo_aa{
  string name;
  bool valid;
  vec* N;
  vec* ca;
  vec* C;
  vec* cb;
  vec* side;
 public:
  pseudo_aa(string iname){
    check(iname);
    name = iname;
    N=NULL;
    C=NULL;
    ca = NULL;
    cb = NULL;
    side = NULL;
  }
  
  pseudo_aa(const aa& iaa){
    name = iaa.getName();
    int len = iaa.getSize();
    valid = iaa.SizeCheck();
    if(valid){
      if(len==1){
	ca = new vec( *(iaa.getAtomAt(0)->getR()) );
	cb = NULL;
	side = NULL;
      }
      else if(len==2){
	ca = new vec( *(iaa.getAtomAt(0)->getR()) );
	cb = new vec( *(iaa.getAtomAt(1)->getR()) );
	side = NULL;
      }
      else{
	ca = new vec( *(iaa.getAtomAt(0)->getR()) );
	cb = new vec( *(iaa.getAtomAt(1)->getR()) ); 
	side = new vec(0, 0, 0);
	for(int i=2; i<len; i++){
	  (*side) += *(iaa.getAtomAt(i)->getR());
	}
	(*side) /= (double)(len-2);
      }
    }
    else{
      ca = NULL;
      cb = NULL;
      side = NULL;
    }
    if(!iaa.getN()){
      N=NULL;
      valid=false;
    }
    else N = new vec(*iaa.getN()->getR());
    if(!iaa.getC()){
      C=NULL;
      valid=false;
    }
    else C = new vec(*iaa.getC()->getR());
  }
  
  ~pseudo_aa(){
    if(N!=NULL)delete N;
    if(C!=NULL)delete C;
    if(ca!=NULL)delete  ca;
    if(cb!=NULL)delete  cb;
    if(side!=NULL)delete side;
  }
  
  void setCa(double x, double y, double z){
    ca = new vec(x, y, z);
  }
  
  void setCb(double x, double y, double z){
    cb = new vec(x, y, z);
  }
  
  void setSide(double x, double y, double z){
    side = new vec(x, y, z);
  }
  
  string getName(){
    return name;
  }
  
  vec* getN(){
    return N;
  }

  vec* getC(){
    return C;
  }

  vec* getCa(){
    return ca;
  }
  
  vec* getCb(){
    return cb;
  }
  
  vec* getSide(){
    return side;
  }
  
  bool getValid(){
    return valid;
  }
  
  void print(ostream& out){
    out << name << endl;
    if(N){
      out << "N\t:"; N->print(out); newLine(out);
    }
    out << "CA\t:"; ca->print(out); newLine(out);
    if(C){
      out << "C\t:"; C->print(out); newLine(out);
    }
    if(cb!=NULL){
      out << "CB\t:"; cb->print(out); newLine(out);
    }
    if(side!=NULL){
      out << "Side\t:"; side->print(out); newLine(out);
    }
  }
};

class protein{
  string PDBCode;
  int origin;           //1: original PDB, 0: reduced PDB
  vector<aa*> list;
  int errStatus;
  //error status of reading PDB:
  //0: everything is good!
  //1: chain broken
  //2: if it is an original PDB but not full presentation
  //     of atom in one of the Amino Acid
  //3: reading error: CA is not the first
  //4: reading error: CB is not the second
  //5: can not open the PDB file
  char** cm;
  string errMsg;
  int offset;
 public:
  protein(string fileName, int isOrigin=0, double cutoff=7.5);
  
  ~protein(){
    int size = list.size();
    for(int i=0; i<size; i++) delete list[i];
    list.clear();
    for(int i=0; i<size; i++) delete [] cm[i];
    delete [] cm;
      //cout << "delete protein" << endl;
  }
  
  int getLength(){
    return list.size();
  }
  
  string getName(){
    return PDBCode;
  }
  
  void addElement(aa* newAA){
    list.push_back(newAA);
  }
  
  int getErrStatus(){
    return errStatus;
  }
  
  aa* getElementAt(int i){
    if(i>=0 && i<list.size()) return list[i];
  }
  
  void print(ostream& out){
    out << "Protein: " << PDBCode << endl;
    out << "Length : " << list.size() << endl;
    vector<aa*>::iterator w;
    for(w=list.begin(); w!=list.end(); w++){
      (*w)->print(out);
    }      
  }

  double getRg(){
    int natom=0;
    vec cm(0,0,0);
    int len = getLength();
    for(int i=0; i<len; i++){
      aa* theAA=list[i];
      atom* theAtom;
      theAtom = theAA->getC();
      if(theAtom){
	natom++;
	cm += *(theAtom->getR());
      }
    }
    cm /= (double)natom;
    double d2=0,d;
    for(int i=0; i<len; i++){
      aa* theAA=list[i];
      atom* theAtom;
      theAtom = theAA->getC();
      d = theAtom->getR()->getDist(cm);
      d2+=d*d;
    }
    d2/=(double)natom;
    return sqrt(d2);
  }
  
  void printCM(ostream& out, char* type, double cutoff){
    int len = getLength();
    if(strcmp(type,"CA")==0){
      for(int i=0; i<len; i++){
	aa* iaa = list[i];
	atom* ia = iaa->getAtomAt(0);
	for(int j=i+2; j<len; j++){
	  aa* jaa = list[j];
	  atom* ja = jaa->getAtomAt(0);
	  if(ia->getDist(*ja) < cutoff )
	    out << i+1 << "\t" << j+1 << endl;
	}
      }
    }
    else if(strcmp(type, "CB")==0){
      for(int i=0; i<len; i++){
	aa* iaa = list[i];
	atom* ia;
	if(iaa->getName()=="GLY") ia = iaa->getAtomAt(0);
	else ia = iaa->getAtomAt(1);
	for(int j=i+2; j<len; j++){
	  aa* jaa = list[j];
	  atom* ja;
	  if(jaa->getName()=="GLY") ja = jaa->getAtomAt(0);
	  else ja = jaa->getAtomAt(1);
	  if(ia->getDist(*ja) < cutoff )
	    out << i+1 << "\t" << j+1 << endl;
	}
      }
    }
  }
  
  void errPrint(ostream& out){
    if(errStatus==1)
      out << PDBCode << ": Main chain is broken " << errMsg 
	  << "@" << offset+getLength() << endl;
    else if(errStatus==5)
      out << PDBCode << ": The file can not be open" << endl;
    else if(errStatus==2)
      out << PDBCode << ": No full presentation of " << errMsg 
	  << "@" << offset+getLength() << endl;
    else if(errStatus==3)
      out << PDBCode << ": CA is not the 1st ATOM of " 
	  << errMsg << "@" << offset+getLength() << endl;
    else if(errStatus==4)
      out << PDBCode << ": CB is not the 2nd ATOM of " 
	  << errMsg << "@" << offset+getLength() << endl;
  }
  
  int getNcAt(int iloc){
    int nc = 0;
    for(int i=0; i<list.size(); i++) nc += cm[iloc][i];
    return nc;
  }

  int getNcAll(){
    int nc = 0;
    for(int i=0; i<list.size(); i++)
      for(int j=i; j<list.size(); j++)
	nc += cm[i][j];
    return nc;
  }

  void makePDB(char* outfile, int shift=0){
    ofstream out(outfile);
    char field[81];
    int iatom=0;
    double x,y,z;
    aa *theAmino;
    for(int i=0; i<getLength(); i++){
      theAmino=getElementAt(i);
      cstr name(theAmino->getName());
      //N
      if(theAmino->getN()){
	iatom++;
	x=theAmino->getN()->getR()->getX();
	y=theAmino->getN()->getR()->getY();
	z=theAmino->getN()->getR()->getZ();
	sprintf(&field[0],"ATOM  %5d  N   %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00           C  ",
		iatom,name.getStr(),i+1+shift,x,y,z);
	out << field << endl;
      }
      //CA
      iatom++;
      x=theAmino->getAtomAt(0)->getR()->getX();
      y=theAmino->getAtomAt(0)->getR()->getY();
      z=theAmino->getAtomAt(0)->getR()->getZ();
      sprintf(&field[0],"ATOM  %5d  CA  %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00           C  ",
	      iatom,name.getStr(),i+1+shift,x,y,z);
      out << field << endl;
      //C
      if(theAmino->getC()){
	iatom++;
	x=theAmino->getC()->getR()->getX();
	y=theAmino->getC()->getR()->getY();
	z=theAmino->getC()->getR()->getZ();
	sprintf(&field[0],"ATOM  %5d  C   %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00           C  ",
		iatom,name.getStr(),i+1+shift,x,y,z);
	out << field << endl;
      }
      //O
      if(theAmino->getO()){
	iatom++;
	x=theAmino->getO()->getR()->getX();
	y=theAmino->getO()->getR()->getY();
	z=theAmino->getO()->getR()->getZ();
	sprintf(&field[0],"ATOM  %5d  O   %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00           C  ",
		iatom,name.getStr(),i+1+shift,x,y,z);
	out << field << endl;
      }
      if(theAmino->getName()!="GLY"){
	//CB
	iatom++;
	x=theAmino->getAtomAt(1)->getR()->getX();
	y=theAmino->getAtomAt(1)->getR()->getY();
	z=theAmino->getAtomAt(1)->getR()->getZ();
	sprintf(&field[0],"ATOM  %5d  CB  %3s  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00           C  ",
		iatom,name.getStr(),i+1+shift,x,y,z);
	out << field << endl;
      }
    }
    out.close();
  }
  
  int getPhiPsi(int i, double& phi, double& psi){
    if(i<0 || i>getLength()-1) return 0;
    aa* next;
    aa* prev;
    if(i<getLength()-1) 
      next = getElementAt(i+1);
    else 
      next=NULL;
    aa* curr = getElementAt(i);
    if(i>0) 
      prev = getElementAt(i-1);
    else 
      prev=NULL;

    vec N2_C1;
    vec C1_Ca1;
    vec Ca1_N1;
    
    C1_Ca1 = *(curr->getAtomAt(0)->getR()) - *(curr->getC()->getR());
    Ca1_N1 = *(curr->getN()->getR()) - *(curr->getAtomAt(0)->getR());
    //PHI
    if(prev){
      vec N1_C0  = *(prev->getC()->getR()) - *(curr->getN()->getR());
      double mod2 = Ca1_N1*Ca1_N1;
      vec r1 = C1_Ca1 - (C1_Ca1*Ca1_N1)*Ca1_N1/mod2;
      vec r2 = N1_C0  - ( N1_C0*Ca1_N1)*Ca1_N1/mod2;
      double cos = r1*r2;
      cos /= sqrt((r1*r1)*(r2*r2));
      double coef = r2*(Ca1_N1^r1);
      if(coef >0) phi = acos(cos);
      else phi = PI*2 - acos(cos);
    }
    else{
      phi=-PI*2;
    }
    //PSI
    if(next){
      N2_C1  = *(curr->getC()->getR()) - *(next->getN()->getR());
      double mod2 = C1_Ca1*C1_Ca1;
      vec r3 = Ca1_N1 - (Ca1_N1*C1_Ca1)*C1_Ca1/mod2;
      vec r4 = N2_C1  -  (N2_C1*C1_Ca1)*C1_Ca1/mod2;
      double cos = r3*r4;
      cos /= sqrt((r3*r3)*(r4*r4));
      double coef = r4*(C1_Ca1^r3);
      if(coef <0) psi = acos(cos);
      else psi = PI*2 - acos(cos);
    }
    else{
      psi=-PI*2;
    }
    return 1;
  }
};

protein::protein(string fileName, int isOrigin, double cutoff){
  char junk[21];
  int len = fileName.length();
  fileName.copy(junk,len); junk[len]='\0';
  ifstream in(junk);
  if(!in){
    errStatus = 5;
    return;
  }
    
  PDBCode = fileName;
  origin = isOrigin;
  errStatus = 0;
  
  string line;
  const string AtomField("ATOM  ");
  string tmp;
  while(getline(in,line)){
    tmp = line.substr(0, 6);
    if(tmp == AtomField) break;
  }                                 //read till ATOM fields
  
  char buff[81];
  int pnum = -999;
  int num;
  int side = 0;
  aa* pt;
  while(1){
    tmp=line.substr(22, 4);         //read the residue seq
    tmp.copy(buff, 4); buff[4]='\0';
    sscanf(buff, "%d", &num);
    if(pnum<0) {                    //first residue
      tmp = line.substr(17, 3);     //read residue NAME
      pt = new aa(tmp);             //create the new "aa" object
      pnum = num;
      offset = num;
    }
    else if(pnum==num);             //the same amino acid, do nothing
    else if(pnum==num-1){           //count a new amino acid
      if( origin==1 && !(pt->SizeCheck()) ){
	errMsg += " " + pt->getName() + "@" + dtostr(pnum); 
	errStatus = 2;
      }
      list.push_back(pt);           //push the old one
      tmp = line.substr(17, 3);     //read residue NAME
      pt = new aa(tmp);             //create the new "aa" object
      pnum = num;
      side = 0;                     //set side chain zero
    }
    else{                           //chain is somewhat broken
      errMsg = pt->getName();
      errStatus = 1;
      break;
    }
    
    double x, y , z;
    tmp=line.substr(30,24);
    tmp.copy(buff, 24); buff[24]='\0';
    sscanf(buff,"%lf%lf%lf", &x, &y, &z);//read x, y, z position
	
    char atomName[5];
    tmp=line.substr(12,4);
    tmp.copy(buff, 4); buff[4]='\0';
    sscanf(buff,"%s", atomName);         //read atom name
    
    string altLoc = line.substr(16,1);
    if( (altLoc == " " || altLoc == "A") && strcmp(atomName,"OXT")!=0){
      if(strcmp(atomName,"N")==0){         //Nitrogen 
	pt->setN(x, y, z);
      }
      else if(strcmp(atomName,"C")==0){    //Carbon
	pt->setC(x, y, z);
      }
      else if(strcmp(atomName,"O")==0){    //Oxygen
	pt->setO(x, y, z);
      }
      else if(strcmp(atomName,"CA")==0){   //C alpha
	if(pt->getSize()==0)               //check wheather it is the first
	  pt->addAtom('C',x,y,z);          //creat new atom object for CA
	else{
	  errMsg = pt->getName();
	  errStatus = 3;
	  break;
	}
      }
      else if(strcmp(atomName,"CB")==0){   //C beta
	if(pt->getSize()==1){              //check wheather it is the second
	  pt->addAtom('C',x,y,z);          //creat new atom object for CB
	  side = 1;                        //set side chain 1
	}
	else{
	  errMsg = pt->getName();
	  errStatus = 4;
	  break;
	}
      }
      else if(side == 1){                  //the atom after CB
	char ach = buff[1];                //atom identifier
	if(ach=='C' || ach=='O' ||
	   ach=='N' || ach=='S'){          //heavyAtom in the side chain(after CB)
	  pt->addAtom(ach, x, y, z);
	}
      }
    }
    
    if(!getline(in,line)){
      list.push_back(pt);
      break;
    }
    string field=line.substr(0,6);
    if(field != AtomField){
      if( origin==1 && !(pt->SizeCheck()) ){
	errMsg += " " + pt->getName(); 
	errStatus = 2;
      }
      else{
	list.push_back(pt);
      }
      break;
    }
  }
  in.close();
  
  //prepare the Contact Map --- CM
  int length = list.size();
  cm = new char*[length];
  for(int i=0; i<length; i++) {
    cm[i]= new char[length];
    for(int j=0; j<length; j++)
      cm[i][j]=0;
  }
  for(int i=0; i<length; i++){
    aa* iaa = list[i];
    atom* ia;
    if(iaa->getName()=="GLY") ia = iaa->getAtomAt(0);
    else ia = iaa->getAtomAt(1);
    for(int j=i+2; j<length; j++){
      aa* jaa = list[j];
      atom* ja;
      if(jaa->getName()=="GLY") ja = jaa->getAtomAt(0);
      else ja = jaa->getAtomAt(1);
      if(ia->getDist(*ja) < cutoff )
	cm[i][j] = cm[j][i] = 1;
    }
  }
}

class pseudoProtein : public protein{ //inheritance
  vector<pseudo_aa*> pseudo_list;
 public:
  pseudoProtein(string iname, int isOrigin=0, double cutoff=7.5)
    :protein(iname, isOrigin, cutoff)
  {
    for(int i=0; i<getLength(); i++){
      pseudo_list.push_back(new pseudo_aa(*getElementAt(i)));
    }
  }
  
  ~pseudoProtein(){
    for(int i=0; i<pseudo_list.size(); i++){
      delete pseudo_list[i];
    }
    pseudo_list.clear();
  }
  
  pseudo_aa* getPseudoElementAt(int i){
    if(i>=0 && i<pseudo_list.size()) return pseudo_list[i];
  }
  
  void pseudoPrint(ostream& out){
    for(int i=0; i<pseudo_list.size(); i++){
      pseudo_list[i]->print(out);
    }
  }  
  
  int getPseudoChi(int i, double& chi){
    if(i<0 || i>=getLength()) return 0;
    pseudo_aa* theAA = getPseudoElementAt(i);
    if(theAA->getName()!="GLY" && 
       theAA->getName()!="ALA" &&
       theAA->getName()!="PRO" &&
       theAA->getValid()){
      vec* N = theAA->getN();
      vec* CA= theAA->getCa();
      vec* CB= theAA->getCb();
      vec* CG= theAA->getSide();
      vec CG_CB = *CB - *CG;
      vec CB_CA = *CA - *CB;
      vec CA_N  = *N  - *CA;
      double mod2 = CB_CA*CB_CA;
      vec r1 = CG_CB - (CG_CB*CB_CA)*CB_CA/mod2;
      vec r2 = CA_N  - (CB_CA*CA_N )*CB_CA/mod2;
      double cos = r1*r2;
      cos /= sqrt((r1*r1)*(r2*r2));
      double coef = (CB_CA^r1)*r2;
      if(coef>0) chi = acos(cos);
      else chi = PI*2 - acos(cos);
      return 1;
    }
    else return 0;
  }

  void printCM(ostream& out, char* type, double cutoff){
    if(strcmp(type,"CA")==0 || strcmp(type,"CB")==0){
      protein::printCM(out, type, cutoff);
    }
    else if(strcmp(type,"Side")==0){
      int len = getLength();
      for(int i=0; i<len; i++){
	pseudo_aa* iaa=getPseudoElementAt(i);
	vec* ia;
	if(iaa->getName()=="GLY")
	  ia = iaa->getCa();
	else if(iaa->getName()=="PRO" || iaa->getName()=="ALA")
	  ia = iaa->getCb();
	else
	  ia = iaa->getSide();
	for(int j=i+2; j<len; j++){
	  pseudo_aa* jaa=getPseudoElementAt(j);
	  vec* ja;
	  if(jaa->getName()=="GLY")
	    ja = jaa->getCa();
	  else if(jaa->getName()=="PRO" || jaa->getName()=="ALA")
	    ja = jaa->getCb();
	  else
	    ja = jaa->getSide();
	  if(ia->getDist(*ja) < cutoff)
	    out << i+1 << "\t" << j+1 << endl;
	}
      }
    }
  }
};
#endif
