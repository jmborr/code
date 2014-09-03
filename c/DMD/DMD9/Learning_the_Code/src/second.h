#ifndef _SECOND_H_
#define _SECOND_H_
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
using namespace std;
/*using the dihedral angles to estimate the secondary structure propensities
  by providing the fy files, the f,y space is devided in to a 2-D matrix with
  bin size of PI/3.
  
  A  G  M  S  m  g  A
  F  L  R  X  n  h  F
  E  K  Q  W  o  i  E
  D  J  P  V  p  j  D
  C  I  O  U  q  k  C
  B  H  N  T  r  l  B
  A  G  M  S  m  g  A

  S   = {A,F,G,L,M,R}
  H   = {O}
  T   = [J,O,P}
  T'  = {j,o,p}
  U   = {M,R}
  U'  = {m,r}

  for jth amino acid, 
  j-1,j,j+1 E H -----> j Helix
  j-1,j,j+1 E S -----> j Strand
  j E T  & j+1 E T  -> j,j+1 Turn I
  j E T' & j+1 E T' -> j,j+1 Turn I'
  j E U  & j+1 E T' -> j,j+1 Turn II
  j E U' & j+1 E T  -> j,j+1 Turn II'
  else  -------------> j Coil
*/

const double inf = 1.0E38;
const double PI1 = atan(1.0)*4.0;

typedef enum{
  nonBin=-1,
    A=0,B,C,D,E,F,
    G,H,I,J,K,L,
    M,N,O,P,Q,R,
    S,T,U,V,W,X,
    m,r,q,p,o,n,
    g,l,k,j,i,h}
bin;

typedef enum{
  nonSec=-1, HELIX, STRAND, TURN1, TURN1P, TURN2, TURN2P, COIL}
seconds;

char second2Char(seconds i){
  switch(i){
  case HELIX:
    return 'H';
    break;
  case STRAND:
    return 'S';
    break;
  case TURN1:
  case TURN1P:
  case TURN2:
  case TURN2P:
    return 'T';
    break;
  default:
    return 'C';
    break;
  }
}

class Dihedrals{
private:
  double phi;
  double psi;
  bin theBin;
public:
  Dihedrals(){
    phi=inf;
    psi=inf;
    theBin=nonBin;
  }

  Dihedrals(double fi, double yi){
    phi = fi-PI1;
    psi = yi-PI1;
    double f = phi;
    if(f>=PI1*5.0/6.0){
      f -= 2*PI1;
    }
    double y = psi;
    if(y>=PI1*5.0/6.0){
      y -= 2*PI1;
    }
    f += PI1*7.0/6.0;
    y += PI1*7.0/6.0;
    int iphi = static_cast<int>(f/(PI1/3.0));
    int ipsi = static_cast<int>(y/(PI1/3.0));
    theBin = static_cast<bin>(iphi*6+ipsi);
  }
  
  void init(double fi, double yi){
    phi = fi-PI1;
    psi = yi-PI1;
    double f = phi;
    if(f>=PI1*5.0/6.0){
      f -= 2*PI1;
    }
    double y = psi;
    if(y>=PI1*5.0/6.0){
      y -= 2*PI1;
    }
    f += PI1*7.0/6.0;
    y += PI1*7.0/6.0;
    int iphi = static_cast<int>(f/(PI1/3.0));
    int ipsi = static_cast<int>(y/(PI1/3.0));
    theBin = static_cast<bin>(iphi*6+ipsi);
  }

  ~Dihedrals(){}
  
  bin getBin(){
    return theBin;
  }

  int isS(){
    if(theBin==A||theBin==F||theBin==G||
       theBin==L||theBin==M||theBin==R) return 1;
    return 0;
  }

  int isH(){
    if(theBin==O) return 1;
    return 0;
  }

  int isT(){
    if(theBin==J||theBin==O||theBin==P)return 1;
    return 0;
  }

  int isTP(){
    if(theBin==j||theBin==o||theBin==p)return 1;
    return 0;
  }

  int isU(){
    if(theBin==M||theBin==R)return 1;
    return 0;
  }
  
  int isUP(){
    if(theBin==m||theBin==r)return 1;
    return 0;
  }

  double getPhi(){return phi;}

  double getPsi(){return psi;}
};

seconds getSeconds(Dihedrals& curr, Dihedrals& prev, Dihedrals& next){
  if(prev.isH()&&curr.isH()&&next.isH()) return HELIX;
  if(prev.isS()&&curr.isS()&&next.isS()) return STRAND;
  if(prev.isT() &&curr.isT() || 
     curr.isT() &&next.isT())  return TURN1;
  if(prev.isTP()&&curr.isTP() || 
     curr.isTP()&&next.isTP()) return TURN1P;
  if(prev.isU() &&curr.isTP() ||
     curr.isU() &&next.isTP()) return TURN2;
  if(prev.isUP()&&curr.isT()  ||
     curr.isUP()&&next.isT())  return TURN2P;
  return COIL;
}

seconds getHeadSecond(Dihedrals& curr, Dihedrals& next, Dihedrals& nnxt){
  if(curr.isH()&&next.isH()&&nnxt.isH()) return HELIX;
  if(curr.isS()&&next.isS()&&nnxt.isS()) return STRAND;
  if(curr.isT() &&next.isT())  return TURN1;
  if(curr.isTP()&&next.isTP()) return TURN1P;
  if(curr.isU() &&next.isTP()) return TURN2;
  if(curr.isUP()&&next.isT()) return TURN2P;
  return COIL;
}

seconds getTailSecond(Dihedrals& curr, Dihedrals& prev, Dihedrals& pprv){
  if(pprv.isH()&&prev.isH()&&curr.isH()) return HELIX;
  if(pprv.isS()&&prev.isS()&&curr.isS()) return STRAND;
  if(prev.isT() &&curr.isT())  return TURN1;
  if(prev.isTP()&&curr.isTP()) return TURN1P;
  if(prev.isU() &&curr.isTP()) return TURN2;
  if(prev.isUP()&&curr.isT())  return TURN2P;
  return COIL;
}
#endif //_SECOND_H_
