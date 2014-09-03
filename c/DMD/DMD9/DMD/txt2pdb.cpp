#include <fstream>
#include <iostream>
#include <string>
#include "pseudoPep.h"
#include "pseudoTXT.h"
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

int main(int argc, char* argv[]){
  if(argc<6){
    cout << "usage: txt2pdb.linux txt inputMethod inputFile nProt pdb" << endl;
    exit(1);
  }
  pseudoPep p(argv[3], static_cast<input_t>(atoi(argv[2])));
  int natoms = p.makeIndex();
  int nProts = atoi(argv[4]);
  
  double** r = (double**)malloc(natoms*nProts*sizeof(double*));
  r[0]=(double*)malloc(3*natoms*nProts*sizeof(double));
  for(int i=1; i<natoms*nProts; i++) r[i]=r[i-1]+3;
  
  /*read txt*/
  ifstream in(argv[1]);
  string line;
  int theKey=0;
  int tmpKey=0;
  int nTotal=0;
  int index=0;
  int iprot=0;
  char buf[1024];
  int iatom, itype;
  double vx,vy,vz,boundary;
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
            str2ch(line, buf);
            sscanf(buf, "%ld", &nTotal);//read nTotal
	    if(nTotal<(natoms*nProts)){
	      cout << "Input $nProts is not bigger than that in TXT file" << endl;
	      exit(1);
	    }
            break;
          case 3://type of atoms
            break;
          case 4: //non-elastic collisions
          case 5: //elastic collisions
          case 6: //bond
          case 7: //reaction
            break;
          case 8://list of atoms
	    if(iprot<nProts){
	      str2ch(line, buf);
	      sscanf(buf, "%ld %ld %lf %lf %lf %lf %lf %lf",
		     &iatom, &itype, &r[index][0], &r[index][1], &r[index][2], 
		     &vx, &vy, &vz);
	      index++;
	      if(index==natoms*(iprot+1)){
		iprot++;
	      }
	    }
            break;
          case 9:
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
  /*considering PBC*/
  for(int i=1; i<natoms*nProts; i++){
    if(r[i][0]-r[i-1][0] > boundary/2.0f ) r[i][0]-=boundary;
    if(r[i][0]-r[i-1][0] < -boundary/2.0f) r[i][0]+=boundary;
    if(r[i][1]-r[i-1][1] > boundary/2.0f ) r[i][1]-=boundary;
    if(r[i][1]-r[i-1][1] < -boundary/2.0f) r[i][1]+=boundary;
    if(r[i][2]-r[i-1][2] > boundary/2.0f ) r[i][2]-=boundary;
    if(r[i][2]-r[i-1][2] < -boundary/2.0f) r[i][2]+=boundary;
  }
  double range = 10;
  double cm[nProts][3];
  int& nres = natoms;
  for(int ip=0; ip<nProts; ip++){
    cm[ip][0]=r[ip*nres][0];
    cm[ip][1]=r[ip*nres][1];
    cm[ip][2]=r[ip*nres][2];
    for(int i=1+ip*nres; i<nres*(ip+1); i++){
      if(r[i][0]-r[i-1][0] > boundary/2.0f  && r[i-1][0]+boundary-r[i][0]<range) r[i][0]-=boundary;
      if(r[i][0]-r[i-1][0] < -boundary/2.0f && r[i][0]+boundary-r[i-1][0]<range) r[i][0]+=boundary;
      if(r[i][1]-r[i-1][1] > boundary/2.0f  && r[i-1][1]+boundary-r[i][1]<range) r[i][1]-=boundary;
      if(r[i][1]-r[i-1][1] < -boundary/2.0f && r[i][1]+boundary-r[i-1][1]<range) r[i][1]+=boundary;
      if(r[i][2]-r[i-1][2] > boundary/2.0f  && r[i-1][2]+boundary-r[i][2]<range) r[i][2]-=boundary;
      if(r[i][2]-r[i-1][2] < -boundary/2.0f && r[i][2]+boundary-r[i-1][2]<range) r[i][2]+=boundary;
      cm[ip][0]+=r[i][0];
      cm[ip][1]+=r[i][1];
      cm[ip][2]+=r[i][2];
    } 
    for(int i=0; i<3; i++)cm[ip][i]/=(double)nres;
  }
  
  int alloc[nProts];
  for(int ip=0; ip<nProts; ip++)alloc[ip]=1;
  int next=0;
  alloc[next]=0;
  double cmCurr[3];

  for(int ncount=1; ncount<nProts; ncount++){
    for(int j=0; j<3; j++) cmCurr[j]=0;
    int ialloc=0;
    for(int i=0; i<nProts; i++)if(!alloc[i]){
      ialloc++;
      for(int j=0;j<3;j++)cmCurr[j]+=cm[i][j];
    }
    for(int j=0; j<3; j++) cmCurr[j]/=(double)ialloc;
    
    double dist=0, min=3000*pow(boundary,2);
    int imin;
    for(int i=0; i<nProts; i++)if(alloc[i]){
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

  {
    ofstream out(argv[5]);
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
      p.printPDB(out);
      if(iProts<nProts-1) out<< "TER" << endl;
    }
    out << "END" << endl;
  }
}
