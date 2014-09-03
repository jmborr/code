#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include "second.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc<4){
    cout << "usage: fy2prop.linux fy_File nRes out_File [startFrame nFrame]" << endl;
    exit(1);
  }
  
  ifstream in(argv[1]);
  int nRes = atoi(argv[2]);
  ofstream out(argv[3]);
  int startFrame=0, nFrame=0;
  if(argc==6){
    startFrame=atoi(argv[4]);
    nFrame=atoi(argv[5]);
  }


  Dihedrals FYs[nRes-2];
  seconds props[nRes-2];
  double seqProp[nRes-2][7];
  for(int i=0; i<nRes-2; i++)
    for(int j=0; j<7; j++)seqProp[i][j]=0;
  int ires;
  double f,y;
  int end=0;
  int iCount=0;
  int iFrame=0;
  int getSecond=0;
  while(in >> ires >> f >> y){
    if(ires!=2){
      cout << "mismatch!" << endl;
      exit(1);
    }
    FYs[0].init(f,y);
    for(int i=3; i<=nRes-1; i++){
      if(in>>ires>>f>>y){
	FYs[i-2].init(f,y);
      }
      else{
	end=1;
	break;
      }
    }
    if(end)break;
    iCount++;

    props[0]=getHeadSecond(FYs[0],FYs[1],FYs[2]);
    props[nRes-3]=getTailSecond(FYs[nRes-3],FYs[nRes-4],FYs[nRes-5]);
    for(int i=1; i<nRes-3; i++){
      props[i]=getSeconds(FYs[i], FYs[i-1], FYs[i+1]);
    }
    /*check HELIX*/
    int pt=0;
    int start=0;
    while(pt<nRes-2){
      if(props[pt]==HELIX){
	if(!start){//begin
	  start=1;
	  int j=1;
	  while(static_cast<seconds>(props[pt-j])==TURN1 && pt>j){
	    props[pt-j]=HELIX;
	    j++;
	  }
	}
      }
      else{
	if(start){//end of helix
	  start=0;
	  int j=0;
	  while(static_cast<seconds>(props[pt+j])==TURN1 && pt+j<nRes-2){
	    props[pt+j]=HELIX;
	    j++;
	  }
	  pt+=j;
	}
      }
      pt++;
    }
    /*check STRAND*/
    pt=0;
    start=0;
    while(pt<nRes-2){
      if(props[pt]==STRAND){
	if(!start){//begin
	  start=1;
	  if(pt>0 && props[pt-1]==COIL)props[pt-1]=STRAND;
	}
      }
      else{
	if(start){//end of strand
	  start=0;
	  if(props[pt]==COIL)props[pt]=STRAND;
	}
      }
      pt++;
    }
    
    if(!startFrame&&!nFrame){
      for(int i=0; i<nRes-2;i++)
	seqProp[i][props[i]]+=1;
    }
    else if(iCount>=startFrame && iFrame<nFrame){
	iFrame++;
	for(int i=0; i<nRes-2;i++)
	  seqProp[i][props[i]]+=1;
    }
  }
  
  out << "# iRes Helix Strand Turn Coil" << endl;
  for(int i=0; i<nRes-2; i++)
    for(int j=0; j<7; j++){
      if(!startFrame&&!nFrame){
	seqProp[i][j]/=static_cast<double>(iCount);
      }
      else{
	seqProp[i][j]/=static_cast<double>(iFrame);
      }
    }
  
  if(!startFrame&&!nFrame){cout<< iCount << endl;}

  for(int i=0; i<nRes-2; i++)
    out<< i+2 << " " << seqProp[i][0] << endl;
  out << "&" << endl;
  
  for(int i=0; i<nRes-2; i++)
    out<< i+2 << " " << seqProp[i][1] << endl;
  out << "&" << endl;
  
  for(int i=0; i<nRes-2; i++)
    out<< i+2 << " " << seqProp[i][2]+ seqProp[i][3]+ seqProp[i][4]+ seqProp[i][5] << endl;
  out << "&" << endl;
  
  for(int i=0; i<nRes-2; i++)
    out<< i+2 << " " << seqProp[i][6] << endl;
  out << "&" << endl;
  
  in.close();
  out.close();
  }
