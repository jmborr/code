#include<iostream>
#include "manageSeqDat.h"
using namespace std;

int main(){
  char *file=new char[200];
  int *beginSec=new int[1000];
  int *resperSec=new int[1000];
  int nSec=0;
  int i,curr,confidence;
  char *aa;
  aa=new char[4];
  sprintf(file,"./10gsA.SEQ");
  printf("%s\n",file);
  nSec=initSecSeg(file,beginSec,resperSec);
  printf("nSec=%d\n",nSec);
  for(i=0;i<nSec;i++){
    printf("%d %d %d\n",i,beginSec[i],resperSec[i]);
  }

  return 0;
}
