#include<iostream>
#include <math.h>
#include "pdbClasses2.h"
#include "readAndWrite.h"
#include "allocate.h"

#define MAX_N_RES 2000

void abort(){
  printf("  checkCACAdis.x [options]\n");
  printf("  Returns all CA_i-CA_{i+1} distances\n");
  printf("  Options:\n");
  printf("   -a: pdbfile\n");
  exit(1);
}

int main(int argc, char **argv){

  int i;
  char *pdbf=NULL;
  int nn; /*sequence length*/
  double **r=alloc2<double>(MAX_N_RES,3); /*coordinates*/
  double dx,dy,dz;

  if(argc<2) abort();

  for(int i=1; i<argc; i++){
    if(argv[i][0]=='-'){
      switch(argv[i][1]){
      case 'a': pdbf=argv[++i]; break;
      case 'h': abort(); break;
      default: printf("Unknown option: %c\n", argv[i][1]); return -1;
      }
    }      
  }

  nn=read_pdb(pdbf,r); /*printf("nn=%d\n",nn);exit(1);*/
  
  for(i=0;i<nn-1;i++){
    dx=r[i+1][0]-r[i][0] ; dy=r[i+1][1]-r[i][1] ; dz=r[i+1][2]-r[i][2] ;
    printf("%5.2lf\n",sqrt(dx*dx+dy*dy+dz*dz));
  }

    
  return 0;
}
