/*  g++ -Wno-deprecated -O3 -o zrrmsd.x zrrmsd.cpp  */
#include <iostream.h>
#include <cstring>
#include <cmath>
#include "readAndWrite.h"
#include "allocate.h"
#include "zrrmsd.h"

#define MAX_N_RES 2000

void abort(){
  printf("  zrrmsd.x [options] model.pdb native.pdb\n");
  printf("  Returns the Z-rRMSD value of CA's between two structures\n");
  printf("  Options:\n");
  printf("   -a outf: pdb file with superposition of model to native\n");
  exit(1);
}

int main(int argc, char **argv){

  char *outf  =NULL;    /*file name to output alignment*/
  char *model =NULL;   /*file name of model pdb       */
  char *native=NULL;  /*file name of native pdb      */

  int nn;                                    /*sequence length*/
  double **mr=alloc2<double>(MAX_N_RES,3);  /*model coordinates*/
  double **nr=alloc2<double>(MAX_N_RES,3); /*native coordinates*/

  double *w=new double[MAX_N_RES];     /*weights when calculating rsmd*/
  double *tmr=new double[MAX_N_RES];  /*traslation of model to origin of coords.*/
  double *tnr=new double[MAX_N_RES]; /*traslation of native to origin of coords.*/
  double **U=alloc2<double>(3,3);   /*Rotation matrix*/

  {/*BLOCK TO READ COMMAND LINE ARGUMENTS*/
    if(argc<3) abort();

    for(int i=1; i<argc; i++){
      if(argv[i][0]=='-'){
	switch(argv[i][1]){
	case 'a': outf=argv[++i]; break;
	case 'h': abort();
	default: printf("Unknown option: %c\n", argv[i][1]); return -1;
	}
      }      
      else{
	if(!model) model=argv[i];
	else native=argv[i];
      }
    }

    if(!model) cerr<<"ERROR: pass a model pdb\n";
    if(!native)cerr<<"ERROR: pass a native pdb\n";
    /*printf("%s %s\n",model,native); exit(1);*/
  }/*end of block ro read input*/



  {/*BLOCK TO READ MODEL AND NATIVE COORDINATES INTO MEMORY*/
    nn=read_pdb(native,nr); /*printf("nn=%d\n",nn);exit(1);*/
    if( nn!=read_pdb(model,mr) ){
      cerr<<"ERROR: native and model have different number of CA's\n";
      return -1;
    }
  }/*end of block to read model and native coordinates into memory*/



  {/*BLOCK TO COMPUTE Z-rRMSD*/
    for(int i=0; i<MAX_N_RES; i++) w[i]=1.0;
    cout<<get_Zrot(mr,tmr,nr,tnr,nn,w,U)<<endl;
    if(outf){
      cerr<<"ERROR:writing output alignment not implemented yet!\n";
    }
  }/*end of block to read compute Z-rRMSD*/

  return 0;
}
