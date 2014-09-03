#include<iostream.h>
#include<math.h>
#include "rmsd.h"
#include "readAndWrite.h"
#include "allocate.h"

#define MAX_N_RES 2000

void abort(){
  printf("  Rg.x [options]\n");
  printf("  Returns the Radius of gyration\n");
  printf("  Options:\n");
  printf("   -a: pdbfile\n");
  printf("   -b: Returns also expected Rg according to Betancourt\n       & Skolnick, Biopolymers 59, 305-309 (2001)\n");
  printf("   -c: trajectory file (will return filtered trajectory\n       file, with suffix \".filtered\"\n");
  exit(1);
}

bool _EXPECTED_=false;

int main(int argc, char **argv){

  char *pdbf=NULL;
  char *tra=NULL;
  int nn; /*sequence length*/
  double **r=alloc2<double>(MAX_N_RES,3); /*coordinates*/

  {/*BLOCK TO READ COMMAND LINE ARGUMENTS*/
    /*printf("argc=%d\n",argc);*/
    if(argc<2) abort();

    for(int i=1; i<argc; i++){
      if(argv[i][0]=='-'){
	switch(argv[i][1]){
	case 'a': pdbf=argv[++i]; break;
	case 'b': _EXPECTED_=true; break;
	case 'c': tra=argv[++i]; /*printf("%s\n",tra);*/ break;
	case 'h': abort(); break;
	default: printf("Unknown option: %c\n", argv[i][1]); return -1;
	}
      }      
    }
  }/*end of block ro read input*/


  if(pdbf){
    {/*BLOCK TO READ PDB FILE COORDINATES INTO MEMORY*/
      nn=read_pdb(pdbf,r); /*printf("nn=%d\n",nn);exit(1);*/
    }/*end of block to read model and native coordinates into memory*/

    {/*BLOCK TO COMPUTE RG*/
      double erg;
      printf("%7.2lf",get_rg(r,nn));
      if(_EXPECTED_){
	printf(" %7.2lf\n",erg=3.08*pow(nn,0.3333));
      }
      else 
	printf("\n");
    }/*end of block to read compute rmsd*/
  }

  if(tra){
    int i;
    int all=0;
    int accepted=0;
    double rgco,rg;
    FILE *ptra=fopen(tra,"r");
    char *outf=new char[512];
    sprintf(outf,"%s.filtered",tra);
    FILE *poutf=fopen(outf,"w");
    char *header=new char[512];
    nn=read_trajectory_snapshot(ptra,r,header);
    rgco=1.20*2.2*pow(nn,0.38); /*printf("rgco=%lf\n",rgco);*/
    while(nn){
      all++;
      rg=get_rg(r,nn); /*printf("rg=%lf\n",rg);*/
      if(rg<rgco){ /*write trajectory*/
	accepted++;
	fprintf(poutf, "%s",header);
	for(i=0;i<nn;i++) fprintf(poutf,"%10.3lf%10.3lf%10.3lf\n",r[i][0],r[i][1],r[i][2]);
      }
      nn=read_trajectory_snapshot(ptra,r,header);
      /*printf("nn=%d\n",nn);*/
    }
    printf("%d snapshots accepted out of %d\n",accepted,all);
    fclose(ptra);
    fclose(poutf);
  }

  return 0;
}
