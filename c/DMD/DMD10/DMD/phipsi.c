#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "phipsi.h"
#include "bcp.h"

typedef struct{
  atom* N;
  atom* C;
  atom* Ca;
}backbone_aa;


static atom* a;
static int natom;
static int naa;
static backbone_aa* protein;
static good;
static char file_name[80];
static FILE* phipsi_out;
static char list_file[80];
static FILE* list_in;
static double PI2;
static double delta_t;

void init_phipsi(){
  int i, in, ica, ic;
  good=0;
  PI2 = atan(1.0)*8.0;
  printf("We do not calculate phi-psi values\n");
  if(yes()) return;
  
  printf("What is the phi-psi storing file?\n");
  scanf("%s",file_name);
  if((phipsi_out=fopen(file_name,"wb"))==NULL) return;

  
  a=get_atom();
  natom = get_atom_number();
  printf("What is the protein list file?\n");
  printf("FORMAT: N_AminoAcids\nN_1 Ca_1 C_1\nN_2 Ca_2 C_2\n.....\n");
  scanf("%s", list_file);
  if((list_in=fopen(list_file,"r"))==NULL) return;
  
  
  fscanf(list_in,"%ld", &naa);
  protein = (backbone_aa*)malloc(naa*sizeof(backbone_aa));
  for(i=0; i<naa; i++){
    fscanf(list_in, "%ld %ld %ld", &in, &ica, &ic);
    if(in>natom) return;
    protein[i].N = a+(in-1);
    
    if(ica>natom) return;
    protein[i].Ca= a+(ica-1);
    
    if(ic>natom) return;
    protein[i].C = a+(ic-1);
  }
  fclose(list_in);
  
  printf("What is time interval to record the phi-psi values\n");
  scanf("%lf", &delta_t);
  
  good=1;
  return;
}


void get_phipsi(backbone_aa* prev, 
		backbone_aa* curr,
		backbone_aa* next,
		double* phi, double* psi)
{
  crd N2_C1, C1_Ca1, Ca1_N1, N1_C0;
  crd r1, r2, r3, r4;
  double mod2, dot, coef, cos;
  
  N2_C1.x = curr->C->q.x - next->N->q.x;
  N2_C1.y = curr->C->q.y - next->N->q.y;
  N2_C1.z = curr->C->q.z - next->N->q.z;
  
  C1_Ca1.x = curr->Ca->q.x - curr->C->q.x;
  C1_Ca1.y = curr->Ca->q.y - curr->C->q.y;
  C1_Ca1.z = curr->Ca->q.z - curr->C->q.z;
  
  Ca1_N1.x = curr->N->q.x - curr->Ca->q.x;
  Ca1_N1.y = curr->N->q.y - curr->Ca->q.y;
  Ca1_N1.z = curr->N->q.z - curr->Ca->q.z;
  
  N1_C0.x = prev->C->q.x - curr->N->q.x;
  N1_C0.y = prev->C->q.y - curr->N->q.y;
  N1_C0.z = prev->C->q.z - curr->N->q.z;

  /*phi*/
  mod2 = Ca1_N1.x*Ca1_N1.x + Ca1_N1.y*Ca1_N1.y + Ca1_N1.z*Ca1_N1.z;
  
  dot = C1_Ca1.x*Ca1_N1.x + C1_Ca1.y*Ca1_N1.y + C1_Ca1.z*Ca1_N1.z;
  coef = dot/mod2;
  r1.x = C1_Ca1.x - Ca1_N1.x*coef;
  r1.y = C1_Ca1.y - Ca1_N1.y*coef;
  r1.z = C1_Ca1.z - Ca1_N1.z*coef;
  
  dot = N1_C0.x*Ca1_N1.x + N1_C0.y*Ca1_N1.y + N1_C0.z*Ca1_N1.z;
  coef = dot/mod2;
  r2.x = N1_C0.x - Ca1_N1.x*coef;
  r2.y = N1_C0.y - Ca1_N1.y*coef;
  r2.z = N1_C0.z - Ca1_N1.z*coef;
  
  cos = (r1.x*r2.x+r1.y*r2.y+r1.z*r2.z);
  cos/= sqrt((r1.x*r1.x+r1.y*r1.y+r1.z*r1.z)*(r2.x*r2.x+r2.y*r2.y+r2.z*r2.z));

  coef = r2.x*(Ca1_N1.y*r1.z-Ca1_N1.z*r1.y);
  coef+= r2.y*(Ca1_N1.z*r1.x-Ca1_N1.x*r1.z);
  coef+= r2.z*(Ca1_N1.x*r1.y-Ca1_N1.y*r1.x);

  if(coef >0) *phi = acos(cos);
  else *phi = PI2 - acos(cos);
  
  /*psi*/
  mod2 = C1_Ca1.x*C1_Ca1.x + C1_Ca1.y*C1_Ca1.y + C1_Ca1.z*C1_Ca1.z;
  
  dot  = Ca1_N1.x*C1_Ca1.x + Ca1_N1.y*C1_Ca1.y + Ca1_N1.z*C1_Ca1.z;
  coef = dot/mod2;
  r3.x = Ca1_N1.x - C1_Ca1.x*coef;
  r3.y = Ca1_N1.y - C1_Ca1.y*coef;
  r3.z = Ca1_N1.z - C1_Ca1.z*coef;
  
  dot  = N2_C1.x*C1_Ca1.x + N2_C1.y*C1_Ca1.y + N2_C1.z*C1_Ca1.z;
  coef = dot/mod2;
  r4.x = N2_C1.x - C1_Ca1.x*coef;
  r4.y = N2_C1.y - C1_Ca1.y*coef;
  r4.z = N2_C1.z - C1_Ca1.z*coef; 

  cos = r3.x*r4.x+r3.y*r4.y+r3.z*r4.z;
  cos/= sqrt((r3.x*r3.x+r3.y*r3.y+r3.z*r3.z)*(r4.x*r4.x+r4.y*r4.y+r4.z*r4.z));

  coef = r4.x*(C1_Ca1.y*r3.z-C1_Ca1.z*r3.y);
  coef+= r4.y*(C1_Ca1.z*r3.x-C1_Ca1.x*r3.z);
  coef+= r4.z*(C1_Ca1.x*r3.y-C1_Ca1.y*r3.x);

  if(coef < 0) *psi = acos(cos);
  else *psi = PI2 - acos(cos);
}

void write_phipsi(){
  double phi, psi;
  int i, nbyte;
  unsigned char s[512];
  moveatoms();
  for(i=1; i<naa-1; i++){
    get_phipsi(&protein[i-1], &protein[i], &protein[i+1],
	       &phi, &psi);
    nbyte=sprintf(s,"%ld %lf %lf\n", i+1, phi, psi);
    if(nbyte<=0){
      fclose(phipsi_out); 
      good=0; return;
    }
    if(fwrite(s, 1, nbyte, phipsi_out)!=nbyte){
      fclose(phipsi_out); 
      good=0; return;
    }
  }
}

void process_phipsi(){
  static double last_time=0;
  if(good){
    double new_time= get_time();
    if((new_time - last_time) > delta_t){
      last_time = new_time;
      write_phipsi();
    }
  }
}
