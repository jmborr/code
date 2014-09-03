#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cell_size.h"
#include "bonds.h"

shell_order_scheme soc;
atom *a; /*initialized first time when "allocate_shell_scheme" is invoked*/
static dimensions *bound; /*initialized first time when set_amso is invoked*/

/*=======================================================================*/
void set_input_cell_size(double cs){ 
  soc.input_cell_size=cs;  /*printf("cell_size=%lf\n",cell_size);*/
}
/*=======================================================================*/
void dump_input_cell_size(void){
  printf("CELL SIZE=%lf\n",soc.input_cell_size);
}
/*=======================================================================*/
double return_input_cell_size(void){
  return soc.input_cell_size;
}
/*=======================================================================*/
void set_input_safe_limit(double sf){ 
  soc.sf=sf;  /*printf("safe_limit=%lf\n",sf);*/
}
/*=======================================================================*/
void dump_input_safe_limit(void){
  printf("SAFE LIMIT=%lf\n",soc.sf);
}
/*=======================================================================*/
void set_dl(double dl){
  soc.dl=dl;
}
/*=======================================================================*/
double return_input_cell(void){
  return soc.input_cell_size;
}
/*=======================================================================*/
void set_amso(double interaction_range)
{/*we assume that the whole system is a cube, thus bound[0].dl==bound[1].dl==bound[2].dl*/
  int i;
  bound=get_bounds();
  soc.amso = (int)(1+ interaction_range/bound[0].dl);/*printf("amso=%d\n",soc.amso);*/
}
/*=======================================================================*/

int allocate_and_init_shell_order_scheme(double ***coldata,atom *sam,int nat){
  int nat2,i,j,l,tot;
  shell_order so;
  double range;
  short **temp,*temp2;
  /*printf("allocate_and_init_shell_order_scheme(...)\n");*/
  a=get_atom();/*this "a" is defined at the begining of this file. This way we
		 initialize this variable whose scope is global within this
		 file. We can use later "a" on any function defined in this
		 file without having to pass it as an argument*/
  /*"bound" has already been allocated when set_amso invoked*/
  soc.nat=nat;
  if(!(temp=(short **)malloc(nat*sizeof(short *)))) return 0;
  temp--;
  l=nat*soc.amso;
  if(!(temp[1]=(short *)malloc(l*sizeof(short)))) return 0;
  for(i=0;i<l;i++)
    temp[1][i]=0;


  if(!(temp2=(short *)malloc(nat*sizeof(short)))) return 0;
  temp2--;

  if(!(soc.mso_tp=(shell_order*)malloc(nat*sizeof(shell_order)))) return 0;
  soc.mso_tp--;/*this way soc.mso_tp[1] corresponds to first memory holder*/
  if(!(soc.mso_tptq=(shell_order**)malloc(nat*sizeof(shell_order*)))) return 0;
  soc.mso_tptq--;
  nat2=nat*nat;
  if(!(soc.mso_tptq[1]=(shell_order*)malloc(nat2*sizeof(shell_order))))
    return 0;
  soc.mso_tptq[1]--;

  for(i=2;i<=nat;i++){
    temp[i]=temp[i-1]+soc.amso;
    soc.mso_tptq[i]=soc.mso_tptq[i-1]+nat;
  }  

  for(i=1;i<=nat;i++){
    temp2[i]=0;
    soc.mso_tp[i]=0;
  }

  tot=soc.nat;
  for(i=1;i<=nat;i++){
    for(j=1;j<=i;j++){/*require i>=j for coldata[i][j]*/
      if(coldata[i][j]){
	l=(int)(coldata[i][j][0]);
	if(l==3) 
	  range=coldata[i][j][1];    /*ellastic collision  */
	else 
	  range=coldata[i][j][l-2]; /*inellastic collision*/
      }
      else
	range=sam[i].s+sam[j].s;
      so=(int)(1+range/soc.dl);
      soc.mso_tptq[i][j]=so;
      soc.mso_tptq[j][i]=so;
      if(soc.mso_tp[i]<so) soc.mso_tp[i]=so;
      if(soc.mso_tp[j]<so) soc.mso_tp[j]=so;
      if(!temp[i][so]){
	temp[i][so]=1;
	temp2[i]++;
	tot++;
      }
    }
  }

  if(!(soc.boundaries=(int **)malloc(nat*sizeof(int*)))) return 0;
  soc.boundaries--;
  if(!(soc.boundaries[1]=(int*)malloc(tot*sizeof(int)))) return 0;
  tot=0;
  for(i=1;i<=nat;i++){
    if(i>1)
      soc.boundaries[i]=soc.boundaries[i-1]+l;
    soc.boundaries[i][0]=temp2[i]; 
    tot++;
    l=1;
    for(j=1;j<=soc.amso;j++){
      if(temp[i][j]){
	soc.boundaries[i][l]=j;
	l++;
	tot++;
      }
    }
  }

  free(temp[1]);
  temp++; free(temp);
  temp2++; free(temp2);
  /*printf("End of allocate_and_init_shell_order_scheme(...)\n");*/
  return 1;
}
/*=======================================================================*/
shell_order return_mso (atom_number p, atom_number q){
  if(isFriend(p,q))/*p and q are bonded*/
    return INF_SO; /*by definition, inifinite shell order bonded atoms*/
  return soc.mso_tptq[a[p].origc][a[q].origc];
}
/*=======================================================================*/
shell_order return_so (atom_number p, atom_number q){
  atom *ap,*aq;
  cell_index di;
  shell_order so,amso=soc.amso;

  ap=a+p;
  aq=a+q;

  so=abs(ap->i.x.i-aq->i.x.i);
  if(so>amso) so-=bound[0].period;

  di=abs(ap->i.y.i-aq->i.y.i); 
  if(di>amso) di-=bound[1].period;
  if(di>so) so=di;

  di=abs(ap->i.z.i-aq->i.z.i); 
  if(di>amso) di-=bound[2].period;
  if(di>so) so=di;

  return so;
}
/*=======================================================================*/
void dump_boundaries(void){
  int i,j,nb,**a=soc.boundaries,*b;
  for(i=1;i<=soc.nat;i++){
    printf("type=%d :",i);
    b=a[i];
    nb=b[0];
    for(j=1;j<=nb;j++)
      printf("%d ",b[j]);
    printf("\n");
  }
}
/*=======================================================================*/
