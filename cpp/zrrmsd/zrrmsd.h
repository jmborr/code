#include "rmsd.h"

/*======================================================================
  compute Z-rRMSD */
double zrrmsd(double **x,double **y,long nn,double *w){
  double coef=0.42-0.05*nn*exp(-nn/4.7)+0.63*exp(-nn/37);
  double invd=1/(0.09+1.16*exp(-nn/1.6)+0.25*exp(-nn/36));
  double rmsd=getrmsd(x,y,nn,w);
  
  double rg1=get_rg(x,nn);
  double rg2=get_rg(y,nn);
  
  return (rmsd/sqrt(rg1*rg1+rg2*rg2-2*coef*rg1*rg2 -1))*invd;
}
/*======================================================================
  compute Z-rRMSD 
  store in tx the translation such that weighted CM of x would be in
  the origin (same for ty). We pass tx,ty,U so that we don't have to
  allocate these variables*/
double get_Zrot(double **x,double *tx,double **y,double *ty,int nn,double *w,double **U ){
  double coef=0.42-0.05*nn*exp(-nn/4.7)+0.63*exp(-nn/37);
  double invd=1/(0.09+1.16*exp(-nn/1.6)+0.25*exp(-nn/36));
  double rmsd=get_rot(x,tx,y,ty,nn,w,U);  /*printf("rmsd=%lf\n",rmsd);exit(1);*/
  
  double rg1=get_rg(x,nn);
  double rg2=get_rg(y,nn);
  
  return (rmsd/sqrt(rg1*rg1+rg2*rg2-2*coef*rg1*rg2) -1)*invd;
}
/*======================================================================*/
