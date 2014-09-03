#ifndef _MY_RMSD_
#define _MY_RMSD_


#include<iostream>
#include <math.h>
#include "allocate.h"

#define zero 0.0

void shift_CM(double **x, double *tx, int nn){
  int i,j;
  double cm;
  for(j=0;j<3;j++){
    cm=0;
    for(i=0;i<nn;i++) cm+=x[i][j]; 
    cm/=nn;
    tx[j]=-cm;
  }
}
/*=========================================================*/
void shift_wCM(double **x, double *tx, int nn, double *w){
  int i,j;
  double cm,W=0;
  for(i=0;i<nn;i++) W+=w[i]; /*printf("W=%lf\n",W);*/
  for(j=0;j<3;j++){
    cm=0;
    for(i=0;i<nn;i++) cm+=w[i]*x[i][j]; 
    cm/=W; 
    tx[j]=-cm; /*printf("cm=%lf,tx[%d]=%lf\n",cm,j,tx[j]);*/
  }
}
/*======================================================================*/
double get_rg(double **x, int nn){
  int i,j;
  double cm[3],d,rg=0;
  for(j=0;j<3;j++){/*calculate center of mass first*/
    cm[j]=0;
    for(i=0;i<nn;i++) cm[j]+=x[i][j];
    cm[j]/=nn;
  }
  for(i=0;i<nn;i++)   /*sqrt( sum_i( (r_i - cm)^2 ) )*/
    for(j=0;j<3;j++){
      d=x[i][j]-cm[j];
      rg+=d*d;
    }
  return sqrt(rg/nn);
}
/*======================================================================*/
double getrmsd(double **x,double **y,int nn,double *w=NULL){
  double r[3][3], rr[3][3], mu[3], p[4],q[4];
  double aa[3][3], b[3][3], u[3][3];
  double rij,cx,cy;
  int i,j,n,m,k;
  int nnn=0;
  double sum=0;
  double detR ;
  double *tx=new double[3],*ty=new double[3];
  
  if(w){ for(k=0;k<nn;k++) if(w[k]>zero) nnn++; }
  else{
    w=new double[nn]; for(i=0;i<nn;i++) w[i]=1.0;
    nnn=nn;
  }

  shift_wCM(x,tx,nn,w); shift_wCM(y,ty,nn,w);

  for(i=0;i<3;i++) /*obtain matrix r[][]*/
    for(j=0;j<3;j++){
      rij=0;
      for(k=0;k<nn;k++)	rij+=w[k]*(y[k][i]+ty[i])*(x[k][j]+tx[j]);
      r[i][j]=rij;
    }
  
  /*find determinant of r[][]*/
  detR =r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1]);
  detR-=r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0]);
  detR+=r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]);

  for(i=0;i<3;i++) /*find matrix rr*/
    for(j=0;j<3;j++){
      rij=0;
      for(k=0;k<3;k++) rij+=r[k][i]*r[k][j];
      rr[i][j]=rij;
    }

  /* characteristic third order polinomial p(mu)= det|rr -mu^2I| */  
  p[0]=1;
  p[1]=-(rr[0][0]+rr[1][1]+rr[2][2]);
  p[2] =rr[1][1]*rr[2][2]-rr[1][2]*rr[2][1];
  p[2]+=rr[2][2]*rr[0][0]-rr[2][0]*rr[0][2];
  p[2]+=rr[1][1]*rr[0][0]-rr[1][0]*rr[0][1];
  p[3] =rr[0][0]*(rr[1][1]*rr[2][2]-rr[1][2]*rr[2][1]);
  p[3]-=rr[0][1]*(rr[1][0]*rr[2][2]-rr[1][2]*rr[2][0]);
  p[3]+=rr[0][2]*(rr[1][0]*rr[2][1]-rr[1][1]*rr[2][0]);
  p[3]=-p[3];
  /*  for(i=0;i<=3;i++)
    printf("%lf\n",p[i]);*/
  
  {
    double z=-p[1];
    for(i=2;i<4;i++)
      if(z<fabs(p[i]))z=fabs(p[i]);
    k=3;
    /* solution of the characteristic equation for eigenvalues */
    for(i=0;i<2;i++){
      double q1;
      double dz;
      do{
	q[0]=p[0];
	q1=p[0];
	for(j=1;j<k;j++){
	  q[j]=q[j-1]*z+p[j];
	  q1=q1*z+q[j];
	}
	q[k]=q[k-1]*z+p[k];
	dz=q[k]/q1;
	z=z-dz;
      }while(dz>1.0e-14*z);
      p[k]=z;
      for(j=0;j<k;j++)
	p[j]=q[j];
      k--;
    }
    p[1]=-p[1]/p[0];
  }

  for(k=0;k<nn;k++)
    for(i=0;i<3;i++){
      cx=static_cast<double>(x[k][i]+tx[i]);
      cy=static_cast<double>(y[k][i]+ty[i]);
      sum+=w[k]*(cx*cx+cy*cy);
    }
  
  sum*=0.5;
  
  for(i=0;i<3;i++){  mu[i]=sqrt(p[i+1]);  sum-=mu[i];  }   

  /*find minimum of mu[] and change sign, to avoid improper rotation*/
  i=1;
  if(detR<0){
    for(j=2;j<=3;j++)
      if(p[j]<p[i])i=j;
    sum+=2*mu[i-1];
    mu[i-1]=-mu[i-1];  
  }

  /*free memory*/
  free(tx); free(ty);
  return sqrt(fabs(sum*2/nnn)); /*divide by number of non-zero weigths*/
}
/*=======================================================================*/
/* store in tx the translation such that weighted CM of x would be in
   the origin (same for ty). We pass tx,ty,U so that we don't have to
   allocate these variables*/
double get_rot(double **x,double *tx,double **y,double *ty,int nn,double *w,double **U ){
  double r[3][3], rr[3][3], mu[3], p[4],q[4];
  double aa[3][3], b[3][3], u[3][3];
  double rij,cx,cy;
  int i,j,n,m,k,nnn=0;
  double sum=0;
  double detR ;

  /*for(int i=0;i<nn;i++)  printf("%lf %lf %lf\n",x[i][0],x[i][1],x[i][2]);
    for(int i=0;i<nn;i++)  printf("%lf %lf %lf\n",y[i][0],y[i][1],y[i][2]); */
  
  shift_wCM(x,tx,nn,w); /*printf("tx=[%lf,%lf,%lf]\n",tx[0],tx[1],tx[2]);*/
  shift_wCM(y,ty,nn,w); /*printf("ty=[%lf,%lf,%lf]\n",ty[0],ty[1],ty[2]);*/

  nnn=0;  for(k=0;k<nn;k++)  if(w[k]>zero)  nnn++;

  for(i=0;i<3;i++) /*obtain matrix r[][]*/
    for(j=0;j<3;j++){
      rij=0;
      for(k=0;k<nn;k++)	rij+=w[k]*(y[k][i]+ty[i])*(x[k][j]+tx[j]);
      r[i][j]=rij;
    }
  
  /*find determinant of r[][]*/
  detR =r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1]);
  detR-=r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0]);
  detR+=r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]);

  for(i=0;i<3;i++) /*find matrix rr*/
    for(j=0;j<3;j++){
      rij=0;
      for(k=0;k<3;k++) rij+=r[k][i]*r[k][j];
      rr[i][j]=rij;
    }

  /*we obtain best rotation matrix U by solving the matrix equation
                      U x ( rr + L ) = r
    where L is an (unknown symmetric matrix of lagrange multipliers. If
    we impose U orthogonal, then we can multiply by traspose of rr (rrT)
                   (rr + L) x (rr + L) = r x rT
    both (r x rT) and (rr + L) are symmetric and positive definite. Thus 
    their eigenvectors are the same, and the positive eigenvalues of
    the first are the square values of the eigenvalues of the second.
    
    let mu[0-2] the three eigenvalues, and aa[0-2][-3] the three 
    eigenvectors. Then
                U[i][j]=sum_k b[k][i] * aa[k][j],
    where
                b[k][i]=sum_l r[k][l] * aa[l][i] / sqrt(mu[k])


  /* characteristic third order polinomial p(mu)= det|rr -mu^2I| */  
  p[0]=1;
  p[1]=-(rr[0][0]+rr[1][1]+rr[2][2]);
  p[2] =rr[1][1]*rr[2][2]-rr[1][2]*rr[2][1];
  p[2]+=rr[2][2]*rr[0][0]-rr[2][0]*rr[0][2];
  p[2]+=rr[1][1]*rr[0][0]-rr[1][0]*rr[0][1];
  p[3] =rr[0][0]*(rr[1][1]*rr[2][2]-rr[1][2]*rr[2][1]);
  p[3]-=rr[0][1]*(rr[1][0]*rr[2][2]-rr[1][2]*rr[2][0]);
  p[3]+=rr[0][2]*(rr[1][0]*rr[2][1]-rr[1][1]*rr[2][0]);
  p[3]=-p[3];
  /*  for(i=0;i<=3;i++)
    printf("%lf\n",p[i]);*/
  
  {
    double z=-p[1];
    for(i=2;i<4;i++)
      if(z<fabs(p[i]))z=fabs(p[i]);
    k=3;
    /* solution of the characteristic equation for eigenvalues */
    for(i=0;i<2;i++)
      {
	double q1;
	double dz;
	do 
	  {
	    q[0]=p[0];
	    q1=p[0];
	    for(j=1;j<k;j++)
	      {
		q[j]=q[j-1]*z+p[j];
		q1=q1*z+q[j];
	      }
	    q[k]=q[k-1]*z+p[k];
	    dz=q[k]/q1;
	    z=z-dz;
	  }while(dz>1.0e-14*z);
	p[k]=z;
	for(j=0;j<k;j++)
	  p[j]=q[j];
	k--;
      }
    p[1]=-p[1]/p[0];
  }
  /*
    for(i=0;i<=3;i++)
    printf("%lf\n",p[i]);
  */
  for(k=0;k<nn;k++)
    for(i=0;i<3;i++){
      cx=x[k][i]+tx[i];
      cy=y[k][i]+ty[i];
      sum+=w[k]*(cx*cx+cy*cy);
    }
  
  sum*=0.5;
  
  for(i=0;i<3;i++)
    {
      mu[i]=sqrt(p[i+1]);
      sum-=mu[i];
    }   

  /*find minimum of mu[] and change sign, to avoid improper rotation*/
  i=1;
  if(detR<0)
    {
      for(j=2;j<=3;j++)
	  if(p[j]<p[i])i=j;
      sum+=2*mu[i-1];
      mu[i-1]=-mu[i-1];  
  }

  for(k=0;k<3;k++)
    {
      int imax;
      double bmax;
      for(i=0;i<3;i++)
	{
	  for(j=0;j<3;j++)
	    b[i][j]=rr[i][j];
	  b[i][i]-=p[k+1];
	}
      bmax=0;
      imax=-1;
      for(j=0;j<3;j++)
	if(fabs(b[j][0])>bmax)
	  {
	    bmax=fabs(b[j][0]);
	    imax=j;
	  }
      if(bmax)
	{
	  if(imax)
	    for(i=0;i<3;i++)
	      {
		double c=b[0][i];
		b[0][i]=b[imax][i];
		b[imax][i]=c;  
	      }
          
	  for(i=1;i<3;i++) 
	    {
	      bmax=-b[i][0]/b[0][0];       
	      for(j=0;j<3;j++)
		b[i][j]+=bmax*b[0][j];        
	    }
	  bmax=0;
	  imax=-1;
	  for(j=1;j<3;j++)
	    if(fabs(b[j][1])>bmax)
	      {
		bmax=fabs(b[j][1]);
		imax=j;
	      }
	  if(bmax)
	    {
	      aa[2][k]=1;
	      aa[1][k]=-b[imax][2]/b[imax][1];
	    }
	  else
	    {
	      aa[2][k]=0;
	      aa[1][k]=1;
	    }
	  aa[0][k]=-(aa[1][k]*b[0][1]+aa[2][k]*b[0][2])/b[0][0];
	}
      else
	{
	  aa[0][k]=1;
	  bmax=0;
	  imax=-1;
	  for(j=0;j<3;j++)
	    if(fabs(b[j][1])>bmax)
	      {
		bmax=fabs(b[j][1]);
		imax=j;
	      }
	  if(bmax)
	    {
	      aa[2][k]=1;
	      aa[1][k]=-b[imax][2]/b[imax][1];
	    }
	  else
	    {
	      aa[2][k]=0;
	      aa[1][k]=1;
	    }
	}
    }
  
  for(k=0;k<3;k++)
    {
      double ak=0;
      for(j=0;j<3;j++)
	ak+=aa[j][k]*aa[j][k];
      ak=1/sqrt(ak);
      for(j=0;j<3;j++)
	aa[j][k]*=ak;
    }


  for(k=0;k<3;k++)/*calculate b[][k]*/
    {
      for(i=0;i<3;i++)
	{
	  double bki=0;
	  for(j=0;j<3;j++)
	    bki+=r[i][j]*aa[j][k];
	  b[i][k]=bki/mu[k];
	}
    }

  /* u is rotatation matrix : x'=ux */
  /* such that rms(x'-y) is minimal */
  for(i=0;i<3;i++)/*calculate u[][]*/
    {
      for(j=0;j<3;j++)
	{
	  double uij=0;
	  for(k=0;k<3;k++)
	    uij+=b[i][k]*aa[j][k];
	  u[i][j]=uij;	  /*printf("u[%d][%d]=%lf\n",i,j,uij);*/
	}
    }
    
  /*rms is equal to the half of the sum of the squares of coordinates 
    y and x minus sum of kabash eigenvalues mu */ 
  for(i=0; i<3; i++){
    for(j=0;j<3;j++){ 
      U[i][j] = u[i][j];
      /*printf("U[%d][%d]=%lf ",i,j,U[i][j]);*/
    } 
    /*printf("\n");*/
  }
  return sqrt(fabs(sum*2/nnn)); /*divide by number of non-zero weigths*/
}
/*=======================================================================*/

#endif /*_MY_RMSD_*/
