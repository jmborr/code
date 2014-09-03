#ifndef _KRMS_H_
#define _KRMS_H_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;
/*---------------------------------------------------------- 
  
To calculate the rmsd between two sets of vectors x[nn][3] and y[nn][3]; the
proper rotation matrix is stored in rotate[3][3]; the residual rmsd for each
residue is stored in dr[nn]; The input sets of vectors should have redueced
the translation, so that the COM of the two sets should sit at the orig.

REF: Kabsch, W., Acta Cryst, 1978, 827-828

implemented by Sergey V. Buldyrev
adding proper rotation matrix correction by Feng Ding
---------------------------------------------------------*/
double get_rms(double **x, double ** y, int nn,
	       double * dr=NULL, double** new_x=NULL, double** rotate=NULL){
  double r[3][3], rr[3][3], mu[3], p[4], q[4];
  double aa[3][3], b[3][3], u[3][3];
  int i,j,n,m,k;
  double sum=0;
  /*Construct the |r|===YX;
    matrix S===XX; where L is a multiplier and L=Trans(L); 
    the least square value exist at u(S+L)=r, Trans(u)*u=1 is the rotation
    matrix*/
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      double  rij=0;
      for(k=0;k<nn;k++){
	rij+=y[k][i]*x[k][j];
      }
      r[i][j]=rij;
    }
  /*R is an asysmmetric matrix, consruct rr===Trans(r)*r=(S+L);*/
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      double  rij=0;
      for(k=0;k<3;k++){
	rij+=r[k][i]*r[k][j];
      }
      rr[i][j]=rij;
    }
  }
  /*rr is a symmetric positive definit matrix, solve the eigen system*/
  
  /*charcteristic third order polinomial p(mu)= det|rr -mu^2I| 
    solving the equation: p[0]*x^3+p[1]*x^2+p[2]*x+p[3]=0*/
  p[0]=1;
  p[1]=-(rr[0][0]+rr[1][1]+rr[2][2]);
  p[2] =rr[1][1]*rr[2][2]-rr[1][2]*rr[2][1];
  p[2]+=rr[2][2]*rr[0][0]-rr[2][0]*rr[0][2];
  p[2]+=rr[1][1]*rr[0][0]-rr[1][0]*rr[0][1];
  p[3] =rr[0][0]*(rr[1][1]*rr[2][2]-rr[1][2]*rr[2][1]);
  p[3]-=rr[0][1]*(rr[1][0]*rr[2][2]-rr[1][2]*rr[2][0]);
  p[3]+=rr[0][2]*(rr[1][0]*rr[2][1]-rr[1][1]*rr[2][0]);
  p[3]=-p[3];
  /*for(i=0;i<=3;i++)
    printf("%lf\n",p[i]);*/
  
  {
    double z=-p[1];
    for(i=2;i<4;i++)
      if(z<fabs(p[i]))z=fabs(p[i]);
    k=3;
    /* solution of the characteristic equation for eigenvalues */
    /*f(x)=P0x^3+P1x^2+P2*x+P3
      f(x')=f(x)+f'(x)*(x'-x)=0
      x'-x=dx=-f(x)/f'(x)
      therefore, it can be solved by a loop.
      once a solution is obtained.
      f(x)/(x-z) = [P0*x^3+P1*x^2+P2*x-(P0*z^3+P1*z^2+P2*z)]/(x-z)
      =P0*(x^2+x*z+z^2) + P1*(x+Z) + P2
      =P0*x^2 + (P0*z+P1)*x + (P0*z^2+P1*z+P2)
      another iteration!!!
    */
    for(i=0;i<2;i++){
      double q1;
      double dz;
      do {
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
  /*
    for(i=0;i<=3;i++)
    printf("%lf\n",p[i]);
  */
  for(k=0;k<nn;k++)
    for(i=0;i<3;i++)
      sum+=x[k][i]*x[k][i]+ y[k][i]*y[k][i]; 
  
  sum*=0.5;
  
  /*rms is equal to the hulf of the sum of the squares of coordinates 
    y and x minus sum of kabash eigenvalues mu 
    ----This statement might not be always true. Because the eigenvalue of 
    r equeals the eigenvalue of rr times an arbitrory sign. 
    --- Sort the eigen values so that p[0]>=p[1]>=p[2]>=0.
    If det(r)>0; rmsd^2=E0/2-sqrt[P0]-sqrt[P1]-sqrt[P2];
    otherwise, rmsd^2 = E0/2-sqrt[P0]-sqrt[P1]+sqrt[P2];
    In summary, rmsd^2 = E0/2 - sqrt[P0]-sqrt[P1]-sgn(determ(r))sqrt[P2];
    infact, only the min(p0,p1,p2) is needed!
  */
  double min=p[3];
  int min_index=3;
  for(int i=2; i>0; i--)
    if(p[i]<min){
      min = p[i];
      min_index=i;
    }
  /*printf("%ld\n", min_index)*/
  double dummy = r[0][0]*r[1][1]*r[2][2]+r[1][0]*r[2][1]*r[0][2]+r[2][0]*r[1][2]*r[0][1];
  dummy -= r[0][0]*r[2][1]*r[1][2]+r[0][1]*r[1][0]*r[2][2]+r[0][2]*r[1][1]*r[2][0];
  int determ_sgn = (dummy>0) ? 1 : -1;
  for(i=0;i<3;i++){
    if((i+1)==min_index && determ_sgn<0)
      mu[i]=-sqrt(p[i+1]);
    else
      mu[i]=sqrt(p[i+1]);
    sum-=mu[i];
  }
  /*find the eigen vectors a*/
  if(dr || new_x || rotate){
    for(k=0;k<3;k++){
      int imax;
      double bmax;
      for(i=0;i<3;i++){
	for(j=0;j<3;j++)
	  b[i][j]=rr[i][j];
	b[i][i]-=p[k+1];
      }
      bmax=0;
      imax=-1;
      for(j=0;j<3;j++)
	if(fabs(b[j][0])>bmax){
	  bmax=fabs(b[j][0]);
	  imax=j;
	}
      if(bmax){
	if(imax)
	  for(i=0;i<3;i++){
	    double c=b[0][i];
	    b[0][i]=b[imax][i];
	    b[imax][i]=c;  
	  }
	
	for(i=1;i<3;i++){
	  bmax=-b[i][0]/b[0][0];       
	  for(j=0;j<3;j++)
	    b[i][j]+=bmax*b[0][j];        
	}
	bmax=0;
	imax=-1;
	for(j=1;j<3;j++)
	  if(fabs(b[j][1])>bmax){
	    bmax=fabs(b[j][1]);
	    imax=j;
	  }
	if(bmax){
	  aa[2][k]=1;
	  aa[1][k]=-b[imax][2]/b[imax][1];
	}
	else{
	  aa[2][k]=0;
	  aa[1][k]=1;
	}
	aa[0][k]=-(aa[1][k]*b[0][1]+aa[2][k]*b[0][2])/b[0][0];
      }
      else{
	aa[0][k]=1;
	bmax=0;
	imax=-1;
	for(j=0;j<3;j++)
	  if(fabs(b[j][1])>bmax){
	    bmax=fabs(b[j][1]);
	    imax=j;
	  }
	if(bmax){
	  aa[2][k]=1;
	  aa[1][k]=-b[imax][2]/b[imax][1];
	}
	else{
	  aa[2][k]=0;
	  aa[1][k]=1;
	}
      }
    }
    /*Normalize*/
    for(k=0;k<3;k++){
      double ak=0;
      for(j=0;j<3;j++)
	ak+=aa[j][k]*aa[j][k];
      ak=1/sqrt(ak);
      for(j=0;j<3;j++)
	aa[j][k]*=ak;
    }
    /*b=|r|*a*/
    for(k=0;k<3;k++){
      for(i=0;i<3;i++){
	double bki=0;
	for(j=0;j<3;j++)
	  bki+=r[i][j]*aa[j][k];
	b[i][k]=bki/mu[k];
      }
    }
    
    
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
	double uij=0;
	for(k=0;k<3;k++)
	  uij+=b[i][k]*aa[j][k];
	u[i][j]=uij;
      }
    }
    /* is rotatation matrix : x'=ux */
    /* such that rms(x'-y) is minimal */
    dummy=0;
    double tdummy=0;
    for(k=0;k<nn;k++){
      if(dr)dr[k]=0;
      for(j=0;j<3;j++){
	double ykj=0;
	for(i=0;i<3;i++)
	  ykj+=u[j][i]*x[k][i];
	if(new_x)new_x[k][j]=ykj;
	if(dr)dr[k]+=(ykj-y[k][j])*(ykj-y[k][j]);
	dummy+=(ykj-y[k][j])*(ykj-y[k][j]);
      }
    }
    /*printf("%lf %lf\n", dummy, sum);*/
    if(fabs(dummy-sum)<1e-4){
      printf("error happens!\n");
      exit(2);
    }
    if(rotate){
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++) rotate[i][j]=u[i][j];
    }
    /* dr[k] accumulates rms for k-yh atom from many calls of
       this function  the rotated coordinates x are ykj , k
       goes from 0 to nn-1, j is 0(x), 1(y), or 2(z). */
  }
  return sum;
}
#endif/*_KRMS_H_*/
