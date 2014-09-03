#include<iostream>
#include <math.h>
#include "rmsd.h"

#define LONGSEQ 10000
#define INF_PASSES 5000
/* 0.0 < TM-score < 0.17, Random predictions                           
   0.4 < TM-score < 1.00, Meaningful predictions  */

int nrms=0; /*number of rmsd's for a TM-score computation*/

/*=======================================================================*/
double d2(double*x,double*y){
  return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]);
}
/*=======================================================================*/
/*d0 eliminates size dependency of rmsd and tm-score*/
double getd0(int nn){
  if(nn<18) return -1;
  else return 1.24*pow(nn-15,0.333333)-1.8;
}
/*=======================================================================*/
/*return percentage of residues with distance below d0. It's assume we did tm before
  optional parameter align[i]==1 iif residues below d0
  see below coverageII when we want to pass d1 but not align
*/
double coverage(double**x,double **y,int n,double*tx,double*ty,double**U,
		int *align=NULL,double d1=0){
  double **x2=alloc2<double>(n,3,x);
  double d_2;
  if(d1>0) d_2=d1*d1;
  else{ d_2=getd0(n);  d_2*=d_2; }
  double a,b,c;
  int i;
  for(i=0;i<n;i++){
    a=x2[i][0]+tx[0]; b=x2[i][1]+tx[1]; c=x2[i][2]+tx[2];
    x2[i][0] = U[0][0]*a + U[0][1]*b + U[0][2]*c - ty[0];
    x2[i][1] = U[1][0]*a + U[1][1]*b + U[1][2]*c - ty[1];
    x2[i][2] = U[2][0]*a + U[2][1]*b + U[2][2]*c - ty[2];
  }
  double cov=0;
  if(!align){
    for(i=0;i<n;i++) if( d2(x2[i],y[i]) < d_2 ) cov++;
  }
  else{
    for(i=0;i<n;i++){
      if( d2(x2[i],y[i]) < d_2 ){
	cov++;
	align[i]=1;
      }
      else align[i]=0;
      /*printf("%d %d, ",i,align[i]);*/
    }
  }
  return cov/n;
}
/*=======================================================================*/
double coverageII(double**x,double **y,int n,double*tx,double*ty,double**U,double d1=0){
  int *pt=new int[n]; /*pass a fake align pointer*/
  double cov=coverage(x,y,n,tx,ty,U,pt,d1);
  delete pt;
  return cov;
}
/*=======================================================================*/
double rmsd_cov(double**x,double **y,int n,double*tx,double*ty,double**U,double d1=0){
  double **x2=alloc2<double>(n,3,x);
  double d_2;
  if(d1>0) d_2=d1*d1;
  else{ d_2=getd0(n);  d_2*=d_2; }
  /*std::cout<<"d1="<<d1<<" d_2="<<d_2<<std::endl;*/
  double a,b,c;
  int i;
  for(i=0;i<n;i++){
    a=x2[i][0]+tx[0]; b=x2[i][1]+tx[1]; c=x2[i][2]+tx[2];
    x2[i][0] = U[0][0]*a + U[0][1]*b + U[0][2]*c - ty[0];
    x2[i][1] = U[1][0]*a + U[1][1]*b + U[1][2]*c - ty[1];
    x2[i][2] = U[2][0]*a + U[2][1]*b + U[2][2]*c - ty[2];
  }
  double *w=new double[n];
  for(i=0;i<n;i++) if( d2(x2[i],y[i]) < d_2 ){
    /*std::cout<<sqrt(d2(x2[i],y[i]))<<std::endl; */
    w[i]=1.0;
  }
  else w[i]=0.0;
  return getrmsd(x,y,n,w);
}
/*=======================================================================*/
/*find rmsd, then update the weights*/
double do_a_pass(double **x, double *tx, double **y, double *ty, long nn,
		  double *di2, double d0, double *w, double **U){
  int i;
  double tm=0,a,b,c,A,B,C;
  double d02=d0*d0;
  double rmsd;
  /*store in tx translation such that center of mass of tx would be in origin. The "mass"
    associated to each position "i" is the weigh w[i]*/
  rmsd=get_rot(x,tx,y,ty,nn,w,U); nrms++; /*printf("rmsd=%lf\n",rmsd);*/
  /*for(i=0;i<3;i++){for(int j=0;j<3;j++)printf("U[%d][%d]=%lf ",i,j,U[i][j]);printf("\n");}*/
  for(i=0;i<nn;i++){
    a=x[i][0]+tx[0]; b=x[i][1]+tx[1]; c=x[i][2]+tx[2];
    A = U[0][0]*a + U[0][1]*b + U[0][2]*c - y[i][0] - ty[0];
    B = U[1][0]*a + U[1][1]*b + U[1][2]*c - y[i][1] - ty[1];
    C = U[2][0]*a + U[2][1]*b + U[2][2]*c - y[i][2] - ty[2];
    di2[i]=A*A+B*B+C*C;   /*printf("di2[%3d]=%lf\n ",i,di2[i]);*/
    a=1/(1+di2[i]/d02);
    w[i]=a*a;  /*update weights*/
    tm+=a;     /*calculate tm*/
  }
  return tm;
}
/*=======================================================================*/
double gettm(double **x, double **y, int nn, double d0=0., double *tx=NULL, double *ty=NULL,
	     double **U=NULL, int Lmin=18, int n_chunks=10, double tmco=0.85){
  int i,j,L,min,max,p,pmax,step;
  double err=0.01*nn; /*error tolerance*/
  double tm,tmM=0,tmprev,nn_tmco=tmco*nn;

  /*allocate once, use in subsequent calls! WH:work-horse*/
  static double *txWH=new double[3];
  static double *tyWH=new double[3];
  static double *wWH=new double[LONGSEQ];
  static double *di2WH=new double[LONGSEQ];
  static double **UWH=new double*[3];
  UWH[0]=new double[9]; /*this is NOT static, thus need to deallocate later*/
  for(i=1;i<3;i++) UWH[i]=UWH[i-1]+3; 

  /*generate defaults arguments, if not passed*/
  bool uf=false, txf=false, tyf=false;
  if(d0<0.01) d0=getd0(nn);
  if(!U){ uf=true; U=new double*[3]; U[0]=new double[9]; for(i=1;i<3;i++) U[i]=U[i-1]+3; }
  if(!tx){ txf=true; tx=new double[3]; }
  if(!ty){ tyf=true; ty=new double[3]; }

  nrms=0;
  if(Lmin<4) Lmin=4; /*minimum chunk length*/
  if(nn<3*Lmin) Lmin=int(nn/3);
  step=nn/n_chunks; if(step<1) step=1; /*jumps in the sequence*/
  /*printf("Lmin=%d,n_chunks=%d,tmco=%lf,step=%d\n",Lmin,n_chunks,tmco,step);*/
  L=nn;
  while(L>=Lmin){
    for(min=0;min<=nn-L;min+=step){
      max=min+L;
      for(i=0;i<nn;i++) wWH[i]=0.;    /*initialize weigths*/
      for(i=min;i<max;i++) wWH[i]=1.;
      tmprev=do_a_pass(x,txWH,y,tyWH,nn,di2WH,d0,wWH,UWH);
      if(tmprev>tmM){
	tmM=tmprev;
	for(i=0;i<3;i++){tx[i]=txWH[i];ty[i]=tyWH[i];}
	for(i=0;i<9;i++)U[0][i]=UWH[0][i];
      }
      tm=do_a_pass(x,txWH,y,tyWH,nn,di2WH,d0,wWH,UWH);
      if(tm>tmM){ /*we got a better tm-score*/
	tmM=tm;
	for(i=0;i<3;i++){tx[i]=txWH[i];ty[i]=tyWH[i];}
	for(i=0;i<9;i++)U[0][i]=UWH[0][i];
      }
      while(fabs(tmprev-tm)>err){
	tmprev=tm;
	tm=do_a_pass(x,txWH,y,tyWH,nn,di2WH,d0,wWH,UWH);
	if(tm>tmM){
	  tmM=tm;
	  for(i=0;i<3;i++){tx[i]=txWH[i];ty[i]=tyWH[i];}
	  for(i=0;i<9;i++)U[0][i]=UWH[0][i];
	}
      }      
      if(tmM>nn_tmco) goto similar; /*structures are already very similar*/
    }
    L/=2;
  }
 similar:
  /*deallocate memory*/
  delete [] UWH[0];
  if(uf) delete [] U;  if(txf) delete [] tx;  if(tyf) delete [] ty;
  /*printf("nrms=%d\n",nrms);*/
  return tmM/nn;
}
/*===========================================================*/
double gettmBAK(double **x, double **y, int nn, double d0=0., double *tx=NULL, double *ty=NULL,
	     double **U=NULL, int Lmin=18, int n_chunks=10, double tmco=0.85){
  int i,j,L,min,max,p,pmax,step;
  double err=0.01*nn; /*error tolerance*/
  double tm,tmM=0,tmprev,nn_tmco=tmco*nn;
  double *w=new double[nn], *di2=new double[nn];

  /*for(int i=0;i<nn;i++)  printf("%lf %lf %lf\n",x[i][0],x[i][1],x[i][2]);
    for(int i=0;i<nn;i++)  printf("%lf %lf %lf\n",y[i][0],y[i][1],y[i][2]); */
  
  /*allocate to store best transformation*/
  double *txb=new double[3], *tyb=new double[3];
  double **Ub=new double*[3]; Ub[0]=new double[9]; for(i=1;i<3;i++) Ub[i]=Ub[i-1]+3; 

  /*generate defaults arguments, if not passed*/
  if(d0<0.01) d0=getd0(nn);
  if(!U){ U=new double*[3]; U[0]=new double[9]; for(i=1;i<3;i++) U[i]=U[i-1]+3; }
  if(!tx) tx=new double[3]; 
  if(!ty) ty=new double[3];

  nrms=0;
  if(Lmin<4) Lmin=4; /*minimum chunk length*/
  step=nn/n_chunks; if(step<1) step=1; /*jumps in the sequence*/
  /*printf("Lmin=%d,n_chunks=%d,tmco=%lf,step=%d\n",Lmin,n_chunks,tmco,step);*/
  L=nn;
  while(L>=Lmin){
    for(min=0;min<=nn-L;min+=step){
      max=min+L;
      for(i=0;i<nn;i++) w[i]=0.;    /*initialize weigths*/
      for(i=min;i<max;i++) w[i]=1.;
      tmprev=do_a_pass(x,tx,y,ty,nn,di2,d0,w,U);
      if(tmprev>tmM){
	tmM=tmprev;
	for(i=0;i<3;i++){txb[i]=tx[i];tyb[i]=ty[i];}
	for(i=0;i<9;i++)Ub[0][i]=U[0][i];
      }
      tm=do_a_pass(x,tx,y,ty,nn,di2,d0,w,U);
      if(tm>tmM){ /*we got a better tm-score*/
	tmM=tm;
	for(i=0;i<3;i++){txb[i]=tx[i];tyb[i]=ty[i];}
	for(i=0;i<9;i++)Ub[0][i]=U[0][i];
      }
      while(fabs(tmprev-tm)>err){
	tmprev=tm;
	tm=do_a_pass(x,tx,y,ty,nn,di2,d0,w,U);
	if(tm>tmM){
	  tmM=tm;
	  for(i=0;i<3;i++){txb[i]=tx[i];tyb[i]=ty[i];}
	  for(i=0;i<9;i++)Ub[0][i]=U[0][i];
	}
      }      
      if(tmM>nn_tmco) goto similar; /*structures are already very similar*/
    }
    L/=2;
  }
 similar:
  for(i=0;i<3;i++){tx[i]=txb[i];ty[i]=tyb[i];}
  for(i=0;i<9;i++)U[0][i]=Ub[0][i];
  /*printf("nrms=%d\n",nrms);*/
  return tmM/nn;
}
/*===========================================================*/
/*TM-alignment with the help of secondary structure chunks*/
double gettmSeqDat(double **x, double **y, int nn, int* beginSec, int* resperSec, int nSec, double d0=0.0, double *tx=NULL, double *ty=NULL, double **U=NULL, int maxnrms=INF_PASSES,int step=1,int Lmin=1,double tmco=0.90){
  int i,j,L,min,max,p,pmax;
  int begin,end;
  double err=0.01*nn; /*error tolerance*/
  double tm,tmM=0,tmprev,nn_tmco=tmco*nn;
  int nrmsperpass,maxnrmsperpass=4;

  /*allocate once, use in subsequent calls! WH:work-horse*/
  static double *txWH=new double[3];
  static double *tyWH=new double[3];
  static double *wWH=new double[LONGSEQ];
  static double *di2WH=new double[LONGSEQ];
  static double **UWH=new double*[3];
  UWH[0]=new double[9]; /*this is NOT static, thus need to deallocate later*/
  for(i=1;i<3;i++) UWH[i]=UWH[i-1]+3; 

  /*generate defaults arguments, if not passed*/
  bool uf=false, txf=false, tyf=false;
  if(d0<0.01) d0=getd0(nn);
  if(!U){ uf=true; U=new double*[3]; U[0]=new double[9]; for(i=1;i<3;i++) U[i]=U[i-1]+3; }
  if(!tx){ txf=true; tx=new double[3]; }
  if(!ty){ tyf=true; ty=new double[3]; }

  nrms=0;
  /*step=1; jump "step" secondary structures*/
  /*Lmin=1; minimum number of secondary structure segments to consider*/
  L=nSec; /*initially consider all secondary structure segments*/
  
  while(L>=Lmin){
    for(min=0;min<=nSec-L;min+=step){
      nrmsperpass=0;
      max=min+L;
      for(i=0;i<nn;i++) wWH[i]=0.;   /*initialize weigths*/
      for(i=min;i<max;i++){        /*go through every secondary structure segment*/
	begin=beginSec[i]-1;       /*where does the segment begins*/
	end=beginSec[i]+resperSec[i]-1; /*where does the segment ends?*/
	for(j=begin;j<end;j++)  wWH[j]=1.;
      }
      tmprev=do_a_pass(x,txWH,y,tyWH,nn,di2WH,d0,wWH,UWH);
      if(tmprev>tmM){
	tmM=tmprev;
	for(i=0;i<3;i++){tx[i]=txWH[i];ty[i]=tyWH[i];}
	for(i=0;i<9;i++)U[0][i]=UWH[0][i];
      }
      tm=do_a_pass(x,txWH,y,tyWH,nn,di2WH,d0,wWH,UWH);
      if(tm>tmM){ /*we got a better tm-score*/
	tmM=tm;
	for(i=0;i<3;i++){tx[i]=txWH[i];ty[i]=tyWH[i];}
	for(i=0;i<9;i++)U[0][i]=UWH[0][i];
      }
      while(fabs(tmprev-tm)>err && nrmsperpass<maxnrmsperpass){
	tmprev=tm;
	tm=do_a_pass(x,txWH,y,tyWH,nn,di2WH,d0,wWH,UWH);
	if(tm>tmM){
	  tmM=tm;
	  for(i=0;i<3;i++){tx[i]=txWH[i];ty[i]=tyWH[i];}
	  for(i=0;i<9;i++)U[0][i]=UWH[0][i];
	}
	nrmsperpass++;
      }
      if(tmM>nn_tmco) goto similar; /*structures are already very similar*/
      if(nrms>maxnrms) break; /*otherwise too many calculations*/
    }
    L/=2;
  }
 similar:
  /*deallocate memory*/
  delete [] UWH[0];
  if(uf) delete [] U;  if(txf) delete [] tx;  if(tyf) delete [] ty;
  /*printf("nsec=%2d Lmin=%d step=%d L=%2d nrmsd=%d\n",nSec,Lmin,step,L,nrms);*/
  return tmM/nn;
}
/*===========================================================*/
/*TM-alignment with the help of secondary structure chunks*/
double gettmSeqDatBAK(double **x, double **y, int nn, int* beginSec, int* resperSec, int nSec, double d0=0.0, double *tx=NULL, double *ty=NULL, double **U=NULL, double tmco=0.90){
  int i,j,L,min,max,p,pmax;
  int Lmin,step;
  int begin,end;
  double err=0.01*nn; /*error tolerance*/
  double tm,tmM=0,tmprev,nn_tmco=tmco*nn;
  double *w=new double[nn], *di2=new double[nn];

  /*allocate to store best transformation*/
  double *txb=new double[3], *tyb=new double[3];
  double **Ub=new double*[3]; Ub[0]=new double[9]; for(i=1;i<3;i++) Ub[i]=Ub[i-1]+3; 


  /*generate defaults arguments, if not passed*/
  bool uf=false, txf=false, tyf=false;
  if(d0<0.01) d0=getd0(nn);
  if(!U){ uf=true; U=new double*[3]; U[0]=new double[9]; for(i=1;i<3;i++) U[i]=U[i-1]+3; }
  if(!tx){ txf=true; tx=new double[3]; }
  if(!ty){ tyf=true; ty=new double[3]; }

  nrms=0;
  step=1; /*jump "step" secondary structures*/
  Lmin=1; /*minimum number of secondary structure segments to consider*/
  L=nSec; /*initially consider all secondary structure segments*/
  while(L>=Lmin){
    for(min=0;min<=nSec-L;min+=step){
      max=min+L;
      for(i=0;i<nn;i++) w[i]=0.;   /*initialize weigths*/
      for(i=min;i<max;i++){        /*go through every secondary structure segment*/
	begin=beginSec[i]-1;       /*where does the segment begins*/
	end=beginSec[i]+resperSec[i]-1; /*where does the segment ends?*/
	for(j=begin;j<end;j++)  w[j]=1.;
      }
      tmprev=do_a_pass(x,tx,y,ty,nn,di2,d0,w,U);
      if(tmprev>tmM){
	tmM=tmprev;
	for(i=0;i<3;i++){txb[i]=tx[i];tyb[i]=ty[i];}
	for(i=0;i<9;i++)Ub[0][i]=U[0][i];
      }
      tm=do_a_pass(x,tx,y,ty,nn,di2,d0,w,U);
      if(tm>tmM){ /*we got a better tm-score*/
	tmM=tm;
	for(i=0;i<3;i++){txb[i]=tx[i];tyb[i]=ty[i];}
	for(i=0;i<9;i++)Ub[0][i]=U[0][i];
      }
      while(fabs(tmprev-tm)>err){
	tmprev=tm;
	tm=do_a_pass(x,tx,y,ty,nn,di2,d0,w,U);
	if(tm>tmM){
	  tmM=tm;
	  for(i=0;i<3;i++){txb[i]=tx[i];tyb[i]=ty[i];}
	  for(i=0;i<9;i++)Ub[0][i]=U[0][i];
	}
      }      
      if(tmM>nn_tmco) goto similar; /*structures are already very similar*/
    }
    L/=2;
  }
 similar:
  /*load best rotation and translations*/
  if(txf && tyf)  for(i=0;i<3;i++){tx[i]=txb[i];ty[i]=tyb[i];}
  if(uf) for(i=0;i<9;i++)U[0][i]=Ub[0][i];
  /*deallocate memory*/
  delete [] w; delete [] di2;  delete [] txb; delete [] tyb; delete [] Ub[0]; delete [] Ub;
  if(uf) delete [] U;  if(txf) delete [] tx;  if(tyf) delete [] ty;
  /*printf("nrmsd=%d\n",nrms);*/
  return tmM/nn;
}
/*===========================================================*/
