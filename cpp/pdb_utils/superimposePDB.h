#include<iostream>
#include "pdbClasses2.h"
#include "allocate.h"

/*=======================================================================*/
/*adapt x into y, but do not modify x*/
double **adapt(double**x,double **y,int n,double*tx,double*ty,double**U){
  double **x2=alloc2<double>(n,3,x);
  double a,b,c;
  for(int i=0;i<n;i++){
    a=x2[i][0]+tx[0]; b=x2[i][1]+tx[1]; c=x2[i][2]+tx[2];
    x2[i][0] = U[0][0]*a + U[0][1]*b + U[0][2]*c - ty[0];
    x2[i][1] = U[1][0]*a + U[1][1]*b + U[1][2]*c - ty[1];
    x2[i][2] = U[2][0]*a + U[2][1]*b + U[2][2]*c - ty[2];
  }
  return x2;
}
/*=======================================================================*/
/*superimpose first structure to second structure*/
void supCAs(char*templ,char*outname,double**x,double **y,int n,
	    double*tx=NULL,double*ty=NULL,double**U=NULL,int shift=5001){
  double **x2;
  if(tx && ty && U) x2=adapt(x,y,n,tx,ty,U);
  else x2=x;
  /*read input structure*/
  PDBchain inchain(templ);
  listPDBatom cachain;
  inchain.createCAchain(cachain);
  /*output the superposition*/
  ofstream outp(outname,ios::app);
  /*first structure*/
  cachain.renumberFullyAtomSerialNumber(1);
  cachain.dump_coordinates_to_chain(y,n);
  outp<<cachain<<"TER\n";
  outp<<cachain.connectRasmol();
  /*second structure*/
  cachain.renumberFully(shift);
  cachain.renumberFullyAtomSerialNumber(shift);
  cachain.dump_coordinates_to_chain(x2,n);
  outp<<cachain<<"TER\n";
  outp<<cachain.connectRasmol();
  outp.close();
}
