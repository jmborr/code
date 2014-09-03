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
/*superimpose first structure (x) to second structure (y) according to translation vectors
tx, ty and rotation matrix U  (x[i] -> U(x[i]+tx) - ty */
void supCAs(char*templx,char*temply,char*outname,double**x,double **y,int nx, int ny,
	    double*tx=NULL,double*ty=NULL,double**U=NULL,int shift=5001){
  double **x2;
  string a("A"),b("B");
  if(tx && ty && U) x2=adapt(x,y,nx,tx,ty,U);
  else x2=x;
  /*read input structure*/
  listPDBatom cachainA, cachainB;
  PDBchain inchainA(templx);
  inchainA.createCAchain(cachainA);  /*first  structure*/
  cachainA.setChainID(a);
  cachainA.renumberFullyAtomSerialNumber(1);
  cachainA.dump_coordinates_to_chain(x2,nx);
  PDBchain inchainA(temply);  
  inchainB.createCAchain(cachainB);  /*second structure*/
  cachainB.setChainID(b);
  cachainB.renumberFully(shift);
  cachainB.renumberFullyAtomSerialNumber(shift);
  cachainB.dump_coordinates_to_chain(y,ny);
  /*output the superposition*/
  ofstream outp(outname,ios::app);
  /*initialize the rasmol command lines*/
  outp<<"load inline\nselect all\nset connect false\nbackbone\ncolor red\n";
  outp<<"select :A\nbackbone 80\nselect :B\nbackbone 30\n";
  outp<<cachainA.blue2redRasmol(); /*rasmol commands to color each atom from blue to red*/
  outp<<cachainB.blue2redRasmol(); /*rasmol commands to color each atom from blue to red*/
  outp<<"select all\nexit\n";
  /*output the pdb files*/
  outp<<cachainA<<"TER\n";
  outp<<cachainA.connectRasmol();
  outp<<cachainB<<"TER\n";
  outp<<cachainB.connectRasmol();
  outp.close();
}
/*=======================================================================*/
