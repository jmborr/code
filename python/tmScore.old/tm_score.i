/*
swig -c++ -python -o tm_score_cpp_wrap.cpp tm_score.i && g++ -fPIC -Wno-deprecated -c pdbClasses2.cpp && gcc -fPIC -c tm_score_cpp_wrap.cpp -o tm_score_cpp_wrap.o -I/usr/include/python2.3/ && g++ -shared tm_score_cpp_wrap.o pdbClasses2.o -o _tm_score_cpp.so
*/

%module tm_score_cpp
%{
#include "allocate.h"
#include "rmsd.h"
#include "tm_score.h"
#include "createCAchain.h"
#include "superimposePDB.h"
%}

%include cpointer.i
%include "std_string.i"

// /gpfs1/active/jose/code/cpp/rmsd.h
extern double getrmsd(double **x,double **y,int nn,double *w=NULL);

// /gpfs1/active/jose/code/cpp/pdb_utils/source/createCAchain.h
extern int getCAcoords2(double **x,char *inpdb);

// /gpfs1/active/jose/code/cpp/tm_score/tm_score.h
extern double getd0(int nn);
extern double coverage(double**x,double **y,int n,double*tx,double*ty,double**U,int *align=NULL,double d1=0);
extern double coverageII(double**x,double **y,int n,double*tx,double*ty,double**U,double d1=0.);
extern double rmsd_cov(double**x,double **y,int n,double*tx,double*ty,double**U,double d1=0.);
extern double gettm(double **x, double **y, int nn, double d0=0., double *tx=NULL, double *ty=NULL, double **U=NULL, int Lmin=18, int n_chunks=10, double tmco=0.85);

// /gpfs1/active/jose/code/cpp/pdb_utils/superimposePDB.h
extern void supCAs(char*templx,char*temply,char*outname,double**x,double **y,int nx,int ny,double*tx=NULL,double*ty=NULL,double**U=NULL,int shift=5001);

/*create a double array from a list*/
%typemap(python,in) (double *) {
  if (PyList_Check($input)) {
    int i,size = PyList_Size($input);
    $1 = new double[size];
    for(i=0; i<size; i++){
      PyObject *o = PyList_GetItem($input,i);
      if (PyFloat_Check(o))
	$1[i] = PyFloat_AsDouble(o);
      else {
	PyErr_SetString(PyExc_TypeError,"list must contain doubles");
	delete $1;
	return NULL;
      }
    }
  }
  else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

%inline %{

  void print_array(int n, int m, double **x){
    printf("init print_array\n");
    int i,j,nm=n*m;
    double **z=new double *[n];
    z[0]=x[0];
    for(i=1;i<n;i++) z[i]=z[i-1]+m;
    for(i=0;i<n;i++){
      for(j=0;j<m;j++) printf(" %lf",z[i][j]);
      printf("\n");
    }
//  printf("finish print_array\n");
  }

// /gpfs1/active/jose/code/cpp/allocate/allocate.h

  double **allocD2(int n1,int n2) { return alloc2<double>(n1,n2); }

  void copyD2(int n1, int n2, double **r, double *x){
    int i,nm;
    nm=n1*n2;
    for(i=0;i<nm;i++){ 
      r[0][i]=x[i]; 
    }
  }


  double *allocD1(int n) { return new double[n]; }

  int *allocI1(int n) { return new int[n]; } 

  std::string lowerNotAligned(std::string seq2,int *alg){
    int i,n=seq2.length();
    std::string seq=seq2;
    for(i=0;i<n;i++)
      if(alg[i]==0){
        seq[i]=tolower(seq[i]);
      }
    return seq;
  }
%}
