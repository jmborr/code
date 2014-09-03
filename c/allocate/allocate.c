#include<stdlib.h>
#include<stdio.h>

/*======================================================================*/
/*               INT TYPE MULTIDIMENSIONAL ARRAYS                       */
/*======================================================================*/
int *allocI1(int a){
  return (int *)malloc(a*sizeof(int));
}
/*======================================================================*/
void deallocI2(int **a){
  free(a[0]);
  free(a);
}
/*======================================================================*/
void dumpI1(int *i1, int a, int I){
  int i;
  char format[124],c='%';
  sprintf(format,"%c%dd",c,I);
  for(i=0;i<a;i++) printf(format,i1[i]);
  printf("\n");
}
/*======================================================================*/
int **allocI2(int a, int b){
  int i;
  int **tmp=(int**)malloc(a*sizeof(int*));
  tmp[0]=(int*)malloc(a*b*sizeof(int));
  for(i=1;i<a;i++) tmp[i]=tmp[i-1]+b;
  return tmp;
}
/*======================================================================*/
void dumpI2(int **d2, int a, int b, int I){
  int i,j;
  char format[124],c='%';
  sprintf(format,"%c%dd",c,I);
  for(i=0;i<a;i++){
    for(j=0;j<b;j++)
      printf(format,d2[i][j]);
    printf("\n");
  }
}
/*======================================================================*/
void dumpsymI2(int **d2, int a, int I){
  int i,j;
  char format[124],c='%';
  sprintf(format,"%c%dd",c,I);
  for(i=0;i<a;i++){
    for(j=0;j<=i;j++)
      printf(format,d2[i][j]);
    printf("\n");
  }
}
/*======================================================================*/
int ***allocI3(int a, int b, int c){
  int i,ab=a*b;
  int ***tmp=(int***)malloc(a*sizeof(int**));
  tmp[0]=(int**)malloc(ab*sizeof(int*));
  for(i=1;i<a;i++) tmp[i]=tmp[i-1]+b;
  tmp[0][0]=(int*)malloc(ab*c*sizeof(int));
  for(i=1;i<ab;i++) tmp[0][i]=tmp[0][i-1]+c;
  return tmp;
}
/*======================================================================*/
void deallocI3(int ***a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}
/*======================================================================*/
/*first index bigger or equal than second index*/
int **allocsymI2(int a){ /*only (i,j) with i<j is filled*/
  int i;
  int n=a*(a+1)/2;
  int **tmp=(int**)malloc((a-1)*sizeof(int*));
  tmp[0]=(int*)malloc(n*sizeof(int));
  for(i=1;i<a;i++) tmp[i]=tmp[i-1]+i;
  return tmp;
}
/*======================================================================*/
void deallocsymI2(int **a){
  free(a[0]);
  free(a);
}
/*======================================================================*/
/*               DOUBLE TYPE MULTIDIMENSIONAL ARRAYS                    */
/*======================================================================*/
double *allocD1(int a){
  return (double *)malloc(a*sizeof(double));
}
/*======================================================================*/
void dumpD1(double *d1, int a, int F, int D){
  int i;
  char format[124],c='%';
  sprintf(format,"%c%d.%dlf",c,F,D);
  for(i=0;i<a;i++) printf(format,d1[i]);
  printf("\n");
}
/*======================================================================*/
double **allocD2(int a, int b){
  int i;
  double **tmp=(double**)malloc(a*sizeof(double*));
  tmp[0]=(double*)malloc(a*b*sizeof(double));
  for(i=1;i<a;i++) tmp[i]=tmp[i-1]+b;
  return tmp;
}
/*======================================================================*/
void deallocD2(double **a){
  free(a[0]);
  free(a);
}
/*======================================================================*/
void dumpD2(double **d2, int a, int b, int F, int D){
  int i,j;
  char format[124],c='%';
  sprintf(format,"%c%d.%dlf",c,F,D);
  for(i=0;i<a;i++){
    for(j=0;j<b;j++)
      printf(format,d2[i][j]);
    printf("\n");
  }
}
/*======================================================================*/
double ***allocD3(int a, int b, int c){
  int i,ab=a*b;
  double ***tmp=(double***)malloc(a*sizeof(double**));
  tmp[0]=(double**)malloc(ab*sizeof(double*));
  for(i=1;i<a;i++) tmp[i]=tmp[i-1]+b;
  tmp[0][0]=(double*)malloc(ab*c*sizeof(double));
  for(i=1;i<ab;i++) tmp[0][i]=tmp[0][i-1]+c;
  return tmp;
}
/*======================================================================*/
void deallocD3(double ***a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}
/*======================================================================*/
/*               CHAR TYPE MULTIDIMENSIONAL ARRAYS                    */
/*======================================================================*/
char *allocC1(int a){
  return (char *)malloc(a*sizeof(char));
}
/*======================================================================*/
char **allocC2(int a, int b){
  int i;
  char **tmp=(char**)malloc(a*sizeof(char*));
  tmp[0]=(char*)malloc(a*b*sizeof(char));
  for(i=1;i<a;i++) tmp[i]=tmp[i-1]+b;
  return tmp;
}
/*======================================================================*/
void deallocC2(char **a){
  free(a[0]);
  free(a);
}
/*======================================================================*/
char ***allocC3(int a, int b, int c){
  int i,ab=a*b;
  char ***tmp=(char***)malloc(a*sizeof(char**));
  tmp[0]=(char**)malloc(ab*sizeof(char*));
  for(i=1;i<a;i++) tmp[i]=tmp[i-1]+b;
  tmp[0][0]=(char*)malloc(ab*c*sizeof(char));
  for(i=1;i<ab;i++) tmp[0][i]=tmp[0][i-1]+c;
  return tmp;
}
/*======================================================================*/
void deallocC3(char ***a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}
/*======================================================================*/
char **allocsymC2(int a){ /*only (i,j) with i<j is filled*/
  int i;
  int n=a*(a-1)/2;
  char **tmp=(char**)malloc((a-1)*sizeof(char*));
  tmp[0]=(char*)malloc(n*sizeof(char));
  for(i=1;i<a-1;i++) tmp[i]=tmp[i-1]+a-i;
  return tmp;
}
/*======================================================================*/
void deallocsymC2(int **a){
  free(a[0]);
  free(a);
}
/*======================================================================*/
