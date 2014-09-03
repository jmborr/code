#ifndef _MY_ALLOCATE_
#define _MY_ALLOCATE_

template<class T>
T **alloc2(int i,int j, T **x=NULL){
  T **pt=new T*[i];  pt[0]=new T[i*j];
  for(int k=1;k<i;k++) pt[k]=pt[k-1]+j;
  if(x!=NULL){int a=i*j;  for(int k=0;k<a;k++) pt[0][k]=x[0][k]; }/*initialize*/ 
  return pt;
}
/*==========================================================*/
/*NOTE: use ALLWAYS dealloc2Half to de-allocate the memory*/
/*pt[i][j] with i<j*/
template<class T>
T **alloc2Half(int n){
  T **pt=new T*[n-1];  pt[0]=new T[n*(n-1)/2];
  int j=n-2;
  pt[0]--;  /*so that pt[0][1] is first allocated memory element*/
  for(int k=1;k<n-1;k++){ pt[k]=pt[k-1]+j; j--; } /*so that pt[i][i+delta] newly mem*/
  return pt;
}
/*==========================================================*/
template<class T>
T **dealloc2Half(T **pt){
  pt[0]++; /*before removing pt[0], we have to undo the shift!!! */
  delete [] pt[0];
  delete [] pt;
}
/*==========================================================*/
template<class T>
void alloc2ref(T **&pt, int i,int j){
  pt=new T*[i];  pt[0]=new T[i*j];
  for(int k=1;k<i;k++) pt[k]=pt[k-1]+j;
}
/*==========================================================*/
template<class T>
void dealloc2(T **&pt){ free(pt[0]); free(pt); }
/*==========================================================*/
template<class T>
T *** alloc3(int a, int b, int c){
  int i,ab=a*b;
  T ***pt;
  pt=new T**[a];
  pt[0]=new T*[ab];
  for(i=1;i<a;i++) pt[i]=pt[i-1]+b;
  pt[0][0]=new T[ab*c];
  for(i=1;i<ab;i++) pt[0][i]=pt[0][i-1]+c;
  return pt;
}
/*======================================================================*/
template<class T>
void alloc3ref(T ***&pt, int a, int b, int c){
  int i,ab=a*b;
  pt=new T**[a];   pt[0]=new T*[ab];
  for(i=1;i<a;i++) pt[i]=pt[i-1]+b;
  pt[0][0]=new T[ab*c];
  for(i=1;i<ab;i++) pt[0][i]=pt[0][i-1]+c;
}
/*======================================================================*/
template<class T>
void dealloc3(T ***&pt){ free(pt[0][0]); free(pt[0]); free(pt); }
/*======================================================================*/

#endif /*_MY_ALLOCATE_*/