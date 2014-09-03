#include<iostream.h>
using namespace std;
#include<stdio.h>
#include<iomanip.h>
#include<string.h>
#include<stdlib.h>
#include<fstream.h>
#include<math.h>
#include "pdbClasses2.h"
/*==============================================================*/
/*Some functions that deal with generalized atom types*/
string gat2str(gen_atom_type type){
  string gat;
  switch(type){
  case _CA_: gat.assign(" CA ");break;
  case _CB_: gat.assign(" CB ");break;
  case _HV_: gat.assign(" HV ");break;
  case _HN_: gat.assign(" H  ");break;
  case _METHYLC_: gat.assign(" MC ");break;
  }
  return gat;
}
/*==============================================================*/
gen_atom_type str2gat(string &type){
  if(type==" CA "||type=="_CA_") return _CA_;
  else if(type==" CB "||type=="_CB_") return _CB_;
  else if(type==" HV "||type=="_HV_") return _HV_;
  else if(type==" MC "||type=="_METHYLC_") return _METHYLC_;
}
/*==============================================================*/
/*Some functions that deal with contact selections*/
string cs2str(cont_sel &type){
  string cs;
  switch(type){
  case _CA_CA_: cs.assign(" CA "); break;
  case _CB_CB_: cs.assign(" CB "); break;
  case _HV_HV_: cs.assign(" HV "); break;
  case _HVSC_HVSC_: cs.assign(" SS "); break;/*S, SC stand for sidechain*/
  case _HVBK_HVBK_: cs.assign(" BB "); break;
  case _HVBK_HVSC_: cs.assign(" BS "); break;
  case _F_HVBK_S_HVSC_: cs.assign("FBSS"); break;
  case _F_HVSC_S_HVBK_: cs.assign("FSSB"); break;
  case _HN_HN_: cs.assign(" HH "); break;
  case _HN_METHYLC_: cs.assign("HNMC"); break;
  case _METHYLC_HN_: cs.assign("MCHN"); break;
  case _METHYLC_METHYLC_: cs.assign("MCMC"); break;
  default:   cs.assign("NONE"); break;
  }
  return cs;
}
/*==============================================================*/
cont_sel str2cs(string &type){
  string cs;
  if(type==" CA " || type=="_CA_CA_") return _CA_CA_;
  else if(type==" CB " || type=="_CB_CB_") return _CB_CB_;
  else if(type==" HV " || type=="_HV_HV_") return _HV_HV_;
  else if(type==" SS " || type=="HVSC_HVSC" || type=="_HVSC_HVSC_")
    return _HVSC_HVSC_;
  else if(type==" BB " || type=="HVBK_HVBK" || type=="_HVBK_HVBK_")
    return _HVBK_HVBK_;
  else if(type==" BS " || type=="HVBK_HVSC" || type=="_HVBK_HVSC_")
    return _HVBK_HVSC_;
  else if(type=="FBSS" || type=="F_HVBK_S_HVSC" || type=="_F_HVBK_S_HVSC_")
    return _F_HVBK_S_HVSC_;
  else if(type=="FSSB" || type=="F_HVSC_S_HVBK" || type=="_F_HVSC_S_HVBK_")
    return _F_HVSC_S_HVBK_;
  else if(type==" HH " || type=="HN_HN" || type=="_HN_HN_")
    return _HN_HN_;
  else if(type=="HNMC" || type=="HN_METHYLC_" || type=="_HN_METHYLC_")
    return _HN_METHYLC_;
  else if(type=="MCHN" || type=="METHYLC_HN" || type=="_METHYLC_HN_")
    return _METHYLC_HN_;
  else if(type=="MCMC" || type=="METHYLC_METHYLC" || type=="_METHYLC_METHYLC_")
    return _METHYLC_METHYLC_;
 
}
/*==============================================================*/
bool symmetric_map(cont_sel &type){
  switch(type){
  case _CA_CA_:
  case _CB_CB_:
  case _HV_HV_:
  case _HVSC_HVSC_:
  case _HVBK_HVBK_:
  case _HVBK_HVSC_: 
  case _HN_HN_:
  case _METHYLC_METHYLC_: return true;
  default: return false;
  }
}
/*==============================================================*/
void print_symmetric_map(ostream &out, double **map, int &N){
  char buff[32];
  for(int i=0;i<N-1;i++)
    for(int j=i+1;j<N;j++)
      if(map[i][j]>0.){
	sprintf(buff,"%10d%10d\n",i+1,j+1);
	out<<buff;
      }
}
/*==============================================================*/
void print_asymmetric_map(ostream &out, double **map, int &N){
  char buff[32];
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      if(map[i][j]>0.){
	sprintf(buff,"%10d%10d\n",i+1,j+1);
	out<<buff;
      }
}
/*==============================================================*/
/*Some functions to zip and unzip files*/
zip_format returnZipFormat(const char *file){
  if(strstr(file,".Z"))
    return _Z_;
  else if(strstr(file,".gz"))
    return _gz_;
  else if(strstr(file,".bz2"))
    return _bz2_;
  else if(strstr(file,".zip"))
    return _zip_;
  else 
    return _NONE_;
}
/*==============================================================*/
string zf2str(zip_format fmt){
  string zf;
  switch(fmt){
  case _Z_:  zf.assign("_Z_");   break;
  case _gz_: zf.assign("_gz_");  break;
  case _bz2_:zf.assign("_bz2_"); break;
  default:   zf.assign("_NONE_");   break;
  }
  return zf;
}
/*==============================================================*/
/*unzips de file (if zipped) and returns the name of the file*/
/*without the zip extension*/
bool unzip(char *uzfile,const char *file){
  /*cout<<"bool unzip(...)\n";*/
  zip_format fmt=returnZipFormat(file); /*cout<<zf2str(fmt)<<endl;*/
  bool r;
  char command[256],util[64];
  switch(fmt){
  case _Z_:
  case _gz_:
    strcpy(command,"gunzip ");strcat(command,file);
    strcpy(util,"gunzip"); r=true;
    break;
  case _bz2_:
    strcpy(command,"bunzip2 ");strcat(command,file);
    strcpy(util,"bunzip2"); r=true;
    break;
  default:
    r=false;
    break;
  }
  if(r){
    /*cout<<"command="<<command<<endl;*/
    int n=(int)(system(command));
    if(n!=0){/*Exception handling*/
      cout<<"ERROR: from zip function:\n";
      if(n==-1)cout<<"       on \"system\" call\n";
      else cout<<"       "<<util<<": command not found\n";
      exit(0);
    }
    /*remove zip extension from filename*/
    char dot='.';
    char *pdot=strrchr(file,(int)(dot)); 
    n=(int)(pdot-file);
    strncpy(uzfile,file,n);
  }
  else
    strcpy(uzfile,file);
  return r;
}
/*==============================================================*/
/*zips the file (if unzipped) and returns the name of the file*/
/*with the desired zip extension*/
bool zip(char *zfile,const char *uzfile,zip_format &fmt){
  bool r;
  char command[256],util[64];
  strcpy(zfile,uzfile);

  /*check if already zipped*/
  zip_format fmt2=returnZipFormat(uzfile);
  if(fmt2!=_NONE_){/*file is already zipped*/
    if(fmt2!=fmt) /*different zipped format than desired*/
      cout<<"WARNING: "<<uzfile<<" already zipped with a different zip utility than desired\n" ;
    return true;
  }

  switch(fmt){
  case _Z_:
    strcpy(command,"compress ");strcat(command,uzfile);
    strcat(zfile,".Z"); strcpy(util,"compress"); r=true;
    break;
  case _gz_:
    strcpy(command,"gzip --best ");strcat(command,uzfile);
    strcat(zfile,".gz"); strcpy(util,"gzip"); r=true;
    break;
  case _bz2_:
    strcpy(command,"bzip2 --best ");strcat(command,uzfile);
    strcat(zfile,".bz2"); strcpy(util,"bzip2"); r=true;
    break;
  default:

    r=false;
  }
  if(r){
    int n=(int)(system(command));
    if(n!=0){/*Exception handling*/
      cout<<"ERROR: from zip function:\n";
      if(n==-1)cout<<"       on \"system\" call\n";
      else cout<<"       "<<util<<": command not found\n";
      exit(0);
    }
  }
  return r;
}
/*==============================================================*/
/*zips the file (if unzipped) with the desired zip extension*/
bool zip(const char *uzfile,zip_format &fmt){
  bool r;
  char command[256],util[64];

  /*check if already zipped*/
  zip_format fmt2=returnZipFormat(uzfile);
  if(fmt2!=_NONE_){/*file is already zipped*/
    if(fmt2!=fmt) /*different zipped format than desired*/
      cout<<"WARNING: "<<uzfile<<" already zipped with a different zip utility than desired\n" ;
    return true;
  }

  switch(fmt){
  case _Z_:
    strcpy(command,"compress ");strcat(command,uzfile);
    strcpy(util,"compress"); r=true;
    break;
  case _gz_:
    strcpy(command,"gzip --best ");strcat(command,uzfile);
    strcpy(util,"gzip"); r=true;
    break;
  case _bz2_:
    strcpy(command,"bzip2 --best ");strcat(command,uzfile);
    strcpy(util,"bzip2"); r=true;
    break;
  default:

    r=false;
  }
  if(r){
    int n=(int)(system(command));
    if(n!=0){/*Exception handling*/
      cout<<"ERROR: from zip function:\n";
      if(n==-1)cout<<"       on \"system\" call\n";
      else cout<<"       "<<util<<": command not found\n";
      exit(0);
    }
  }
  return r;
}
/*==============================================================*/
/*Some funtions to convert strings to numbers and viceversa*/

void int2string5( const int & number , string &target ){
  char tmpChar[6] ; //5 plus null terminator character
  sprintf( tmpChar , "%5d" , number ) ;
  string tmpString( tmpChar ) ;
  target = tmpString ;
}

void int2string4( const int & number , string &target ){
  char tmpChar[5] ; //4 plus null terminator character
  sprintf( tmpChar , "%4d" , number ) ;
  string tmpString( tmpChar ) ;
  target = tmpString ;
}

void string2int( const string &theString ,  int &number ){
  const char *pt = theString.c_str() ;
  number  = atoi( pt ) ;
}

void double2string8_3( const double &number , string &target){
  char tmpChar[9] ; //8 plus null terminator character
  sprintf( tmpChar , "%8.3lf" , number ) ;
  string tmpString( tmpChar ) ;
  target = tmpString ;
}

void double2string6_2( const double &number , string &target){
  char tmpChar[7] ; //6 plus null terminator character
  sprintf( tmpChar , "%6.2lf" , number ) ;
  string tmpString( tmpChar ) ;
  target = tmpString ;
}

void string2double( const string &theString , double &number){
  const char *pt = theString.c_str() ;
  number = atof( pt ) ;
}

//----some handy functions-------------------
bool is_ATOM( string &theString){
  string check = "ATOM ";
  if( check.compare( theString.substr(0,5) ) == 0 ){ return true ; }
  else{ return false ; }
}

bool is_ANISOU( string &theString){
  string check = "ANISOU " ;
  if( check.compare( theString.substr(0,7) ) == 0 ){ return true ; }
  else{ return false ; }
}
/*==================================================================*/
bool is_TER( string &theString){
  string check = "TER";
  if( check.compare( theString.substr(0,3) ) == 0 ){ return true ; }
  else{ return false ; }
}
/*==================================================================*/
bool is_END( string &theString){
  string check = "END";
  if( check.compare( theString.substr(0,3) ) == 0 ){ return true ; }
  else{ return false ; }
}
/*==================================================================*/
bool is_END( ifstream &PDB ){
  /*cout<<"function bool is_END( ifstream & )\n" ;*/
  int position = PDB.tellg( ) ; 
  string line ; /*cout<<"line="<<line<<endl;*/
  bool flag=false ;
  getline( PDB, line ) ;
  if( is_END( line ) ){ flag = true ; }
  PDB.seekg( position ) ;
  return flag ;
}
/*==================================================================*/
bool is_HETATM( string &theString ){
  if(theString.find("HETATM") == string::npos) return false ;
  return true ;
}
/*==================================================================*/
bool is_SIGATM( string &theString ){
  if(theString.find("SIGATM") == string::npos) return false ;
  return true ;
}

bool is_AMINOACID( string &theString ){
  string query=theString.substr(17,3) ;
  if( all_aa.find(query) == string::npos ) return false ;
  return true ;
}
/*=================================================================*/
string three_letter_code_to_one_letter_code( const string &name){
  string tmp("-") ;
  if( name=="GLY" ){ tmp.assign("G") ; }
  else if( name=="ALA" ){ tmp.assign("A") ; }
  else if( name=="VAL" ){ tmp.assign("V") ; }
  else if( name=="LEU" ){ tmp.assign("L") ; }
  else if( name=="ILE" ){ tmp.assign("I") ; }
  else if( name=="SER" ){ tmp.assign("S") ; }
  else if( name=="THR" ){ tmp.assign("T") ; }
  else if( name=="CYS" ){ tmp.assign("C") ; }
  else if( name=="MET" ){ tmp.assign("M") ; }
  else if( name=="PRO" ){ tmp.assign("P") ; }
  else if( name=="ASP" ){ tmp.assign("D") ; }
  else if( name=="ASN" ){ tmp.assign("N") ; }
  else if( name=="GLU" ){ tmp.assign("E") ; }
  else if( name=="GLN" ){ tmp.assign("Q") ; }
  else if( name=="LYS" ){ tmp.assign("K") ; }
  else if( name=="ARG" ){ tmp.assign("R") ; }
  else if( name=="HIS" ){ tmp.assign("H") ; }
  else if( name=="PHE" ){ tmp.assign("F") ; }
  else if( name=="TYR" ){ tmp.assign("Y") ; }
  else if( name=="TRP" ){ tmp.assign("W") ; }
  return tmp;
}
/*=================================================================*/
/*from bigger to smaller*/
void bubble_sort( double r[], const int &n ){
  int k ;
  double tmp ;
  for( int j=1; j<n; j++ ){
    k=j ;
    while( r[k-1]<r[k] && k>0 ){
      tmp=r[k-1]; r[k-1]=r[k]; r[k]=tmp ;
      k--;
    }
  }
}
/*=================================================================*/
/*from bigger to smaller, reordering stored in index[]*/
void bubble_sort( double r[], const int &n, int index[] ){
  int k, i_tmp ;
  double tmp ;
  for( int j=0; j<n; j++ ){ index[j]=j ; }
  for( int j=1; j<n; j++ ){
    k=j ;
    while( r[k-1]<r[k] && k>0 ){
      tmp=r[k-1]; r[k-1]=r[k]; r[k]=tmp ;
      i_tmp=index[k-1]; index[k-1]=index[k]; index[k]=i_tmp ;
      k--;
    }
  }
}
/*=================================================================*/
void dump_array( double **r, int n_rows, int n_cols ){
  for(int i=0; i<n_rows; i++){
      for(int j=0; j<n_cols; j++){
	cout << r[i][j] << "  " ;
      }
      cout << endl ;
  }
}
/*=================================================================*/
void dump_array2( double **r, int n_rows, int n_cols ){
  for(int i=0; i<n_rows; i++)
    for(int j=0; j<n_cols; j++)
      cout <<i<<" "<<j<<" "<<r[i][j]<<"\n";
}
/*=======================================================================*/
double **alloc_array( int n_rows, int n_cols ){
  double **ptr = new double* [n_rows] ;
  ptr[0] =  new double[n_rows*n_cols] ;
  for( int i=1; i<n_rows; i++ ){ ptr[i] = ptr[ i-1 ] + n_cols ; }
  return ptr ;
}
/*=======================================================================*/
void init_array(const double &x, double **A, int &n){
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      A[i][j]=x;
}
/*=======================================================================*/
void assign_array(double **B,double **A,int &n){
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      B[i][j]=A[i][j];
}
/*=======================================================================*/
void mult_arrayBA(double **B, double **A, int &n){
  double **C=alloc_array(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      C[i][j]=0.;
      for(int k=0;k<n;k++)
	C[i][j]+=B[i][k]*A[k][j];
    }
  }
  assign_array(B,C,n);
}
/*=======================================================================*/
void mult_arrayCBA(double **C, double **B, double **A, int &n){
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      C[i][j]=0.;
      for(int k=0;k<n;k++)
	C[i][j]+=B[i][k]*A[k][j];
    }
  }
}
/*=======================================================================*/
int **alloc_int_array( int n_rows, int n_cols ){
  int **ptr = new int* [n_rows] ;
  ptr[0] =  new int[n_rows*n_cols] ;
  for( int i=1; i<n_rows; i++ ){ ptr[i] = ptr[ i-1 ] + n_cols ; }
  return ptr ;
}
/*=======================================================================*/
void init_int_array(int &x, int **A, int &n){
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      A[i][j]=x;
}
/*=======================================================================*/
void assign_array(int **B,int **A,int &n){
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      B[i][j]=A[i][j];
}
/*=======================================================================*/
void mult_int_arrayBA(int **B, int **A, int &n){
  int **C=alloc_int_array(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      C[i][j]=0;
      for(int k=0;k<n;k++)
	C[i][j]+=B[i][k]*A[k][j];
    }
  }
  assign_array(B,C,n);
}
/*=======================================================================*/
void mult_int_arrayCBA(int **C, int **B, int **A, int &n){
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      C[i][j]=0;
      for(int k=0;k<n;k++)
	C[i][j]+=B[i][k]*A[k][j];
    }
  }
}
/*=======================================================================*/
double det_3x3( double **m ){
  return m[0][0]*( m[1][1]*m[2][2] - m[1][2]*m[2][1] ) -
    m[1][0]*( m[0][1]*m[2][2] - m[0][2]*m[2][1] ) +
    m[2][0]*( m[0][1]*m[1][2] - m[0][2]*m[1][1] ) ;
}
/*===========================================================*/
void remove_bound_cond(double **r,const int &n,const double &box){
  double prev[3],dr[3],hb=box/2.;
  for(int i=0;i<3;i++){prev[i]=r[0][i];}
  for(int j=1;j<n;j++){ for(int i=0;i<3;i++){
      dr[i]=r[j][i]-prev[i];
      if(dr[i]>=hb){dr[i]-=box;}  else if(dr[i]<=-hb){dr[i]+=box;}
      prev[i]=r[j][i];   r[j][i]=r[j-1][i]+dr[i];
  }  }
}
/*=======================================================================
  store the rotation matrix, rot_matrix, that is to be performed on the set
  of x coordinates so that rmsd between x and y is minimized. We assume that
  the center of mass of x and y coincide.
         |x'[0]|     |u[0][0] ... u[0][2]|  |x[0]|
	 | .   |     |  .           .    |  | .  |
	 | .   | =   |  .           .    |  | .  |       
	 | .   |     |  .           .    |  | .  |       
         |x'[0]|     |u[0][0] ... u[0][2]|  |x[0]|

*/
double get_rot_matrix_by_rms( double **x, double **y, long nn,
			      double **rot_matrix ){
  double r[3][3], rr[3][3], mu[3], p[4],q[4];
  double aa[3][3], b[3][3], u[3][3];
  int i,j,n,m,k;
  double sum=0;
  double *dr = new double[nn];
  double detR ;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      double  rij=0;
      for(k=0;k<nn;k++)
	{
	  rij+=y[k][i]*x[k][j];
	}
      r[i][j]=rij;
    }

   detR =r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1]);
   detR-=r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0]);
   detR+=r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	double  rij=0;
	for(k=0;k<3;k++)
	  {
	    rij+=r[k][i]*r[k][j];
	  }
	rr[i][j]=rij;
      }
  /* charcteristic third order polinomial p(mu)= det|rr -mu^2I| */
  
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
    for(i=0;i<3;i++)
      sum+=x[k][i]*x[k][i]+ y[k][i]*y[k][i]; 
  
  sum*=0.5;
  
  for(i=0;i<3;i++)
    {
      mu[i]=sqrt(p[i+1]);
      sum-=mu[i];
    }   

  /*find minimum of mu[] and change sign*/
  i=1;
  if(detR<0)
    {
      for(j=2;j<=3;j++)
	  if(p[j]<p[i])i=j;
      sum+=2*mu[i-1];
      mu[i-1]=-mu[i-1];  
  }

  /*rms is equal to the half of the sum of the squares of coordinates 
    y and x minus sum of kabash eigenvalues mu 
  */
  
  
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

  for(i=0;i<3;i++)/*calculate u[][]*/
    {
      for(j=0;j<3;j++)
	{
	  double uij=0;
	  for(k=0;k<3;k++)
	    uij+=b[i][k]*aa[j][k];
	  u[i][j]=uij;
	}
    }
  
  
  /* u is rotatation matrix : x'=ux */
  /* such that rms(x'-y) is minimal */
  for(int i=0; i<3; i++){ for(int j=0;j<3;j++){ rot_matrix[i][j] = u[i][j]; } }

  sum = sqrt( fabs( sum * 2.0 / (double)nn ) ) ;
  return sum;
}
/*=======================================================================
  The tranformation of x coordinates so that rmsd between x and y is minimized
  is stored in z coordinates. We assume that the center of mass of x and y
  coincide.*/
double get_rmsII(double **x,double ** y, double ** z, long nn)
{
  double r[3][3], rr[3][3], mu[3], p[4],q[4];
  double aa[3][3], b[3][3], u[3][3];
  int i,j,n,m,k;
  double sum=0;
  double *dr = new double[nn];
  double detR;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
    {
      double  rij=0;
      for(k=0;k<nn;k++)
	{
	  rij+=y[k][i]*x[k][j];
	}
      r[i][j]=rij;
    }

 /*computing the determinant of matrix R~R */
 detR =r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1]);
 detR-=r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0]);
 detR+=r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	double  rij=0;
	for(k=0;k<3;k++)
	  {
	    rij+=r[k][i]*r[k][j];
	  }
	rr[i][j]=rij;
      }
  /* charcteristic third order polinomial p(mu)= det|rr -mu^2I| */
  
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
    double z=-p[1];  /*z = min( |p[i]| ), i=0..2*/
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
 
 /*rms is equal to the half of the sum of the squares of coordinates y[][] and
   x[][], minus sum of kabash eigenvalues mu[]  */
  for(k=0;k<nn;k++)
    for(i=0;i<3;i++)
      sum+=x[k][i]*x[k][i]+ y[k][i]*y[k][i]; 
  sum*=0.5;

  for(i=0;i<3;i++){ mu[i]=sqrt(p[i+1]);  sum-=mu[i]; }
  
  i=1;
  if(detR<0){ /*the minimum rms includes inversion. Let's avoid this*/
    /*find minimum of mu[] and change sign. [] stores square values of mu[]*/
    for(j=2;j<=3;j++){ if(p[j]<p[i]){ i = j ; } }
    sum+=2*mu[i-1]; mu[i-1]=-mu[i-1];  
  }
  
  for(k=0;k<3;k++) /*obtain eigenvalues a[][]*/
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
  
  for(k=0;k<3;k++) /*eigenvectors aa[][k] are orthonormal*/
    {
      double ak=0;
      for(j=0;j<3;j++)
	ak+=aa[j][k]*aa[j][k];
      ak=1/sqrt(ak);
      for(j=0;j<3;j++)
	aa[j][k]*=ak;
    }
  
  
  for(k=0;k<3;k++) /*obtain b[][k]=r[][]*aa[][k]*/
    {
      for(i=0;i<3;i++)
	{
	  double bki=0;
	  for(j=0;j<3;j++)
	    bki+=r[i][j]*aa[j][k];
	  b[i][k]=bki/mu[k];/*if detR<0, then there's one negative mu[]*/
	}
    }
  
  
  for(i=0;i<3;i++) /*orthogonal matrix u[][]=b[][]*aa[][]*/
    {
      for(j=0;j<3;j++)
	{
	  double uij=0;
	  for(k=0;k<3;k++)
	    uij+=b[i][k]*aa[j][k];
	  u[i][j]=uij;
	}
    }
  
  
  /* is rotatation matrix : x'=ux */
  /* such that rms(x'-y) is minimal */
  for(k=0;k<nn;k++)
    {
      for(j=0;j<3;j++)
	{
	  double ykj=0;
	  for(i=0;i<3;i++)
	    ykj+=u[j][i]*x[k][i];
	  z[k][j] = ykj ; //store transformed coordinates.
	  p[j]=ykj;
	  dr[k]+=(ykj-y[k][j])*(ykj-y[k][j]);
	}
    }
  /* dr[k] accumulates rms for k-yh atom from many calls of
     this function  the rotated coordinates x are ykj , k
     goes from 0 to nn-1, j is 0(x), 1(y), or 2(z). */
  sum = sqrt( fabs( sum * 2.0 / (double)nn ) ) ;
  return sum;
}
/*=======================================================================
  /*We assume that the center of mass of x and y coincide*/
double get_rms(double **x,double ** y,long nn)
{
  double r[3][3], rr[3][3], mu[3], p[4],q[4];
  int i,j,n,m,k;
  double sum=0;
  double detR;

  for(i=0;i<3;i++)    /*fill matrix r[][]*/
    for(j=0;j<3;j++)
    {
      double  rij=0;
      for(k=0;k<nn;k++)
	{
	  rij+=y[k][i]*x[k][j];
	}
      r[i][j]=rij;
    }

  detR =r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1]);
  detR-=r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0]);
  detR+=r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]);

  for(i=0;i<3;i++)   /*fill matrix rr[][] = r[][]^t * r[][] */
    for(j=0;j<3;j++)
      {
	double  rij=0;
	for(k=0;k<3;k++)
	  {
	    rij+=r[k][i]*r[k][j];
	  }
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

  for(k=0;k<nn;k++)
    for(i=0;i<3;i++)
      sum+=x[k][i]*x[k][i]+ y[k][i]*y[k][i]; 
  
  sum*=0.5;
  
  for(i=0;i<3;i++)
    {
      mu[i]=sqrt(p[i+1]);
      sum-=mu[i];
    }    

  i=1;
  if(detR<0)
    {
      for(j=2;j<=3;j++)/*find minimum of mu[] and change sign*/
	  if(p[j]<p[i])i=j;
      sum+=2*mu[i-1]; /*correct sum*/
      mu[i-1]=-mu[i-1];  /*change sign of minimum mu[]*/
  }

  sum = sqrt( fabs( sum * 2.0 / (double)nn ) ) ;
  return sum;
}


//------------ class PDBvector ------------- 
string  &operator<<(string  &output, const PDBvector &coords){
  string tmp ; double2string8_3( coords.x , tmp ) ;  output = tmp ;
  double2string8_3( coords.y , tmp ) ; output += tmp ;
  double2string8_3( coords.z , tmp ) ; output += tmp ;

  return output ;
}
/*============================================================*/
ostream &operator<<(ostream &OUT, const PDBvector &coords){
  string tmp ;  tmp << coords ;  OUT << tmp ;
  return OUT ;
}
/*============================================================*/
PDBvector::PDBvector( ){ x = 0.0 ; y=0.0 ; z = 0.0 ; }
/*============================================================*/
PDBvector::PDBvector( double const &xprime, double  const &yprime, 
		      double  const &zprime){
  x = xprime ; y = yprime ;z = zprime ;
}
/*============================================================*/
PDBvector::PDBvector( double *coords ){
  x=coords[0]; y=coords[1]; z=coords[2];
}
/*============================================================*/
void PDBvector::assign( double const & xprime, double const & yprime,
			double const & zprime){
  x = xprime ; y = yprime ;z = zprime ;
}
/*============================================================*/
void PDBvector::assign( double *coord ){
  x = coord[0] ; y = coord[1] ; z=  coord[2] ;
}
/*============================================================*/
double &PDBvector::operator[]( const int &index ){
  if(index<0 || index>2 ){ 
    cout<<"\nERROR: a PDBvector can't have index<0||index>2\n";
    exit(1) ;
  }
  if(index==0){ return x ; }
  else if(index==1){ return y ; }
  else{ return z ; }
}
/*============================================================*/
double PDBvector::getComponent( const int &index ){
  double zero=0.0;
  if(index==0) return x;
  else if(index==1) return y;
  else if(index==2) return z;
  else return zero;
}
/*============================================================*/
double *PDBvector::getComponents( ){
  double *tmp = new double[3] ;
  tmp[0] = x ; tmp[1] = y ; tmp[2] = z ; 
  return tmp ;
}
/*============================================================*/
void PDBvector::getComponentsII( double *r ){
  r[0] = x ; r[1] = y ; r[2] = z ; 
}
/*============================================================*/
PDBvector &PDBvector::operator=( const PDBvector &theVector ){
  x = theVector.x ;  y = theVector.y ;  z = theVector.z ;
  return *this ;
}
/*============================================================*/
PDBvector &PDBvector::operator+=( const PDBvector &theVector ){
  x += theVector.x ;  y += theVector.y ;  z += theVector.z ;
  return *this ;
}
/*============================================================*/
PDBvector PDBvector::operator+( const PDBvector &right ){
  PDBvector tmp;
  tmp.x=x+right.x; tmp.y=y+right.y; tmp.z=z+right.z;
  return tmp;
}
/*============================================================*/  
PDBvector &PDBvector::operator-=( const PDBvector &theVector ){
  x -= theVector.x ;  y -= theVector.y ;  z -= theVector.z ;
  return *this ;
}
/*============================================================*/
PDBvector PDBvector::operator-( const PDBvector &right ){
  PDBvector tmp;
  tmp.x=x-right.x; tmp.y=y-right.y; tmp.z=z-right.z;
  return tmp;
}
/*============================================================*/
PDBvector &PDBvector::operator/=( const double &number ){
  if( number != 0.0 ){ x /= number ; y /= number ; z /= number ; }
  return *this ;
}
/*============================================================*/
PDBvector PDBvector::operator/( const double &number ){
  if(number != 0){
    PDBvector tmp;
    tmp.x = x/number; tmp.y=y/number; tmp.z=z/number;
    return tmp;
  }
  else{
    cout<<"ERROR MESSAGE in 'PDBvector PDBvector::operator/( const double &number )'\nYou cannot divide by zero\n";
    exit(1);
  }
}
/*============================================================*/
PDBvector &PDBvector::operator*=( const double &number ){
  x *= number ; y *= number ; z *= number ; 
  return *this ;
}
/*============================================================*/
PDBvector PDBvector::operator*( const double &number ){
  PDBvector tmp;
  tmp.x = x*number; tmp.y = y*number; tmp.z = z*number;
  return tmp;
}
/*============================================================*/
double PDBvector::operator*( const PDBvector &v2 ){
  return x*v2.x + y*v2.y + z*v2.z ;
}
/*============================================================*/
double PDBvector::d2( const PDBvector &theVector ){
  double d2=(x-theVector.x)*(x-theVector.x) +
         (y-theVector.y)*(y-theVector.y) +
         (z-theVector.z)*(z-theVector.z) ;
  //  cout << x <<" "<< y <<" "<< z <<" "<< theVector.x <<" "<< theVector.y
  //     <<" "<< theVector.z <<" "<< d2 << endl ;
  return d2 ;
}
/*============================================================*/
double PDBvector::norm2( ){ return x*x+y*y+z*z ; }
/*============================================================*/
void PDBvector::normalize( ){ 
  double norm2 = this->norm2( ) ;
  if( norm2 > 0 ){ (*this) /= sqrt( norm2 ) ; }
}
/*============================================================*/
double PDBvector::angle_with( PDBvector &v2){
  return acos( (*this)*v2 / sqrt( this->norm2() * v2.norm2() ) ) ;
}
/*============================================================*/
PDBvector PDBvector::operator^( const PDBvector &right){
  PDBvector tmp;
  tmp.x = y*right.z - z*right.y; 
  tmp.y = z*right.x - x*right.z;
  tmp.z = x*right.y - y*right.x;
  return tmp;
}
/*============================================================*/
void PDBvector::rotate( const PDBvector &axis, const double &angle){
  double phi = 3.1415927*angle/180.0;
  PDBvector n = axis;
  n.normalize( );
  *this=(*this)*cos(phi)+ n*((n*(*this))*(1-cos(phi))) - ((*this)^n)*sin(phi);
}
/*============================================================*/
PDBvector PDBvector::rotateII( const PDBvector &axis, const double &angle ){
  PDBvector tmp = *this;
  tmp.rotate( axis, angle ) ;
  return tmp;  
}
/*============================================================*/
void PDBvector::rotate( const PDBvector &new_origin, const PDBvector &axis,
			const double &angle){
  PDBvector O=new_origin;
  *this =  O + (*this - O).rotateII( axis, angle ) ;
}
/*============================================================*/
PDBvector PDBvector::rotateII( const PDBvector &new_origin, 
			       const PDBvector &axis, const double &angle){
  PDBvector tmp = *this;
  tmp.rotateII( new_origin, axis, angle ) ;
  return tmp;  
}
/*============================================================*/
void PDBvector::rotate( double ** rotation_matrix ){
  double zero=0.0 ;
  PDBvector tmp( zero, zero, zero ) ;
  tmp.x = rotation_matrix[0][0] * x + rotation_matrix[0][1] * y 
    + rotation_matrix[0][2] * z ;
  tmp.y = rotation_matrix[1][0] * x + rotation_matrix[1][1] * y 
    + rotation_matrix[1][2] * z ;
  tmp.z = rotation_matrix[2][0] * x + rotation_matrix[2][1] * y 
    + rotation_matrix[2][2] * z ;
  *this = tmp ;
}
/*============================================================*/
PDBvector PDBvector::rotateII( double ** rotation_matrix ){
  double zero=0.0 ;
  PDBvector tmp(  zero, zero, zero ) ;
  tmp.x = rotation_matrix[0][0] * this->x + rotation_matrix[0][1] * this->y 
    + rotation_matrix[0][2] * this->z ;
  tmp.y = rotation_matrix[1][0] * this->x + rotation_matrix[1][1] * this->y 
    + rotation_matrix[1][2] * this->z ;
  tmp.z = rotation_matrix[2][0] * this->x + rotation_matrix[2][1] * this->y 
    + rotation_matrix[2][2] * this->z ;
  return tmp;
}
/*============================================================
/*radius of the circle that goes through the three vectors*/
double PDBvector::radius(PDBvector&left, PDBvector&right){
  PDBvector A,B,O=*this,C;
  A=left-O; B=right-O; C=A-B;
  double A2=A.norm2(),B2=B.norm2(),C2=C.norm2(),AB=A*B;
  return sqrt(0.5*(A2*B2*C2)/(A2*B2-AB*AB));
}
/*============================================================*/

//----------------------------------------------------------------
//check if the string has ATOM information
bool isPDBatom( string &theString ){
  string tmpString;
  tmpString.assign( theString , 0 , 6 ) ;
  if( tmpString.find( "ATOM  " ) == string::npos ){ return false ; }
  else{ return true ; }
}

//----------------------------------------------------------------

//------------ class PDBatom  --------------
string  &operator<<(string  &output, const PDBatom &PDBatomprime){
  output.assign("ATOM  "); 
  string tmp ; 
  int2string5( PDBatomprime.serial , tmp ) ; output += tmp ;
  output += " " ;
  output += PDBatomprime.name ; output += PDBatomprime.altLoc ;
  output += PDBatomprime.resName ; output += " " ; output += PDBatomprime.chainID ;
  int2string4( PDBatomprime.resSeq , tmp )  ; output += tmp ;
  output += PDBatomprime.iCode ; output += "   " ;
  tmp << PDBatomprime.r ; output +=  tmp ;
  double2string6_2( PDBatomprime.occupancy , tmp ) ; output += tmp ;
  double2string6_2( PDBatomprime.tempFactor , tmp ) ; output += tmp ;
  output += "      " ; output += PDBatomprime.segID ;
  output += PDBatomprime.element ; output += PDBatomprime.charge ;

  return output ;
}

ostream &operator<<( ostream &OUT, const PDBatom &thePDBatom ){
  string tmp ;
  tmp << thePDBatom ;  OUT << tmp ;
  return OUT ;
}

string &operator>>(string &input,  PDBatom &thePDBatom){
  string tmpString; double x , y , z ;
  if(input.length()!=80){
    cout<<"ATOM line:\n"<<input<<"\ndoes not contain 80 characters!\n";
    cout<<"Abort!\n";exit(0);
  }
  tmpString.assign( input , 6 , 5 ) ; 
  string2int( tmpString , thePDBatom.serial ) ;
  (thePDBatom.name).assign( input , 12 , 4 ) ;
  (thePDBatom.altLoc).assign( input , 16 , 1 ) ; 
  (thePDBatom.resName).assign( input , 17 , 3 ) ;
  (thePDBatom.chainID).assign( input , 21 , 1 ) ;
  tmpString.assign( input , 22 , 4 ) ;
  string2int(  tmpString , thePDBatom.resSeq ) ;
  (thePDBatom.iCode).assign( input , 26 , 1 ) ;
  tmpString.assign( input , 30 , 8 ) ;
  string2double( tmpString , x ) ; 
  tmpString.assign( input , 38 , 8 ) ;
  string2double( tmpString , y ) ; 
  tmpString.assign( input , 46 , 8 ) ;
  string2double( tmpString , z ) ;
  PDBvector rprime2( x , y , z ) ; 
  tmpString << rprime2 ; 
  thePDBatom.r = rprime2 ;
  tmpString.assign( input , 54 , 6 ) ;
  string2double( tmpString , thePDBatom.occupancy ) ;
  tmpString.assign( input , 60 , 6 ) ;
  string2double( tmpString , thePDBatom.tempFactor ) ;
  (thePDBatom.segID).assign( input , 72 , 4 ) ;
  (thePDBatom.element).assign( input , 76 , 2 ) ;
  (thePDBatom.charge).assign( input , 78 , 2 ) ;
  thePDBatom.record = true ;

  return input ;
}
/*=====================================================================*/
ifstream &operator>>(ifstream &PDB, PDBatom &thePDBatom){
  string entry ;   getline( PDB, entry ) ; /*cout<<"entry="<<entry<<endl;*/
  if( !is_ATOM(entry) ){ return PDB ;}
  thePDBatom.empty() ;  entry >> thePDBatom ;
  return PDB ;
}
/*=====================================================================*/
PDBatom::PDBatom( ){ ( *this ).empty( ) ; }
/*=====================================================================*/
PDBatom::PDBatom( const PDBvector &coordprime ){
  PDBatom tmp ;   tmp.r = coordprime ;  (*this) = tmp ;
}
/*=====================================================================*/
PDBvector & PDBatom::getCoord( ){ return r ; }
/*=====================================================================*/
bool PDBatom::is_recorded( ){ return record ; }
/*=====================================================================*/
bool PDBatom::is_gen_atom_type( gen_atom_type type ){
  string aa,at;
  aa=this->getResName();  this->getAtomName(at);
  switch(type){
  case _CA_: if(at==" CA ") return true; break;
  case _CB_: if(at==" CB ") return true; break;
  case _HV_: if(at.substr(1,1) != "H") return true;  break;
  case _HN_: if(at==" H  ") return true;  break;
  case _METHYLC_:
    if(aa=="ALA" && at==" CB ") return true;
    if(aa=="VAL"){ if(at==" CG1"||at==" CG2") return true; }
    if(aa=="LEU"){ if(at==" CD1"||at==" CD2") return true; }
    if(aa=="ILE"){ if(at==" CD "||at==" CG2") return true; }
    break;
  default: return false;
  }
  return false;
}
/*=====================================================================*/
bool PDBatom::is_backbone( ){
  string backbonePDBnames( " N  , CA , C  , O  , H  " ) ;
  if( backbonePDBnames.find( name ) != string::npos ){ return true ; }
  return false ;
}

bool PDBatom::matchresSeq( const PDBatom &entryAtom ){
  if( resSeq == entryAtom.resSeq ){ return true ; }
  else { return false ; }
}

void PDBatom::empty( ){
  string nullString("    ") ; PDBvector nullPDBvector(0.0, 0.0, 0.0) ;
  serial = 0 ; 
  name.assign( nullString , 0 , 4 ) ; 
  altLoc.assign( nullString , 0 , 1 ) ;
  resName.assign( nullString , 0 , 3 ) ; 
  chainID.assign( nullString , 0 , 1 ) ; 
  resSeq = 0 ;
  iCode.assign( nullString , 0 , 1 ) ; 
  r = nullPDBvector ;     
  occupancy = 0.00 ; 
  tempFactor = 0.00 ;
  segID.assign( nullString , 0 , 4 ) ; 
  element.assign( nullString , 0 , 2 )  ;
  charge.assign( nullString , 0 , 2 ) ; 
  record = false ;
}

string *PDBatom::printTER( ){
  serial++ ; name.assign( "    " ) ;
  string *tmpString = new string(" " );
  *tmpString <<  *this  ; ( *tmpString ).erase( 26 ) ;
  *tmpString += "                                                      " ;
  ( *tmpString ).replace( 0 , 6 , "TER   " ) ;
  return tmpString ;
}

bool PDBatom::isHeavy( ){
  string heavyAtoms( "C N O S" ), tmpString ; 
  tmpString = name.substr( 1 , 1 ) ;
  if( heavyAtoms.find( tmpString ) != string::npos ){
    return true ; 
  }
  else{ return false ; }
}
/*================================================================*/
bool PDBatom::doContact( const PDBatom &otherAtom , const double &cutOff ){
  double c2 = cutOff * cutOff ;
  if( r.d2( otherAtom.r ) < c2 ){ return true ; }
  else{ return false ; }
}
/*================================================================*/
double PDBatom::d2( PDBatom &otherAtom ){
  return r.d2( otherAtom.r ) ;
}
/*================================================================*/
int PDBatom::getResSeq( ){ return resSeq ; }
/*================================================================*/
string *PDBatom::printName( ){
  string *tmpString = new string ;
  *tmpString = name ;
  return tmpString ;
}
/*================================================================*/
void PDBatom::getAtomName( string &atom_name ){
  atom_name = name ;
}
/*================================================================*/
string *PDBatom::printResName( ){
  string *tmpString = new string ;
  *tmpString = resName ;
  return tmpString ;
}
/*================================================================*/
string PDBatom::getResName( ){ return resName ; }
/*================================================================*/
string *PDBatom::printChainID( ){
  string *tmpString = new string ;
  *tmpString = chainID ;
  return tmpString ;
}
/*================================================================*/
void PDBatom::changeResSeq( const int &theResSeq ){
  resSeq = theResSeq ;
}
/*=====================================================================*/
void PDBatom::printChainID( string &theID ){ theID = chainID ; }
/*=====================================================================*/
void PDBatom::assignCoord( double *coord ){ r.assign( coord ) ; }

void PDBatom::assignCoord(PDBvector &the_coord ){ r = the_coord ; }

bool PDBatom::is( const string &theName ){ return name==theName ; }

double *PDBatom::getCoordToDoubleArray( ){
  double *tmp = new double[3] ;  tmp = r.getComponents( ) ;  return tmp ;
}

void PDBatom::getCoordToDoubleArrayII(double *r2 ){ r.getComponentsII( r2 ) ; }

void PDBatom::changeAtomSerialNumber( const int &theSerial ){serial=theSerial;}

int PDBatom::getSerialNumber( ){ return serial ; }

void PDBatom::assignName( string &the_name ){ name = the_name; }

void PDBatom::assignResName( string &the_resName ){ resName = the_resName; }

double PDBatom::angle_with( PDBatom& a1, PDBatom& a2 ){
  /*We define angle between atoms a1,*this,a2 with example: a1<---*this-->a2
    vector from *this to a1 makes 180deg with vector from *this to a2*/
  PDBvector v1=this->PDBvector_to(a1);
  PDBvector v2=this->PDBvector_to(a2);
  return v1.angle_with(v2);
}
/*============================================================*/
PDBvector PDBatom::PDBvector_to( PDBatom &other_atom ){
  PDBvector tmp= other_atom.getCoord( ) ;
  tmp -= r ;
  return tmp ;
}
/*============================================================*/
void PDBatom::translate( const PDBvector &shift ){
  r += shift;
}
/*============================================================*/
void PDBatom::rotate( const PDBvector &axis, const double &angle ){
  r.rotate( axis, angle) ;
}
/*============================================================*/
void PDBatom::rotate( const PDBvector &new_origin, const PDBvector &axis,
		      const double &angle ){
    r.rotate( new_origin, axis, angle) ;
}
/*============================================================*/
void PDBatom::rotate( double ** rotation_matrix ){
  this->r.rotate( rotation_matrix ) ;
}
/*============================================================*/
PDBatom PDBatom::rotateII( double ** rotation_matrix ){
  PDBatom tmp = *this ;
  tmp.r.rotate( rotation_matrix ) ;
  return tmp;
}
/*============================================================*/
bool PDBatom::is_hydrogen( ){
  if( name.substr( 1, 1 ) == "H" ){ return true ; }
  return false ;
}
/*============================================================*/
string PDBatom::get_resName_name( ){ return resName + " " + name ; }
/*============================================================*/
string PDBatom::get_resName_nameII( ){/*correct for bakcbone atoms*/
  string residue = resName ;
  if( this->is_backbone( ) && (residue != "PRO") ){ 
    residue.assign( "BKB") ; 
  }
  return residue + " " + name ; 
}
/*============================================================*/
PDBvector PDBatom::normal_to_plane( PDBatom &a1, PDBatom &a2 ){
  PDBvector v ;
  v = (this->PDBvector_to( a1 ))^(this->PDBvector_to( a2 ));
  v.normalize( ) ;
  return v ;
}
/*============================================================*/
double PDBatom::radius(PDBatom &left,PDBatom &right){
  PDBvector O,A,B;
  O=this->getCoord( ); A=left.getCoord( ); B=right.getCoord( );
  /*cout<<O<<endl<<A<<endl<<B<<endl;*/
  return O.radius(A,B);
}
/*============================================================*/
/*============================================================*/
/*============================================================*/

//----------- class nodePDBatom -----------
nodePDBatom::nodePDBatom( const PDBatom &atomEntry ){
  theAtom = atomEntry ; next = NULL ;
}
//-----------------------------------------


/*----------- class listPDBatom -----------*/

ostream &operator<<(ostream &output, const listPDBatom &theList ){
  nodePDBatom *currentPtr = theList.firstPtr ;
  string tmpString ;
  while( currentPtr != NULL){
    tmpString << (*currentPtr).theAtom ;
    output << tmpString << endl ;
    currentPtr = (*currentPtr).next ;
  }
  return output ;
}
/*=============================================================*/
ifstream &operator>>( ifstream &input, listPDBatom &the_list ){
  the_list.importList( input ) ;
  return input ;
}
/*=============================================================*/
listPDBatom::listPDBatom( ){ firstPtr = lastPtr = NULL ; }
/*=============================================================*/
listPDBatom::listPDBatom(  const char *input_name ){
  ifstream input(input_name);
  /*string line;  getline(input,line);cout<<"line="<<line<<endl;exit(1);*/
  firstPtr = lastPtr = NULL ;
  this->importList( input ) ;
  input.close( ) ;
}
/*=============================================================*/
/*copy constructor*/
listPDBatom::listPDBatom( const listPDBatom &clone ){ 
  /*cout<<"listPDBatom copy constructor\n";*/
  /*cout<<"clone=\n"<<clone;exit(0);*/  
  firstPtr = lastPtr = NULL ;
  nodePDBatom *ptr = clone.firstPtr ;
  while( ptr ){
    /*cout<<ptr->theAtom<<endl;*/
    this->insertAtBack( ptr->theAtom ) ;
    ptr = ptr->next ;
  }
  /*cout<<"this=\n";cout<<firstPtr->theAtom<<endl;
  ptr=firstPtr;
  while( ptr ){cout<<ptr->theAtom<<endl;ptr=ptr->next;}exit(0);
  cout<<*this;exit(0);
  cout<<"finished listPDBatom copyconstructor\n";*/
}
/*=============================================================*/
listPDBatom::~listPDBatom( ){ ( *this ).empty( ) ; }
/*=============================================================*/
void listPDBatom::empty( ){
  if( !isEmpty( ) ){
    while( firstPtr != NULL ){ ( *this ).rmFromFront( ); }
  }
}
/*=============================================================*/
nodePDBatom *listPDBatom::getNewNode( const PDBatom &entryAtom ){
  nodePDBatom *ptr = new nodePDBatom( entryAtom ) ;
  assert( ptr != 0 ) ;
  return ptr ;
}
/*=============================================================*/
bool listPDBatom::isEmpty( ){ return firstPtr == NULL ; }
/*===========================================================*/
void listPDBatom::insertAtFront( const PDBatom &entryAtom ){
  nodePDBatom *newPtr = getNewNode( entryAtom ) ;
  if( isEmpty() ){ firstPtr = lastPtr = newPtr ;}
  else{ ( *newPtr ).next = firstPtr  ; firstPtr = newPtr ; }
}
/*===========================================================*/
void listPDBatom::rmFromFront( ){
  nodePDBatom *currentPtr = firstPtr ;
  firstPtr = ( *currentPtr ).next ; delete currentPtr ;
}
/*===========================================================*/
void listPDBatom::insertAtBack( const PDBatom &entryAtom ){
  nodePDBatom *newPtr = getNewNode( entryAtom ) ;
  if( isEmpty() ){ firstPtr = lastPtr = newPtr ; }
  else{ lastPtr->next = newPtr ; lastPtr = newPtr ; }
}
/*===========================================================*/
void listPDBatom::insertAtBack( const listPDBatom &list2 ){
  nodePDBatom *Ptr2 = list2.firstPtr ;
  while( Ptr2 ){
    this->insertAtBack( Ptr2->theAtom ) ;
    Ptr2 = Ptr2->next ;
  }
}
/*===========================================================*/
int listPDBatom::length( ){
  int count = 0;
  nodePDBatom *currentPtr = firstPtr ;
  while(currentPtr){ 
    count++ ;
    currentPtr=currentPtr->next ;
  }
  return count ;
}
/*=====================================================*/
int listPDBatom::importList( ifstream &PDB ){
  string line ;
  int position;
  PDBatom theAtom ;
  
  this->empty();

  do{/*go to the first atom*/
    position=PDB.tellg( ) ; getline( PDB, line ) ; 
  }while( !( is_ATOM(line) ) && !( PDB.eof() ) ) ; 
  if( PDB.eof() ){ return 1 ; }
  else{ PDB.seekg( position ) ; } /*go back one line to read the atom*/
  bool next=true; /*cout<<"line="<<line<<endl;exit(1);*/
  do{
    PDB>>theAtom; /*cout<<"line="<<line<<endl;*/
    this->insertAtBack( theAtom ) ;
    position=PDB.tellg( ) ; getline( PDB, line ) ;
    if(PDB.eof()||!is_ATOM(line)){next=false;}
    PDB.seekg(position);
  }while(next);
  return 0 ;
}
/*=========================================================*/
PDBatom *listPDBatom::getAtomAt( const int & index ){
  PDBatom *tmpAtom = new PDBatom ;
  if( index > this->length( ) ){ return tmpAtom ; } 
  nodePDBatom *currentPtr = firstPtr ;  int count = 1 ;
  while( count < index ){
    currentPtr = currentPtr->next ;  count++ ;
  }
  *tmpAtom = currentPtr->theAtom ;
  return tmpAtom ;
}
/*=========================================================*/
PDBatom listPDBatom::getAtomAtII( const int & index ){
  PDBatom *tmp = this->getAtomAt( index ) ;
  return *tmp ;
}
/*=========================================================*/
string *listPDBatom::printAtomAt( const int &index ){
  string *tmpPtr = new string ;
  *tmpPtr << *( ( *this ).getAtomAt( index ) );
  return tmpPtr ;
}
/*=========================================================*/
bool listPDBatom::matchResSeq( const PDBatom &entryAtom ){
  if( isEmpty( ) ){ return true ; }
  if( ( ( *( *this ).firstPtr ).theAtom ). matchresSeq( entryAtom ) ){
    return true ;
  }
  else{  return false ; }
}
/*==========================================================*/
listPDBatom & listPDBatom::operator=( listPDBatom &right ){
  /*cout<<"listPDBatom assignment operator\n";*/
  if( ! this->isEmpty() ){  this->empty( ) ; }
  if( ( &right != this )  && !( right.isEmpty() ) ) {
    firstPtr = getNewNode( right.firstPtr->theAtom ) ;
    nodePDBatom *currentRightPtr=right.firstPtr->next, *currentPtr=firstPtr ;
    while( currentRightPtr != NULL ){
      currentPtr->next = getNewNode( currentRightPtr->theAtom ) ;
      currentPtr = currentPtr->next ;
      currentRightPtr = currentRightPtr->next ;
    }
    lastPtr = currentPtr ;
  }
  return *this ;
}
/*==========================================================*/
listPDBatom & listPDBatom::operator=( const listPDBatom &right ){
  /*cout<<"listPDBatom & listPDBatom::operator=(const listPDBatom &)\n";*/
  /*cout<<right; exit(0);*/
  listPDBatom tmp( right ) ;  /*cout<<tmp;exit(0);*/
  *this = tmp ; /*previous operator assignment definition*/
  /*cout<<(*this);exit(0);*/
  return *this ;
}
/*==========================================================*/
PDBatom *listPDBatom::heavyGeomCent( ){
  PDBatom *tmpAtom = new PDBatom ;
  string *tmpString ;
  PDBvector gm( 0.0 , 0.0 , 0.0 );
  //  string tmp2 ; tmp2 << gm ; cout << tmp2 << endl ;
  int count = 0 ;
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    if( ( ( *currentPtr ).theAtom ).isHeavy( ) ){
      // string tmp2 ; tmp2 << ( *currentPtr ).theAtom ; cout << tmp2 << endl ; 
      gm += ( ( *currentPtr ).theAtom ).r ;  count++ ;
    }
    currentPtr = ( *currentPtr ).next ;
  }
  gm /= (double)count ;
  *tmpAtom = ( *firstPtr ).theAtom ;
  ( *tmpAtom ).r = gm ;
  return tmpAtom ;
}

bool listPDBatom::doContactHeavyAtom( PDBatom &otherAtom , double &cutOff ){
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    if( ( ( *currentPtr ).theAtom ).isHeavy( ) ){
      if( ( (*currentPtr).theAtom ).doContact( otherAtom , cutOff ) ){
	return true ; 
      }
    }
    currentPtr = ( *currentPtr ).next ;
  }
  return false ;
}

bool listPDBatom::doContactHeavyAtom( PDBatom &otherAtom , 
				      double &cutOff , int &BeginIndex ){
  nodePDBatom *currentPtr = firstPtr ;
  int count = 1 ;
  while( (currentPtr != NULL) && (count < BeginIndex) ){
    currentPtr = ( *currentPtr ).next ; count ++ ;
  }
  while( currentPtr != NULL ){
    if( ( ( *currentPtr ).theAtom ).isHeavy( ) ){
      if( ( (*currentPtr).theAtom ).doContact( otherAtom , cutOff ) ){
	return true ; 
      }
    }
    currentPtr = ( *currentPtr ).next ;
  }
  return false ;
}
/*========================================================================*/
bool listPDBatom::doContactHeavyAtom( listPDBatom &otherList , 
				      double &cutOff ){
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
   if( ( ( *currentPtr ).theAtom ).isHeavy( ) ){
     if( otherList.doContactHeavyAtom( (*currentPtr).theAtom, cutOff) ){
       return true ;
     }
   }
   currentPtr = ( *currentPtr ).next ;
  }
  return false ;
}
/*=========================================================================*/
bool listPDBatom::doContactAnyone(listPDBatom &otherList,double &cutOff){
  nodePDBatom *pt=firstPtr,*pt2;
  PDBatom at;
  while(pt){
    at=pt->theAtom;
    pt2=otherList.firstPtr;    
    while(pt2){
      if( at.doContact(pt2->theAtom,cutOff) ) return true;
      pt2=pt2->next;
    }
    pt=pt->next;
  }
  return false;
}
/*=========================================================================*/
int listPDBatom::numberAnyContacts(listPDBatom &otherList,double &cutOff){
  nodePDBatom *pt=firstPtr,*pt2;
  PDBatom at;
  int n=0;
  while(pt){
    at=pt->theAtom;
    pt2=otherList.firstPtr;    
    while(pt2){
      if( at.doContact(pt2->theAtom,cutOff) ) n++;
      pt2=pt2->next;
    }
    pt=pt->next;
  }
  return n;
}
/*=========================================================================*/
void listPDBatom::changeResSeq( const int &theResSeq ){
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr != NULL ){
    ( ( *currentPtr ).theAtom ).changeResSeq( theResSeq ) ;
    currentPtr = ( *currentPtr ).next ;
  }
}
/*=========================================================================*/
void listPDBatom::renumberFully( const int &beginIndex ){
  nodePDBatom *currentPtr = firstPtr ;
  int index = beginIndex ;
  while( currentPtr != NULL ){
    ( ( *currentPtr ).theAtom ).changeResSeq( index ) ;
    index++ ;
    currentPtr = ( *currentPtr ).next ;
  }
}

int listPDBatom::renumberFullyAtomSerialNumber( const int &firstIndex ){
  int index = firstIndex ;
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr ){
    ( ( *currentPtr ).theAtom ).changeAtomSerialNumber( index ) ;
    index++ ; currentPtr = ( *currentPtr ).next ;
  }
  return index ;
}

bool listPDBatom::isThereAtom( const string &theName ){
  nodePDBatom *currentPtr = firstPtr ;
  string *currentName ;
  while( currentPtr != NULL ){
    currentName = ( ( *currentPtr ).theAtom ).printName( ) ;
    if( theName == *currentName){ 
      return true ; 
    }
    delete currentName ;
    currentPtr = ( *currentPtr ).next ;
  }
  return false ;
}
/*==============================================================*/
PDBatom *listPDBatom::getAtomFromName( const string &theName ){
  nodePDBatom *currentPtr = firstPtr ;
  PDBatom *tmp = NULL ;
  string *currentName ;
  while( currentPtr != NULL ){
    currentName = currentPtr->theAtom.printName( ) ;
    if( theName == *currentName){
      tmp = new PDBatom ;
      *tmp =  currentPtr->theAtom ;
      return tmp ;
    }
    delete currentName ;
    currentPtr = currentPtr->next ;
  }
  return tmp ;
}
/*==============================================================*/
PDBatom *listPDBatom::pointToAtomFromName( const string &theName ){
  nodePDBatom *currentPtr = firstPtr ;
  string currentName ;
  while( currentPtr != NULL ){
    currentPtr->theAtom.getAtomName(currentName) ;
    if( theName == currentName){ return &(currentPtr->theAtom) ; }
    currentPtr = currentPtr->next ;
  }
  return NULL ;
}
/*==============================================================*/
void listPDBatom::printChainID( string &theID ){
  if( !this->isEmpty( ) ){
    firstPtr->theAtom.printChainID( theID ) ;
  }
}
/*==============================================================*/
string *listPDBatom::printTER( void ){
  string *tmpString = ( (*lastPtr).theAtom ).printTER() ;
  return tmpString ;
}
/*==================================================================*/
void listPDBatom::changeCoordAt( const int &index , double *coord){
  if( index > this->length( ) ){
    cout<<"\nERROR: the listPDBatom has less than "<<index<<" entries\n" ;
    exit(1) ;
  }
  nodePDBatom *currentPtr = firstPtr ;
  for( int i=1 ; i< index ; i++){ currentPtr = currentPtr->next; }
  currentPtr->theAtom.assignCoord( coord ) ;
}
/*==================================================================*/
void listPDBatom::changeCoordAt( const int &index , PDBvector &theVector){
  double *coord = new double[3] ;
  for( int i=0; i<3; i++ ) { coord[ i ] = theVector[ i ] ; }
  return changeCoordAt( index, coord ) ;
}
/*==================================================================*/
bool listPDBatom::changeCoordOfAtomName( const string &theName,
					 double *coord ){
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr ){
    if( ( (*currentPtr).theAtom ).is( theName ) ){
      ( (*currentPtr).theAtom ).assignCoord( coord ) ;
      return true ;
    }
    currentPtr = (*currentPtr).next ;
  }
  return false ;
}
/*================================================================*/
double *listPDBatom::getCoordOfAtomName( const string &theName ){
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr ){
   if( currentPtr->theAtom.is( theName ) ){
     return currentPtr->theAtom.getCoordToDoubleArray( ) ;
   }
   currentPtr = currentPtr->next ;
  }
  return NULL ;
}
/*================================================================*/
int listPDBatom::getSerialNumberOfAtomName( const string &AtomName ){
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr ){
    if( ( (*currentPtr).theAtom ).is( AtomName ) ){
      return ( (*currentPtr).theAtom ).getSerialNumber( ) ;
    }
    currentPtr = (*currentPtr).next ;
  }
  return 0 ;
}
/*================================================================*/
char * listPDBatom::printContiguousDist( void ){
  nodePDBatom *currentPtr = firstPtr, *nextPtr ;
  if( currentPtr ){ nextPtr = (*currentPtr).next ; }
  else{ return NULL ; }

  char *list = (char*) malloc( sizeof( char ) * (*this).length( ) * 20 ) ; 
  char line[20] ;  double dist ;  int k = 1 ;

  while( nextPtr ){
    dist = sqrt( ((*currentPtr).theAtom).d2( (*nextPtr).theAtom ) ) ;
    sprintf( line, "%4d %4d %8.3f\n", k , k+1, dist ) ; k++ ; 
    strcat( list, line ) ;
    currentPtr = nextPtr ; nextPtr = (*nextPtr).next ;
  }

  return list ;
}

char * listPDBatom::printNext_N_NeighDist( int &N ){
  nodePDBatom *currentPtr = firstPtr, *nextPtr = currentPtr ;

  if( (*this).length() <= N ){ return NULL ; }
  for(int i=0; i<N; i++){ nextPtr = (*nextPtr).next ; }

  char *list = (char*) malloc( sizeof( char ) * (*this).length( ) * 20 ) ; 
  char line[20] ;  double dist ;  int k = 1 ;

  while( nextPtr ){
    dist = sqrt( ((*currentPtr).theAtom).d2( (*nextPtr).theAtom ) ) ;
    sprintf( line, "%4d %4d %8.3f\n", k , k+N, dist ) ; k++ ; 
    strcat( list, line ) ;
    currentPtr = (*currentPtr).next ; nextPtr = (*nextPtr).next ;
  }

  return list ;
}
/*=============================================================*/
int listPDBatom::createContactMap( ostream &map, double &cutOff ){
  nodePDBatom *atom1Ptr = firstPtr ;
  nodePDBatom *atom2Ptr = ( *(*atom1Ptr).next ).next ;//avoid next neighbors 
  double c2=cutOff*cutOff, d2;
  int l = (*this).length( ) ;

  for(int i=1; i<l-1; i++){
    for(int j=i+2; j<=l; j++){ //avoid next neighbors
      d2=( (*atom1Ptr).theAtom ).d2( (*atom2Ptr).theAtom ) ;
      if( d2<c2 ){ map << i <<" "<< j <<endl ;  }
      atom2Ptr = (*atom2Ptr).next ;
    }
    atom1Ptr = (*atom1Ptr).next ;
    atom2Ptr = ( *(*atom1Ptr).next ).next ;
  }
  return 0 ;
}
/*=============================================================*/
int listPDBatom::createContactMap( double **map, double &cutOff ){
  /*cout<<"listPDBatom::createContactMap\n";*/
  nodePDBatom *atom1Ptr = firstPtr ;
  nodePDBatom *atom2Ptr = ( *(*atom1Ptr).next ).next ;//avoid next neighbors 
  double c2=cutOff*cutOff, d2;
  int l = (*this).length( ) ;

  /*avoid next neighbors, and first index bigger than second index*/
  for(int i=1; i<l-1; i++){
    for(int j=i+2; j<=l; j++){
      d2=( (*atom1Ptr).theAtom ).d2( (*atom2Ptr).theAtom ) ;
      if( d2<c2 ){ 
	/*printf("before map[%d][%d]=%3.1f\n",j-1,i-1,map[j-1][i-1]);*/
	map[j-1][i-1]++ ; /*note first index bigger than second index*/ 
	/*printf("after map[%d][%d]=%3.1f\n",j-1,i-1,map[j-1][i-1]);*/
      }
      atom2Ptr = (*atom2Ptr).next ;
    }
    atom1Ptr = (*atom1Ptr).next ;
    atom2Ptr = ( *(*atom1Ptr).next ).next ;
  }
  return 0 ;
}
/*=============================================================*/
int listPDBatom::createContactMap( listPDBatom &list2, double **map, 
				   double &cutOff ){
  /*note that a contact can have a value of 0, 1, or 2. For example, atom 7 of
    list contacts atom 3 of list2, but also atom 7 of list2 contacts atom 3 of
    list*/
  nodePDBatom *atom1Ptr = firstPtr ;
  nodePDBatom *atom2Ptr = NULL ;
  double c2=cutOff*cutOff, d2;
  int l = this->length( ) ;
  int l2 = list2.length( ) ;

  for(int i=1; i<=l; i++){
    atom2Ptr = list2.firstPtr ;
    for(int j=1; j<=l2; j++){
      d2=( atom1Ptr->theAtom ).d2( atom2Ptr->theAtom ) ;
      if( d2<c2 ){ /*first index bigger or equal than second index*/
	/*printf("%d %d %f\n",i,j,sqrt(d2));*/
	if(i>j){ 
	  /*printf("before map[%d][%d]=%3.1f\n",i-1,j-1,map[i-1][j-1]);*/
	  map[i-1][j-1]++ ; 
	  /*printf("after map[%d][%d]=%3.1f\n",i-1,j-1,map[i-1][j-1]);*/
	}
	else{ map[j-1][i-1]++ ; }
      }
      atom2Ptr = atom2Ptr->next ;
    }
    atom1Ptr = atom1Ptr->next ;
  }
  return 0 ;
}
/*=============================================================*/
double listPDBatom::drms( listPDBatom &list2 ){
  int l = list2.length() ;
  if( l != (*this).length() ){ 
    cout <<"chains have different lenght\n"; exit(0);
  }
  else if( l<3 ){cout<<"chain too short!\n"; exit(0); }
  //We obviate next neighbors in computing the drms.
  nodePDBatom *chain1Ptr1, *chain2Ptr1, *chain1Ptr2, *chain2Ptr2 ;
  
  chain1Ptr1 = firstPtr ;
  chain2Ptr1 = list2.firstPtr ;
  chain1Ptr2 = ( *(*chain1Ptr1).next ).next ;
  chain2Ptr2 = ( *(*chain2Ptr1).next ).next ;

  double d=0.0, d1, d2 ;
  for( int i=1; i<l-1; i++){
    for(int j=i+2; j<=l; j++){
      d1 = ( (*chain1Ptr1).theAtom ).d2( (*chain1Ptr2).theAtom );
      d1 = sqrt( d1 );
      d2 = ( (*chain2Ptr1).theAtom ).d2( (*chain2Ptr2).theAtom );
      d2 = sqrt( d2 );
      d += (d1-d2)*(d1-d2) ;
      chain1Ptr2 = (*chain1Ptr2).next ;
      chain2Ptr2 = (*chain2Ptr2).next ;
    }
    chain1Ptr1 = (*chain1Ptr1).next ;
    chain2Ptr1 = (*chain2Ptr1).next ;
    chain1Ptr2 = ( *(*chain1Ptr1).next ).next ;
    chain2Ptr2 = ( *(*chain2Ptr1).next ).next ;
  }
  //there are l(l-1)/2 - (l-1) pairs, excluding next neighbors.  
  return sqrt(2*d/(l*l));
}

double listPDBatom::rmsd( listPDBatom &list2 ){
  int l = list2.length( ) ;
  if( l != (*this).length() ){ return -1 ; }
  PDBvector cm1, cm2 ; cm1 = (*this).get_CM( ) ;  cm2 = list2.get_CM( ) ;
  cm2 -= cm1 ;
  listPDBatom list3; list3 = *this; list3.translate( cm2 ) ;
  double **r2, **r3 ;
  r2 = list2.dump_coordinates_to_array( ) ;
  r3 = list3.dump_coordinates_to_array( ) ;
  return get_rms( r3, r2, l ) ;
}
/*===================================================================*/
/*transform list2 so that adapts to *this by minimizing rmsd.*/
double listPDBatom::adapt_to_this_chain_by_rmsd( listPDBatom &list2 ){
  int l = list2.length( ) ;
  if( l != (*this).length() ){ return -1 ; }
  PDBvector cm1, cm2 ;  cm1 = (*this).get_CM( ) ; cm2 = list2.get_CM( ) ;
  cm1 -= cm2 ;
  listPDBatom list3; list3 = list2; list3.translate( cm1 ) ;
  double **r1, **r3, **r4 ;
  r1 = (*this).dump_coordinates_to_array( ) ;
  r3 = list3.dump_coordinates_to_array( ) ; 
  r4 = alloc_array( l, 3 ) ;
//transform r3 coordinates so that rmsd between r1 and r3 is minimized, and 
//store transformed coordinates in r4.
  double rmsd = get_rmsII( r3, r1, r4, l ) ; 
  list2.dump_coordinates_to_chain( r4 ) ;
  return rmsd ;
}
/*===================================================================*/
double listPDBatom::output_tranformed_chain_by_rmsd(listPDBatom &list2, 
						    listPDBatom &list3 ){
  list3 = list2 ;
  return (*this).adapt_to_this_chain_by_rmsd( list3 ) ;
}
/*===================================================================
  rot_matrix is the rotation needed on *this to adapt to list2 by 
  minimizing the rmsd*/
double listPDBatom::rot_matrix_to_adapt_by_rmsd_to( listPDBatom &list2,
						    double **rot_matrix ){
  int l = list2.length( ) ;
  if(  (this->isEmpty( ) || list2.isEmpty( ) ) || ( this->length( ) != l )  ){
    return DBL ;
  }
  PDBvector cm1, cm2 ;  cm1 = this->get_CM( ) ; cm2 = list2.get_CM( ) ;
  cm1 -= cm2 ;
  listPDBatom tmp;
  tmp = *this; 
  tmp.translate( cm1 ) ; 
  double **r1, **r2 ;
  r1 = tmp.dump_coordinates_to_array( ) ;
  r2 = list2.dump_coordinates_to_array( ) ;
  return get_rot_matrix_by_rms( r1, r2, l, rot_matrix ) ;/*returns the rmsd*/
}
/*===================================================================*/
int listPDBatom::diffContactMap( ifstream &listOfContacts, listPDBatom &list2,
				 double &cutOff){
  /*Reads file listOfContacts, each row is a contact, then finds the differences between list2 and (*this). For each contact, if the contact is present in one chain, and missing in the other, then we assign one, and zero otherwise. Thus the total number of differences can be reported.
   */
  if( list2.length() != (*this).length() ){ 
    cout<<"EXIT MESSAGE: Chains of different lenghts when listPDBatom::diffContactMap called !"<<list2.length()<<" "<<(*this).length()<<" \n";
    exit(0);
  }
  PDBatom *atom1, *atom2 ;
  int r1, r2, d, diff=0 ;
  double c2 = cutOff*cutOff ;

  //record actual position of file listOfContacts
  int position=listOfContacts.tellg( ) ;
  while( !listOfContacts.eof() ){
    listOfContacts >> r1 >> r2 ;
    atom1 = (*this).getAtomAt( r1 ) ; atom2 = (*this).getAtomAt( r2 ) ; 
    if( (*atom1).d2(*atom2) < c2 ){ d = 1 ;}
    else{ d = 0 ; }
    delete atom1 ; delete atom2 ;

    atom1 = list2.getAtomAt( r1 ) ; atom2 = list2.getAtomAt( r2 ) ;
    if( (*atom1).d2(*atom2) < c2 ){ d -= 1 ;}
    delete atom1 ; delete atom2 ;

    diff += abs( d ) ;
  }
  //rewind to original position, so that file can be used more than once.
  if( listOfContacts.eof() ){ listOfContacts.clear( ios::goodbit ) ; }
  listOfContacts.seekg( 0 ) ;
  
  return diff ;
}


double listPDBatom::NDrms( ifstream &listOfContacts, listPDBatom &list2 ){
  /*Reads file listOfContacts, each row is a native contact, then finds the differences between list2 and (*this). For each native contact, it calculates the distance in both chains, and calculates the drms. */
  int nn = 0 ;
  double ndRms=0.0, d1, d2, d ;

  if( list2.length() != (*this).length() ){ 
    cout<<"EXIT MESSAGE: Chains of different lenghts when listPDBatom::NDrms called !"<<list2.length()<<" "<<(*this).length()<<" \n";
    exit(0);
  }

  PDBatom *atom1, *atom2 ;
  int r1, r2 ;

  //record actual position of file listOfContacts
  int position=listOfContacts.tellg( ) ;

  while( !listOfContacts.eof() ){
    listOfContacts >> r1 >> r2 ; nn++;

    atom1 = (*this).getAtomAt( r1 ) ; atom2 = (*this).getAtomAt( r2 ) ; 
    d1 = (*atom1).d2( *atom2 ) ;
    delete atom1 ; delete atom2 ;

   atom1 = list2.getAtomAt( r1 ) ; atom2 = list2.getAtomAt( r2 ) ;
   d2 = (*atom1).d2( *atom2 ) ;
   delete atom1 ; delete atom2 ;

   d = sqrt(d1) - sqrt(d2) ; ndRms += d*d ;
  }

  //rewind to original position, so that file can be used more than once.
  if( listOfContacts.eof() ){ listOfContacts.clear( ios::goodbit ) ; }
  listOfContacts.seekg( position ) ;

  return sqrt(ndRms/nn) ;
}

int listPDBatom::numContacts( ifstream &listOfContacts, double &cutOff ){
  PDBatom *atom1, *atom2 ;
  int n=0, r1, r2 ;
  double c2=cutOff*cutOff ;

  //record actual position of file listOfContacts
  int position=listOfContacts.tellg( ) ;
  
  while( !listOfContacts.eof() ){
    listOfContacts >> r1 >> r2 ;
    atom1 = (*this).getAtomAt( r1 ) ; atom2 = (*this).getAtomAt( r2 ) ; 
    if( (*atom1).d2(*atom2) < c2 ){ n++ ; }
  }
  //rewind to original position, so that file can be used more than once.
  if( listOfContacts.eof() ){ listOfContacts.clear( ios::goodbit ) ; }
  listOfContacts.seekg( position ) ;  
  return n ;
}

PDBvector listPDBatom::get_CM( ){
  PDBvector tmp ;
  if( !this ){ cout<<"EXIT MESSAGE: list is empty" ; return tmp ; }
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr != NULL){ 
    tmp += currentPtr->theAtom.r ; 
    currentPtr = currentPtr->next ;
  }
  tmp /= (*this).length( ) ;
  //  cout << tmp << endl ;
  return tmp ;
}
/*============================================================*/
void listPDBatom::translate( const PDBvector &the_vector ){
  if( !this ){ cout<<"EXIT MESSAGE: list is empty" ; return ; }
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr != NULL){ 
    currentPtr->theAtom.r +=  the_vector ;
    currentPtr = currentPtr->next ;
  }
  return ;
}
/*============================================================*/
void listPDBatom::rotate( const PDBvector &axis, const double &angle ){
  if( !this ){ cout<<"EXIT MESSAGE: list is empty" ; return ; }
   nodePDBatom *currentPtr = firstPtr ;
   while( currentPtr != NULL){ 
    currentPtr->theAtom.r.rotate( axis, angle ) ;
    currentPtr = currentPtr->next ;
  }
  return ; 
}
/*============================================================*/
void listPDBatom::rotate( const PDBvector &new_origin, const PDBvector &axis, 
	     const double &angle ){
  if( !this ){ return ; }
   nodePDBatom *currentPtr = firstPtr ;
   while( currentPtr != NULL){ 
    currentPtr->theAtom.r.rotate( new_origin, axis, angle ) ;
    currentPtr = currentPtr->next ;
  }
  return ; 
}
/*============================================================*/
void listPDBatom::rotate( double ** rotation_matrix ){
  if( !this ){ cout<<"EXIT MESSAGE: list is empty" ; return ; }
   nodePDBatom *currentPtr = firstPtr ;
   while( currentPtr != NULL){ 
    currentPtr->theAtom.r.rotate( rotation_matrix ) ;
    currentPtr = currentPtr->next ;
  }
  return ; 
}
/*============================================================*/

//==obtain the radius of the sphere with center in the center of mass that
//==encloses the chain.
double listPDBatom::get_surrounding_sphere( ){
  double radius=0.0, dd ;
  listPDBatom tmp ; 
  tmp = *this ;
  PDBvector cm = tmp.get_CM( ) ; cm *= -1 ;
  tmp.translate( cm ) ; 
  nodePDBatom *currentPtr = tmp.firstPtr ;
  while( currentPtr != NULL){ 
    dd = ( currentPtr->theAtom.r).norm2( ) ;
    if( dd > radius ){ radius = dd ; }
    //    cout << currentPtr->theAtom << endl ;
    currentPtr = currentPtr->next ;
  }
  return sqrt( radius ) ; 
}
/*================================================================*/
double listPDBatom::get_distance_between( int&i1, int&i2 ){
  PDBatom *atom1, *atom2 ;
  atom1 = getAtomAt( i1 ) ; atom2 = getAtomAt( i2 ) ;
  return sqrt( (*atom1).d2( *atom2 ) );
}
/*================================================================*/
void listPDBatom::get_coord_at_atom_number( int&position, double*r ){
  if( position >= (*this).length() ){ return ; }
  PDBatom *atom ; atom = getAtomAt( position ) ;  
  (*atom).getCoordToDoubleArrayII( r ) ;
}
/*================================================================*/
bool listPDBatom::getCoordOfAtomWithAtomIndex( const int &index, double *coords ){
  if( this->isEmpty( ) ){ return false ; }
  nodePDBatom *currentPtr = firstPtr ;
  while(currentPtr){
    if( currentPtr->theAtom.getSerialNumber( ) == index ){
       currentPtr->theAtom.getCoordToDoubleArrayII( coords ) ;
      return true ; 
    }
    currentPtr = currentPtr->next ;
  }
  return false ;
}
/*================================================================*/
bool listPDBatom::getNameOfAtomWithAtomIndex( const int &index,
					      string &atom_name ){
  if( this->isEmpty( ) ){ return false ; }
  nodePDBatom *currentPtr = firstPtr ;
  while(currentPtr){
    if( currentPtr->theAtom.getSerialNumber( ) == index ){
      currentPtr->theAtom.getAtomName( atom_name ) ;
      return true ; 
    }
    currentPtr = currentPtr->next ;
  }
  return false ;
}
/*================================================================*/
double **listPDBatom::dump_coordinates_to_array( ){
  int l = (*this).length( ) ;
  PDBatom *atom ; 
  double **r = new double* [l] ;
  r[0] = new double [3*l] ;
  for( int i=1; i<l; i++ ){ r[i] = r[i-1] + 3 ; }
  for( int i=0; i<l; i++ ){ 
    atom = getAtomAt( i+1 ) ;  
    (*atom).getCoordToDoubleArrayII( r[i] ) ;
    delete atom ;
  }
  return r ;
}
/*============================================================*/
PDBatom* listPDBatom::returnAdressOfAtomWithAtomIndex( const int &index ){
  PDBatom* tmp=NULL ;
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr ){
    if( currentPtr->theAtom.getSerialNumber( ) == index ){
      return &( currentPtr->theAtom ) ;
    }
    currentPtr = currentPtr->next ;
  }
  return tmp ;
}
/*============================================================*/
void listPDBatom::dump_coordinates_to_chain( double **r ){
  int l=(*this).length( ) ;
  for(int i=1; i<=l; i++ ){ changeCoordAt( i, r[i-1] ) ; }
}
/*============================================================*/
void listPDBatom::filterOnly( const string *names, const int &n ){
  if( this->isEmpty( ) ){ return ; }
  bool flag=false;
  listPDBatom filtered ;/*filtered is empty right now*/
  string atom_name ;
  nodePDBatom *currentPtr = firstPtr ;
  while(currentPtr){
    for( int i=0; i<n; i++ ){
      atom_name = names[ i ] ;
      if(atom_name.length()!=4){ 
	cout<<"ERROR in listPDBatom::filterOnly.";
	cout<<"Atom name exceeded 4 characters\n";
	exit(0);
      }
      if( currentPtr->theAtom.is( names[ i ] ) ){ 
	flag = true ;
	break ;
      }
    }
    if( flag ){ filtered.insertAtBack( currentPtr->theAtom ) ; }
    currentPtr = currentPtr->next ;
    flag = false ;
  }
  this->empty( ) ;
  *this = filtered ;/*if all atoms were left out, then *this is an empty list*/
}
/*============================================================*/
double listPDBatom::minD( listPDBatom &neighborList ){
  if( this->isEmpty( ) || neighborList.isEmpty( ) ){
    return DBL; /*Huge dist!*/
    cout<<"One or both chains are empty!"<<endl;
  }
  double mind2 = DBL, d2 ;
  
  nodePDBatom *l1Ptr = firstPtr, *l2Ptr ;
  while(l1Ptr){
    l2Ptr = neighborList.firstPtr ;
    while(l2Ptr){
      d2 = l1Ptr->theAtom.d2( l2Ptr->theAtom ) ;
      if( mind2 > d2 ){ mind2 = d2 ; }
      l2Ptr = l2Ptr->next ;
    }
    l1Ptr = l1Ptr->next ;
  }
  return sqrt( mind2 ) ;
}
/*============================================================*/
PDBatom listPDBatom::getAtomWithAtomIndex( const int &theIndex ){
  PDBatom tmp ;
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr ){
    if( currentPtr->theAtom.getSerialNumber( ) == theIndex ){
      return currentPtr->theAtom ;
    }
    currentPtr = currentPtr->next ;
  }
  return tmp ;
}
/*============================================================*/
bool listPDBatom::hasAtomWithAtomIndex( const int &theIndex ){
  nodePDBatom *currentPtr = firstPtr ;
  while( currentPtr ){
    if( currentPtr->theAtom.getSerialNumber( ) == theIndex ){
      return true ;
    }
    currentPtr = currentPtr->next ;
  }
  return false ;    
}
/*============================================================*/
void listPDBatom::remove_all_hydrogens( ){
  nodePDBatom *currentPtr = firstPtr ;
  nodePDBatom *prevPtr = NULL ;
  nodePDBatom *removed ;
  while( currentPtr ){
    if( currentPtr->theAtom.is_hydrogen( ) ){
      removed = currentPtr ;
      if( prevPtr ){
	prevPtr->next = currentPtr->next ;
	currentPtr = currentPtr->next ;
	delete removed ;
      }
      else{
	firstPtr = firstPtr->next ;
	currentPtr = firstPtr ;
	delete removed ;
      }
    }
    else{
      prevPtr = currentPtr ;
      currentPtr = currentPtr->next ;
    }
  }
}
/*============================================================*/
int listPDBatom::extract(gen_atom_type type, listPDBatom &list){
  nodePDBatom *pt=this->firstPtr;
  int n=0;
  while(pt){
    if(pt->theAtom.is_gen_atom_type(type)){ 
      list.insertAtBack(pt->theAtom); 
      n++;
    }
    pt=pt->next;
  }
  return n;
}
/*============================================================*/
/*============================================================*/
/*============================================================*/
/*============================================================*/


//----------- class nodelistPDBatom -----------
nodelistPDBatom::nodelistPDBatom( const listPDBatom &listEntry,
				  const string &theListID      ){
  /*cout<<"nodelistPDBatom(..)\n";cout<<listEntry; exit(0);*/
  theList=listEntry ; /*cout<<theList;exit(0);*/
  listID=theListID;  /*cout<<"listID=\""<<listID<<"\"\n";exit(0);*/
  next=NULL;
}
//-----------------------------------------


//------------ class listsPDBatom -----------*/

ostream &operator<<(ostream &output, const listsPDBatom &theLists ){
  nodelistPDBatom *currentPtr = theLists.firstPtr ;
  string tmpString ;
  while( currentPtr != NULL){
    output<<currentPtr->theList;
    currentPtr = currentPtr->next ;
    if(currentPtr)
      if(currentPtr->next) output<<"TER   \n";
      else output<<"END   \n";
  }
  return output ;
}
/*=============================================================*/
ifstream &operator>>( ifstream &input, listsPDBatom &the_lists ){
  cout<<"Not yet implemented!\n";
  /*the_lists.importLists( input ) ; */
  return input ;
}
/*=============================================================*/
nodelistPDBatom *listsPDBatom::getNewNode( const listPDBatom &entryList,
					   const string &theListID      ){
  /*cout<<"listsPDBatom::getNewNode(..)\n";cout<<entryList; exit(0);*/
  nodelistPDBatom *ptr = new nodelistPDBatom(entryList,theListID);
  assert(ptr!=0); /*cout<<ptr->theList;*/
  return ptr ;
}
/*=============================================================*/
listsPDBatom::listsPDBatom( ){ firstPtr = lastPtr = NULL ; }
/*=============================================================*/
bool listsPDBatom::isEmpty( ){ return firstPtr == NULL ; }
/*===========================================================*/
void listsPDBatom::rmFromFront( ){
  nodelistPDBatom *currentPtr = firstPtr ;
  firstPtr = currentPtr->next ;
  delete currentPtr;/*remove listPDBatom and listID associated to the node*/
}
/*===========================================================*/
void listsPDBatom::empty( ){
  if( !isEmpty( ) )
    while( firstPtr != NULL ){ this->rmFromFront( ); }
}
/*=============================================================*/
listsPDBatom::~listsPDBatom( ){ this->empty(); }
/*=============================================================*/
void listsPDBatom::insertAtBack( const listPDBatom &entryList,
				 const string &theListID       ){
  /*cout<<"listsPDBatom::insertAtBack(..)\n"; cout<<entryList; exit(0);*/
  nodelistPDBatom *newPtr = getNewNode( entryList, theListID ) ;
  if( isEmpty() ){ firstPtr = lastPtr = newPtr ; }
  else{ lastPtr->next = newPtr ; lastPtr = newPtr ; }
}
/*===========================================================*/
int listsPDBatom::length( ){
  int count = 0;
  nodelistPDBatom *currentPtr = firstPtr ;
  while(currentPtr){ 
    count++ ;
    currentPtr=currentPtr->next ;
  }
  return count ;
}
/*=====================================================*/
listPDBatom &listsPDBatom::operator[]( const int &index ){
  if(index<1 || index>this->length() ){ 
    cout<<"\nERROR: index out of range\n";
    exit(1) ;
  }
  nodelistPDBatom *currentPtr = firstPtr ;
  int n=index-1;
  while(n){ currentPtr=currentPtr->next; n--; } 
  return currentPtr->theList;
}
/*============================================================*/



//------------ class PDBamino_acid  --------------

ostream &operator<<(ostream &output , const PDBamino_acid &theAA ){
  output << theAA.backbone << theAA.residue ;
  return output ;
}
/*=================================================================*/
ifstream &operator>>(ifstream &PDB, PDBamino_acid &theAminoAcid){
  int position , resSeq ;
  string line ; //string junkLine ;
  PDBatom atomEntry ;

  theAminoAcid.empty( ) ;

  position = PDB.tellg( ) ; getline( PDB, line );

  while( !(is_ATOM(line)) || !(is_AMINOACID(line)) ){/*find first ATOM entry*/
    //cout <<"EXIT MESSAGE: We are reading a non-ATOM line\n" ;
    if( PDB.eof() ) return PDB ; 
    if( is_TER(line) || is_END(line) ){ PDB.seekg( position ); return PDB ; }
    position = PDB.tellg( ) ; getline( PDB, line );
  }

  line >> atomEntry ; resSeq = atomEntry.getResSeq( ) ; //initialize

  while( atomEntry.getResSeq( ) == resSeq || is_TER(line) ){
    if( is_ATOM(line) ){ theAminoAcid.addAtom( atomEntry ) ; }
    position = PDB.tellg( ) ; getline( PDB, line ) ;
    if( is_ANISOU(line) || is_SIGATM(line) || is_HETATM(line) ){ continue ; }
    if( PDB.eof() ){ return PDB ; }
    if( !(is_ATOM(line)) ){ PDB.seekg( position ); return PDB ; }
    line >> atomEntry ; 
  }
  
  PDB.seekg( position ); return PDB ;
}
/*=================================================================*/
/*return a pointer to appropriate doContact member function*/
 bool (PDBamino_acid::*sel_doCont(cont_sel type))(PDBamino_acid &,double &){
  switch(type){/*select type of contact*/
  case _CA_CA_:
    cout<<"ERROR: from sel_doCont(..)\nNot yet implemented!\n";exit(0);
  case _CB_CB_:
    cout<<"ERROR: from sel_doCont(..)\nNot yet implemented!\n";exit(0);
  case _HV_HV_:
    return &PDBamino_acid::doContact;
  case _HVSC_HVSC_: /*contact between heavy atoms in the side-chain*/
    return &PDBamino_acid::doContactResidues;
  case _HVBK_HVBK_: /*contact between heavy atoms in the backbone*/
    return &PDBamino_acid::doContactBackbones;
    /*contact between heavy atoms, one in the backbone, the other in the 
     residue. Order is unimportant*/
  case _HVBK_HVSC_: 
    return &PDBamino_acid::doContactResBack;
    /*contact between heavy atoms, one in the backbone, the other in the 
      residue. The atom in the backbone belongs to first amino acid.*/
  case _F_HVBK_S_HVSC_:
    return &PDBamino_acid::doContactFirstBackSecondRes;
    /*contact between heavy atoms, one in the backbone, the other in the 
      residue. The atom in the residue belongs to first amino acid.*/
  case _F_HVSC_S_HVBK_: 
    return &PDBamino_acid::doContactFirstResSecondBack;
  case _HN_HN_:
    return &PDBamino_acid::doContactHN_HN;
  case _METHYLC_METHYLC_:
    return &PDBamino_acid::doContactMETHYLC_METHYLC;
  }
}
/*=================================================================*/
PDBamino_acid::PDBamino_acid( ){ name.assign( " " ) ; resSeq = 0 ;}
/*=================================================================*/
PDBamino_acid::PDBamino_acid( string &res_name ){
  string sN(" N  "), sCA(" CA "), sC(" C  "), sO(" O  ") ;
  /*  string sCB(" CB ");*/
  double theta,phi;
  PDBvector n (N_CA*cos(init_N_CA_x*rpi), -N_CA*sin(init_N_CA_x*rpi), 0.0);
  PDBvector ca( 0.000,  0.000,  0.000);
  PDBvector c (CA_C*cos(init_C_CA_x*rpi), CA_C*sin(init_C_CA_x*rpi),0.0); 
  theta = init_C_CA_x+CA_C_O-180.0; 
  PDBvector o (CA_C*cos(init_C_CA_x*rpi) + C_O*cos(theta*rpi),
	       CA_C*sin(init_C_CA_x*rpi) + C_O*sin(theta*rpi), 0.0);

  theta=(180-CB_N_CA_C)*rpi; 
  phi = (360.0 - N_CA_CBXY - init_N_CA_x)*rpi;


  PDBvector cb(CA_CB*cos(theta)*cos(phi), 
	       CA_CB*cos(theta)*sin(phi), CA_CB*sin(theta));

  PDBatom N (n) ;  N.assignName(sN );  N.assignResName(res_name);
  PDBatom CA(ca); CA.assignName(sCA); CA.assignResName(res_name);
  PDBatom C (c) ;  C.assignName(sC );  C.assignResName(res_name);
  PDBatom O (o) ;  O.assignName(sO );  O.assignResName(res_name);
  /*PDBatom CB(cb); CB.assignName(sCB); CB.assignResName(res_name);*/

  name = res_name ;
  resSeq = 0 ;

  backbone.insertAtBack(N) ;
  backbone.insertAtBack(CA);
  backbone.insertAtBack(C) ;
  backbone.insertAtBack(O) ;

  PDBatom sd ;
  for(int i=0; i<Nside_chain_atoms; i++ ){
    if( side_chain_atoms[i].find( res_name ) != string::npos ){
      side_chain_atoms[i] >> sd ;
      residue.insertAtBack( sd ) ;
    }
  }
  /*residue.insertAtBack(CB) ;*/
}

PDBamino_acid::~PDBamino_acid( ){ ( *this ).empty( ) ; }

bool PDBamino_acid::isEmpty(){
  if( backbone.isEmpty() && residue.isEmpty() ){ return true ; }
  else return false ;
}

void PDBamino_acid::empty( ){
  if( !( *this ).isEmpty( ) ){ 
    name.assign( " " ) ; resSeq = 0 ;
    backbone.empty( ) ; residue.empty( ) ; 
  }
}

void PDBamino_acid::assign( listPDBatom &theBackbone , 
			    listPDBatom &theResidue ){
  backbone = theBackbone ; residue = theResidue ;
}

void PDBamino_acid::assign( const string &theName ){ name = theName ; }
/*=======================================================================*/
string *PDBamino_acid::printName( ){
  string *tmp = new string ;
  ( *tmp ).assign( name ) ;
  return tmp ;
}
/*=======================================================================*/
string PDBamino_acid::getName( ){ return name ; }
/*=======================================================================*/
PDBamino_acid &PDBamino_acid::operator=( PDBamino_acid &aa ){
  name = aa.name ; resSeq = aa.resSeq ;
  backbone = aa.backbone ; residue = aa.residue ;
  return *this ;
}

void PDBamino_acid::addAtom( PDBatom &entryAtom ){
  if( ( *this ).isEmpty( ) ){ 
    resSeq = entryAtom.getResSeq( ) ; 
    string *tmpString ; tmpString = entryAtom.printResName( ) ;
    name = *tmpString ;
  }
  if( entryAtom.is_backbone( ) ){ backbone.insertAtBack( entryAtom ) ; }
  else{ residue.insertAtBack( entryAtom ) ; }
}

bool PDBamino_acid::matchResSeq( PDBatom &entryAtom ){
  if( resSeq == entryAtom.getResSeq( ) ){ return true ; }
  else{ return false ; }
}

string *PDBamino_acid::printLastAtom( ){
  if( residue.isEmpty( ) ){
    return  backbone.printAtomAt( backbone.length( ) ) ;
  }
  else{
    return residue.printAtomAt( residue.length( ) );
  }
}

listPDBatom *PDBamino_acid::getResidue( ){
  listPDBatom *tmpList = new listPDBatom ; 
  if( !(*this).isEmpty( ) ){  *tmpList = residue ; }
  return tmpList ;
}
/*=======================================================================*/
bool PDBamino_acid::doContactResidues( PDBamino_acid &otherAA , 
				       double &cutOff ) {
  return residue.doContactHeavyAtom( otherAA.residue , cutOff ) ;
}
/*=======================================================================*/
bool PDBamino_acid::doContactBackbones( PDBamino_acid &otherAA , 
				        double &cutOff ) {
  return backbone.doContactHeavyAtom( otherAA.backbone , cutOff ) ;
}
/*=======================================================================*/
bool PDBamino_acid::doContactFirstResSecondBack( PDBamino_acid &otherAA , 
				      double &cutOff ) {
  /*printf("bool PDBamino_acid::doContactFirstResSecondBack(..)\n");*/
  return residue.doContactHeavyAtom( otherAA.backbone , cutOff ) ;
}
/*=======================================================================*/
bool PDBamino_acid::doContactFirstBackSecondRes( PDBamino_acid &otherAA , 
				      double &cutOff ) {
  /*printf("bool PDBamino_acid::doContactFirstBackSecondRes(..)\n");*/
  return backbone.doContactHeavyAtom( otherAA.residue , cutOff );
}
/*=======================================================================*/
bool PDBamino_acid::doContactResBack( PDBamino_acid &otherAA, double &cutOff){
  return ( backbone.doContactHeavyAtom( otherAA.residue , cutOff ) ||
	   residue.doContactHeavyAtom( otherAA.backbone , cutOff ) ) ;
}
/*=========================================================================*/
bool PDBamino_acid::doContactHN_HN( PDBamino_acid &otherAA, double &cutOff){
  PDBatom *hn1=backbone.pointToAtomFromName(" H  "),
    *hn2=otherAA.pointToAtomFromName(" H  ");
  if(hn1 && hn2) return hn1->doContact(*hn2,cutOff);
  return 0;
}
/*=========================================================================*/
bool PDBamino_acid::doContactMETHYLC_METHYLC( PDBamino_acid &otherAA, double &cutOff){
  string methyls("ALA VAL LEU ILE");
  string name1=this->getName( ), name2=otherAA.getName();
  if(methyls.find(name1)!=string::npos && methyls.find(name2)!=string::npos){
    listPDBatom l1,l2;
    if( this->extract(_METHYLC_, l1)==0 || otherAA.extract(_METHYLC_, l2)==0 ){
      cerr<<"name1=|"<<name1<<"|\n"<<*this;
      cerr<<"name2=|"<<name2<<"|\n"<<otherAA;
      cerr<<"ERROR in doContactMETHYLC_METHYLC: no carbons to extract\n";
    }
    return l1.doContactAnyone(l2,cutOff);
  }
  return false;
}
/*=========================================================================*/
bool PDBamino_acid::doContact( PDBamino_acid &otherAA , 
			       double &cutOff ) {
  return ( ( *this ).doContactResidues( otherAA , cutOff ) ||
	   ( *this ).doContactBackbones( otherAA , cutOff ) ||
	   ( *this ).doContactResBack( otherAA , cutOff )     ) ;
}
/*=========================================================================*/
int PDBamino_acid::doContact( PDBamino_acid &otherAA, gen_atom_type type1,
			       gen_atom_type type2, double &cutOff){
  /*cout<<"type1=\""<<gat2str(type1)<<"\" type2=\""<<gat2str(type2)<<"\"\n";*/
  listPDBatom l1,l2;
  if( this->extract(type1, l1)!=0 && otherAA.extract(type2, l2)!=0 )
    return l1.numberAnyContacts(l2,cutOff);
  return 0;
}
/*=========================================================================*/
int PDBamino_acid::getResSeq( ){ return resSeq ; }
/*=========================================================================*/
void PDBamino_acid::changeFullyResSeq( const int &theResSeq ){
  resSeq = theResSeq ;
  backbone.changeResSeq( theResSeq ) ;
  residue.changeResSeq( theResSeq ) ;
}
/*=========================================================================*/
bool PDBamino_acid::is( const string &type ){ return name == type ; }
/*================================================================*/
PDBatom *PDBamino_acid::getAtomFromName( const string &theName ){
  PDBatom *tmp ;
  string tmpString ;
  if( backbone.isThereAtom( theName ) ){
    tmp = backbone.getAtomFromName( theName ) ;
  }
  else if( residue.isThereAtom( theName ) ){
    tmp = residue.getAtomFromName( theName ) ;
  }
  else{ 
    tmp = new PDBatom ;
  }
  
  return tmp ;
}
/*================================================================*/
PDBatom *PDBamino_acid::pointToAtomFromName( const string &theName ){
  string tmpString ;
  if( backbone.isThereAtom( theName ) ){
    return backbone.pointToAtomFromName( theName ) ;
  }
  else if( residue.isThereAtom( theName ) ){
    return residue.pointToAtomFromName( theName ) ;
  }
  return NULL ;
}
/*================================================================*/
PDBatom PDBamino_acid::getAtomWithIndex( const int &theIndex ){
  PDBatom tmp ;
  if( backbone.hasAtomWithAtomIndex( theIndex ) ){
    tmp = backbone.getAtomWithAtomIndex( theIndex ) ;
  }
  else if( residue.hasAtomWithAtomIndex( theIndex ) ){
    tmp = residue.getAtomWithAtomIndex( theIndex ) ;
  }
  return tmp ;
}
/*================================================================*/
bool PDBamino_acid::hasAtomWithIndex( const int &theIndex ){
  if( backbone.hasAtomWithAtomIndex( theIndex ) ){ return true ; }
  else if( residue.hasAtomWithAtomIndex( theIndex ) ){ return true ; }
  return false ;
}
/*================================================================*/
void PDBamino_acid::printChainID( string &theID ){
  if(!backbone.isEmpty()){ backbone.printChainID( theID ); return; }
  residue.printChainID(theID);
}
/*================================================================*/
bool PDBamino_acid::hasAtomWithName( const string &theName ){
  if( backbone.isThereAtom( theName ) || residue.isThereAtom( theName )){
    return true ;
  }
  return false ;
}

bool PDBamino_acid::changeCoordOfAtomName( const string &theName,
					   double *coord ){
  if( backbone.isThereAtom( theName ) ){
    return backbone.changeCoordOfAtomName( theName , coord ) ;
  }
  else if( residue.isThereAtom( theName ) ){
    residue.changeCoordOfAtomName( theName , coord ) ;
  }
  else{ return false ; }  
}
/*================================================================*/
double * PDBamino_acid::getCoordOfAtomName( const string &AtomName ){
  if( backbone.isThereAtom( AtomName ) ){
    return backbone.getCoordOfAtomName( AtomName ) ;
  }
  else if( residue.isThereAtom( AtomName ) ){
    return residue.getCoordOfAtomName( AtomName ) ;
  }
  else{ return NULL ; }
}
/*================================================================*/
void PDBamino_acid::getCoordOfAtomName( const string &AtomName, PDBvector &the_vec ){
  the_vec.assign( this->getCoordOfAtomName(AtomName) ) ;
}
/*================================================================*/
bool PDBamino_acid::getCoordOfAtomWithAtomIndex( const int &index, 
						 double *coords ){
  if( this->backbone.getCoordOfAtomWithAtomIndex( index, coords ) ||
      this->residue.getCoordOfAtomWithAtomIndex( index, coords ) ){ return true; }
  return false ;
}
/*================================================================*/
bool PDBamino_acid::changeCoordOfAtomWithIndex( const int &index, double *newCoord ){
  PDBatom *tmp ;
  listPDBatom *list=NULL ;
  if( this->backbone.hasAtomWithAtomIndex( index ) ){
    list = &backbone ;
  }
  else if( this->residue.hasAtomWithAtomIndex( index ) ){
    list = &residue ;
  }
  if( list ){
    tmp = list->returnAdressOfAtomWithAtomIndex( index ) ;
    tmp->assignCoord( newCoord ) ;
    return true ;
  }
  return false ;
}
/*================================================================*/
bool PDBamino_acid::getNameOfAtomWithAtomIndex( const int &index,
						string &atom_name ){
  if( this->backbone.getNameOfAtomWithAtomIndex( index, atom_name ) ||
      this->residue.getNameOfAtomWithAtomIndex( index, atom_name)){ return true; }
  return false ;
}
/*================================================================*/
int PDBamino_acid::numberOfAtoms( ){
  return backbone.length( ) + residue.length( ) ;
}
/*================================================================*/
int PDBamino_acid::renumberFullyAtomSerialNumber( const int &firstIndex ){
  int index = backbone.renumberFullyAtomSerialNumber( firstIndex ) ;
  return residue.renumberFullyAtomSerialNumber( index ) ;
}

int PDBamino_acid::getSerialNumberOfAtomName( const string &AtomName ){
  if( backbone.isThereAtom( AtomName ) ){
    return backbone.getSerialNumberOfAtomName( AtomName ) ;
  }
  else if( residue.isThereAtom( AtomName ) ){
    return residue.getSerialNumberOfAtomName( AtomName ) ;
  }
  else{ return 0 ; }
}

bool PDBamino_acid::bonds_its_amide_H_to_the_backbone_O_of( PDBamino_acid &other_aa ){
  const string H_name=" H  " ;
  if(this->hasAtomWithName( H_name ) ){

    const string O_name=" O  ";
    PDBatom *H=this->getAtomFromName( H_name ) ;
    PDBatom *O=other_aa.getAtomFromName( O_name ) ;
    double d2 = H->d2(*O);
    //printf("d_HO=%lf\n", sqrt(d2));
    if( d2<2.89 || d2>4.41) return false ; //1.9A<distance<2.1A

    const string N_name=" N  " ;
    PDBatom *N=this->getAtomFromName( N_name ) ;
    double angle=H->angle_with(*N,*O);
    //printf("ang_NHO=%lf\n",angle*180/PI);
    if( angle<2.443461 ) return false ; //140<angle<180

    const string C_name=" C  " ;
    PDBatom *C=other_aa.getAtomFromName( C_name ) ;
    angle=O->angle_with(*H,*C);
    //printf("ang_OHC=%lf\n",angle*180/PI);
    //if( angle<2.268928 || angle>2.9670597 ) //130<angle<170
    if(angle<2.268928)//130<angle<180
      return false ;
    else return true ;
  }
  return false;
}
/*============================================================*/
void PDBamino_acid::translate( const PDBvector &shift ){
  backbone.translate(shift);
  residue.translate(shift);
}
/*============================================================*/
void PDBamino_acid::rotate( const PDBvector &axis, const double &angle ){
  backbone.rotate( axis, angle );
  residue.rotate( axis, angle );
}
/*============================================================*/
void PDBamino_acid::rotate( const PDBvector &new_origin, const PDBvector &axis,
			    const double &angle ){
  backbone.rotate( new_origin, axis, angle );
  residue.rotate( new_origin, axis, angle );
}
/*============================================================*/
void PDBamino_acid::rotate( double ** rotation_matrix ){
  backbone.rotate( rotation_matrix );
  residue.rotate( rotation_matrix );
}
/*============================================================*/
PDBvector PDBamino_acid::get_next_N_coords(void){
  PDBvector ca, c, o;
  this->getCoordOfAtomName(" CA ", ca ) ;
  this->getCoordOfAtomName(" C  ", c  ) ;
  this->getCoordOfAtomName(" O  ", o  ) ;

  PDBvector c_ca = ca - c ;  c_ca.normalize( ) ;
  PDBvector c_o  =  o - c ;   c_o.normalize( ) ;

  double CA_C_N = 116.2*PI/180.0, CA_C_O = 120.8*PI/180.0;
  double O_C_N = 2*PI - CA_C_N - CA_C_O;
  double a = (cos(CA_C_N)-cos(O_C_N)*cos(CA_C_O))/(sin(CA_C_O)*sin(CA_C_O));
  double b = (cos(O_C_N)-cos(CA_C_N)*cos(CA_C_O))/(sin(CA_C_O)*sin(CA_C_O));

  double C_N = 1.329;
  return c + (c_ca * a + c_o * b) * C_N ;
}
/*============================================================*/
PDBvector PDBamino_acid::get_nextCA(void){
  PDBvector tmp, x, y, z;
  this->peptidePlaneAxis(x, y, z);
  double theta = (C_CA0_CA + init_C_CA_x)*rpi;
  this->getCoordOfAtomName(" CA ", tmp);
  tmp += x*(CA_CA*cos(theta)) +y*(CA_CA*sin(theta));
  return tmp;
}
/*============================================================*/
void PDBamino_acid::peptidePlaneAxis(PDBvector &x_axis, PDBvector &y_axis,
				     PDBvector &z_axis){
  PDBvector o, ca, c, c_o, c_ca;

  this->getCoordOfAtomName(" O  ", o);
  this->getCoordOfAtomName(" CA ", ca);
  this->getCoordOfAtomName(" C  ", c);

  c_o  = o  - c ;
  c_ca = ca - c ;

  c_o.normalize( ) ;
  c_ca.normalize( );

  y_axis = c_o*(-sin(init_y_C_CA*PI/180.0)/sin(CA_C_O*PI/180.0)) +
    c_ca*(-sin(init_O_C_y*PI/180.0)/sin(CA_C_O*PI/180.0));
  z_axis= c_ca^c_o;
  z_axis.normalize();
  x_axis=y_axis^z_axis;
  return;
}
/*============================================================*/
void PDBamino_acid::intrinsicAxis(PDBvector &x_axis, PDBvector &y_axis,
				  PDBvector &z_axis){
/*Determine the intrinsic axis of the amino acids
  FOR INITIALIZING AMIN ACID STRUCTURE
  
  Two Streched amino acids(psi = -180, and phi = 180)
                     (ix)
       Ca--(iy)   N   |    C
     / a|  \    /  \ a|  /
    N   |    C       Ca--(iy)
       (ix)   
  angle : a is such that the two intrinsic x axis are parallel
  ------------------------------*/
  PDBvector n, ca, c, ca_n, ca_c;;
  this->getCoordOfAtomName(" N  ", n);
  this->getCoordOfAtomName(" CA ", ca);
  this->getCoordOfAtomName(" C  ", c);

  ca_n = n - ca;
  ca_c = c - ca; 

  ca_n.normalize( ) ;
  ca_c.normalize( ) ;

  double n_ca_c = N_CA_C*rpi;
  double n_ca_x = init_N_CA_x*rpi;
  double c_ca_x = init_C_CA_x*rpi;

  double dist_ca_n = sin(c_ca_x)/sin(n_ca_c);
  double dist_ca_c = sin(n_ca_x)/sin(n_ca_c);

  x_axis = ca_n*dist_ca_n + ca_c*dist_ca_c;  
  x_axis.normalize();
  z_axis = ca_n^ca_c;
  z_axis.normalize();
  y_axis = z_axis^x_axis;
}
/*============================================================*/
bool PDBamino_acid::align_intrinsic_axis_to_system_axis( ){
/*Two Streched amino acids(psi = -180, and phi = 180)
                     (ix)
       Ca--(iy)   N   |    C
     / a|  \    /  \ a|  /
    N   |    C       Ca--(iy)
       (ix)   
  angle : a is such that the two intrinsic x axis are parallel
  ------------------------------*/
  if( !(this->hasAtomWithName(" C  ")) || !(this->hasAtomWithName(" CA "))
      || !(this->hasAtomWithName(" N  ")) ){ return false ; }
  
  /*put CA in origin*/
  PDBatom *CA = this->getAtomFromName( " CA " ) ;
  PDBvector vCA = CA->getCoord( ) ;
  vCA *= -1 ;  this->translate( vCA ) ;
  /*rotate to align system and intrinsic axes*/
  PDBvector xi, yi, zi ; 
  this->intrinsicAxis( xi, yi , zi ) ;
  double **rot = alloc_array( 3, 3 ) ;
  rot[0][0] = xi[0] ;  rot[0][1] = xi[1] ;  rot[0][2] = xi[2] ;
  rot[1][0] = yi[0] ;  rot[1][1] = yi[1] ;  rot[1][2] = yi[2] ;
  rot[2][0] = zi[0] ;  rot[2][1] = zi[1] ;  rot[2][2] = zi[2] ;
  this->rotate( rot ) ;
  /*take back translation*/
  vCA *= -1 ;  this->translate( vCA ) ;
  return true ;
}
/*============================================================*/
void PDBamino_acid::join_to( PDBamino_acid &prev_aa ){
  PDBvector ca, prev_ca;
  this->getCoordOfAtomName(" CA ", ca);
  prev_aa.getCoordOfAtomName(" CA ", prev_ca);
  PDBvector shift = prev_aa.get_nextCA( ) - ca ;
  
  //shift *= -1 ;
  PDBvector axisbar[3];
  prev_aa.peptidePlaneAxis(axisbar[0], axisbar[1], axisbar[2] );
  axisbar[0] *= -1.0;
  axisbar[2] *= -1.0;
  double** rotation_matrix = alloc_array( 3, 3);
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){ rotation_matrix[i][j] = axisbar[i][j]; }
  }
  this->rotate( rotation_matrix ) ;
  this->translate( shift ) ;
}
/*============================================================*/
void PDBamino_acid::filterOnly( const string *names, const int &n ){
  backbone.filterOnly( names, n ) ;
  residue.filterOnly( names, n ) ;
}
/*============================================================*/
double PDBamino_acid::minD( PDBamino_acid &neighborAA ){
  /*cout<<"double PDBamino_acid::minD( PDBamino_acid & )\n";*/
  double mind = residue.minD( neighborAA.residue ) ;
  double d = backbone.minD( neighborAA.residue ) ; if( d< mind ){ mind = d ; }
  d = residue.minD( neighborAA.backbone ) ; if( d< mind ){ mind = d ; }
  d = backbone.minD( neighborAA.backbone ) ; if( d< mind ){ mind = d ; }
  /*cout<<"mind="<<mind<<endl;*/
  return mind ;
}
/*============================================================*/
void PDBamino_acid::remove_all_hydrogens( ){
  this->backbone.remove_all_hydrogens( ) ;
  this->residue.remove_all_hydrogens( ) ;
}
/*============================================================*/
void PDBamino_acid::fuseToSingleList( listPDBatom &theList ){
  if( !theList.isEmpty( ) ){ theList.empty( ) ; }
  theList.insertAtBack( backbone ) ;
  theList.insertAtBack( residue  ) ;
}
/*============================================================*/
int PDBamino_acid::extract(gen_atom_type type, listPDBatom &list){
  if(!list.isEmpty()) list.empty();
  return backbone.extract(type,list) + residue.extract(type,list);
}
/*============================================================*/
/*============================================================*/
/*============================================================*/
/*============================================================*/
/*============================================================*/


//------------ class nodePDBamino_acid  --------------
nodePDBamino_acid::nodePDBamino_acid( PDBamino_acid &aa , 
				      const int &theResSeq){
  theAA =  aa ; resSeq = theResSeq ; next = NULL ;
}

void nodePDBamino_acid::changeFullyResSeq( const int &theResSeq ){
  resSeq = theResSeq ;
  theAA.changeFullyResSeq( theResSeq ) ;
}

//------------ class PDBchain ------------------------

ostream &operator<<(ostream &output, const PDBchain &theList ){
  nodePDBamino_acid *currentPtr = theList.firstPtr ;
  while( currentPtr != NULL){
    output << (*currentPtr).theAA ;
    currentPtr = (*currentPtr).next ;
  }
  return output ;
}

ifstream &operator>>(ifstream &input, PDBchain &the_chain ){
  the_chain.importChain( input ) ;
  return input ;
}

PDBchain::PDBchain( ){ firstPtr = lastPtr = NULL ; }
/*===================================================================*/
PDBchain::PDBchain( ifstream &PDB ){ 
  string line ;
  int position, resSeq=1 ;
  PDBamino_acid theAminoAcid ;

  firstPtr = lastPtr = NULL ; 
  do{ position=PDB.tellg( ) ; getline( PDB, line ) ; }while( !is_ATOM(line) ) ;
  while( !is_TER(line) && !is_END(line) && !PDB.eof() ){
    PDB.seekg( position ) ; //go back one line
    PDB >> theAminoAcid ; 
    //cout << theAminoAcid << endl ; 
    (*this).insertAtBack( theAminoAcid , resSeq ) ; resSeq++ ;
    position = PDB.tellg( ) ; getline( PDB, line ) ;
  }
}
/*===================================================================*/
PDBchain::PDBchain( char *file ){
  firstPtr = lastPtr = NULL ;
  /*in case file is zipped*/
  char uzfile[512];
  zip_format fmt=returnZipFormat(file);  unzip(uzfile,file);

  ifstream PDB(uzfile); PDB>>*this;
  PDB.close();  zip(uzfile,fmt);
}
/*===================================================================*/
PDBchain::~PDBchain( ){ ( *this ).empty( ) ; }

nodePDBamino_acid *PDBchain::getNewNode( PDBamino_acid &aa ,
					 const int &theResSeq ){
  nodePDBamino_acid *ptr = new nodePDBamino_acid( aa , theResSeq );
  assert( ptr != 0 ) ;
  return ptr ; 
}

void PDBchain::empty( ){
 if( !this->isEmpty( ) ){
   while( firstPtr != NULL ){ ( *this ).rmFromFront( ); }
 }
}

bool PDBchain::isEmpty( ){ return firstPtr == NULL ; }

void PDBchain::rmFromFront( ){
  nodePDBamino_acid *currentPtr = firstPtr ;
  firstPtr = ( *currentPtr ).next ; delete currentPtr ;
}
/*================================================================*/
void PDBchain::insertAtBack( PDBamino_acid &aa , const int &theResSeq ){
  nodePDBamino_acid *newPtr = getNewNode( aa , theResSeq ) ;
  if( this->isEmpty() ){ firstPtr = lastPtr = newPtr ; }
  else{ ( *lastPtr ).next = newPtr ; lastPtr = newPtr ; }
}
/*================================================================*/
void PDBchain::insertAtBack( PDBamino_acid &aa , 
			     const double prev_psi, const double  phi){
  PDBamino_acid tmp;
  tmp = aa ;
  int theResSeq = 1 ;
  if( !this->isEmpty( ) ){ 
    theResSeq = lastPtr->resSeq + 1 ; 
    /*adapt tmp to lastPtr->theAA, but do NOT add it yet to the chain*/
    tmp.join_to( lastPtr->theAA );
  }
  nodePDBamino_acid *newPtr = getNewNode( tmp , theResSeq ) ;
  int l = this->length( ) ;/*lenght before adding the amino acid*/
  if( this->isEmpty( ) ){ firstPtr = lastPtr = newPtr ; }
  else{ 
    lastPtr->next = newPtr ; lastPtr = newPtr ;
    /*first Psi, then Phi*/
    this->rotatePsiBy( prev_psi, l ) ;
  }
  this->rotatePhiBy( phi, l+1 ) ;
}
/*================================================================*/
int PDBchain::length( ){
  int count = 0; nodePDBamino_acid  *currentPtr = firstPtr ;
  while( currentPtr != NULL ){ 
    count++ ;  currentPtr = (*currentPtr).next ;
  }
  return count ;
}
/*=====================================================*/
int PDBchain::importChain( ifstream &PDB ){  
  string line,chainID,old_chainID("@");
  bool f=false;
  int position, resSeq=1 ;
  PDBamino_acid theAminoAcid ;
  
  (*this).empty(); /*empty the chain*/
  do{//go to the first amino acid, unless TER or END 
    position=PDB.tellg( ) ; 
    getline( PDB, line ) ;
    if(is_TER(line) || is_END(line)){/*the chain is empty*/
      PDB.seekg( position ) ;
      return 1;
    }
  }while( !( is_ATOM(line) && is_AMINOACID(line) ) && !( PDB.eof() ) ) ; 

  if( PDB.eof() ){ return 1 ; }
  else{ PDB.seekg( position ) ; } //go back one line to read the amino acid

  do{
    PDB >> theAminoAcid ;
    if( !theAminoAcid.isEmpty() ){
      theAminoAcid.printChainID(chainID);
      if(chainID != old_chainID){
	if(f){ PDB.seekg( position ) ; break ; }
	else{ old_chainID=chainID; f=true; }/*initialization*/
      }
      (*this).insertAtBack( theAminoAcid , resSeq ) ; 
      resSeq++ ; 
    }
    if( PDB.eof() ){ return 1 ; }
    position=PDB.tellg( ) ;
    getline( PDB, line ) ;
    /*cout<<line<<endl;*/
    if( is_TER(line) ){       /*cout <<"is TER!\n" ;*/
      PDB.seekg( position ) ; break ; 
    }
    if( is_END(line) ){      /*cout <<"is END!\n" ;*/
      PDB.seekg( position ) ; break ; 
    }
    PDB.seekg( position ) ;
  }while(1) ;

  //cout << "out of the do loop\n" ; 
  return 0 ;
}
/*=====================================================*/
bool PDBchain::importChain(char *file,const string &chainID){
  PDBchains chains(file); /*cout<<chains;exit(0);*/
  string none("_");
  if(chainID==none){/*only one chain*/
    /*cout<<"Only one chain\n";*/
    return chains.getPDBchainFromIndex(*this,1);
  }
  else{
    /*cout<<"More than one chain\n";cout<<chains;*/
    return chains.getPDBchainFromChainID(*this,chainID);/*cout<<*this;*/
  }
}
/*=====================================================*/
PDBchain &PDBchain::operator=(PDBchain &right){
  if(this==&right){cout<<"ERROR:self-assignment!\n";exit(0);}
  this->empty();
  nodePDBamino_acid *currentPtr=right.firstPtr ;
  int i=1;
  while(currentPtr){
    this->insertAtBack(currentPtr->theAA,i); i++;
    currentPtr=currentPtr->next;
  }
  return *this;
}
/*=====================================================*/
PDBamino_acid *PDBchain::getAminoFromIndex( const int &index ){
  PDBamino_acid *tmpAmino = NULL ;
  if( (index > this->length( )) || (index<1) ){ return tmpAmino ; }
  nodePDBamino_acid *currentPtr = firstPtr ; 
  int count = 1 ;
  while( count < index ){ currentPtr = currentPtr->next ;  count++ ; }
  tmpAmino = new PDBamino_acid ;  *tmpAmino = currentPtr->theAA ;
  return tmpAmino ;
}
/*=====================================================*/
PDBamino_acid *PDBchain::pointToAminoFromIndex( const int &index ){
  if( (index > this->length( )) || (index<1) ){ return NULL ; }
  nodePDBamino_acid *currentPtr = firstPtr ; 
  int count = 1 ;
  while( count < index ){ currentPtr = currentPtr->next ;  count++ ; }
  return &(currentPtr->theAA);
}
/*=====================================================*/
PDBamino_acid *PDBchain::getAminoFromResSeq( const int &theResSeq ){
  PDBamino_acid *tmpPtr = new PDBamino_acid ;
  nodePDBamino_acid *currentPtr = firstPtr ; 
  while( currentPtr != NULL ){
    if( ( *currentPtr ).resSeq == theResSeq ){
      *tmpPtr = ( *currentPtr ).theAA ;
      return tmpPtr ;
    }
    currentPtr = (*currentPtr).next ;
  }
  return tmpPtr ;
}

string *PDBchain::printTER( ){
  string *tmpString = new string( " " ) ;
  if( ( *this ).isEmpty( ) ){ return tmpString ; }
  PDBamino_acid *tmpAmino = NULL ;
  tmpAmino = ( *this ).getAminoFromIndex( ( *this ).length( ) ) ;
  tmpString = ( *tmpAmino ).printLastAtom( ) ;
  PDBatom tmpAtom ; *tmpString >> tmpAtom ;
  return tmpAtom.printTER( ) ;
}
/*=====================================================================*/
string PDBchain::printChainID(){
  string id(" ");
  firstPtr->theAA.printChainID(id);
  return id;
}
/*========================================================================*/
int PDBchain::getIndexFromResSeq( const int &theResSeq ){
  nodePDBamino_acid *currentPtr = firstPtr ; 
  int count = 0 ;
  while( currentPtr != NULL ){
    count++;
    if( ( *currentPtr ).resSeq == theResSeq ){ return count ; }
    currentPtr = (*currentPtr).next ;
  }
  return count ;
}
/*=======================================================================*/
int PDBchain::getResSeqFromIndex( const int &index ){
  nodePDBamino_acid *currentPtr = firstPtr ;
  int count = 1 ;
  while( (count < index) && (currentPtr != NULL) ){ 
    currentPtr = (*currentPtr).next ;
    count++ ; //must go after previous statement, not before
  }
  return ( *currentPtr ).resSeq ;
}
/*=============================================================*/
string *PDBchain::printNameFromIndex( const int &index ){
  string *tmpString ;
  PDBamino_acid *tmpAmino ;
  tmpAmino = getAminoFromIndex( index ) ;
  tmpString = tmpAmino->printName( ) ;
  return tmpString ;
}
/*=============================================================*/
string PDBchain::getNameFromIndex( const int &index ){
  if( index<1 || index>this->length( ) ){
    cerr<<"ERROR from PDBchain::getNameFromIndex( const int &):\n";
    cerr<<"     residue index out of bounds\n";
  }
  PDBamino_acid *tmpAmino ;
  tmpAmino = getAminoFromIndex( index ) ;
  return tmpAmino->getName( ) ;
}
/*=============================================================*/
void PDBchain::renumberFully( const int &firstIndex ){
  nodePDBamino_acid *currentPtr = firstPtr ;
  int index = firstIndex ;
  while( currentPtr != NULL ){
    ( *currentPtr ).changeFullyResSeq( index ) ;
    index++ ;
    currentPtr = (*currentPtr).next ;
  }
}

int PDBchain::renumberFullyAtomSerialNumber( const int &firstIndex ){
  nodePDBamino_acid *currentPtr = firstPtr ;
  int index = firstIndex ;
  while( currentPtr ){
    index = ( ( *currentPtr ).theAA ).renumberFullyAtomSerialNumber( index ) ;
    currentPtr = ( *currentPtr ).next ;
  }
  return index ;
}

bool PDBchain::isResNameAtIndex( const string &theResName, 
				 const int &theIndex ){
  nodePDBamino_acid *currentPtr = firstPtr ;
  int index=1 ;
  while( index < theIndex ){ currentPtr = (*currentPtr).next ; index++ ; }
  return ( (*currentPtr).theAA ).is( theResName ) ;
}

bool PDBchain::changeCoordOfAtomNameAtIndex( const string &AtomName, 
					     const int &theIndex ,
					     double * coord ){
  if( (*this).length( ) < theIndex ){ return false ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  int index=1 ;
  while( index < theIndex ){ currentPtr = (*currentPtr).next ; index++ ; }
  if(!(  ( (*currentPtr).theAA ).hasAtomWithName( AtomName ) )){ return false; }
  ( (*currentPtr).theAA ).changeCoordOfAtomName( AtomName , coord ) ;
  return true ;
}

double *PDBchain::getCoordOfAtomNameAtIndex( const string &AtomName, 
					     const int &theIndex){
  if( (*this).length( ) < theIndex ){ return NULL ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  int index=1 ;
  while( index < theIndex ){ currentPtr = (*currentPtr).next ; index++ ; }
  if(!(  ( (*currentPtr).theAA ).hasAtomWithName( AtomName ) )){ return NULL ; }
  return ( (*currentPtr).theAA ).getCoordOfAtomName( AtomName ) ;
}
/*================================================================*/
bool PDBchain::getCoordOfAtomNameAtIndex_II( const string &AtomName,
					  const int &theIndex,
					  double *r ){
  double *tmp = (*this).getCoordOfAtomNameAtIndex( AtomName, theIndex ) ;
  if( r == NULL || tmp== NULL ){ return false ; }
  else{ 
    for( int i=0 ; i<3 ; i++ ){ r[i] = tmp[i] ; }
    return true ;
  }
}
/*================================================================*/
bool PDBchain::getCoordOfAtomWithAtomIndex( const int &index, double *coords ){
  if( this->isEmpty( ) ){ return false ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    if( currentPtr->theAA.getCoordOfAtomWithAtomIndex( index, coords ) ){
      return true ;
    }
    currentPtr = currentPtr->next;
  }    
  return false ;
}
/*================================================================*/
int PDBchain::getSerialNumberOfAtomNameAtIndex( const string &AtomName,
						const int &theIndex ){
  if( (*this).length( ) < theIndex ){ return 0 ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  int index=1 ;
  while( index < theIndex ){ currentPtr = (*currentPtr).next ; index++ ; }
  return ( (*currentPtr).theAA ).getSerialNumberOfAtomName( AtomName ) ;
}
/*================================================================*/
PDBatom PDBchain::getAtomWithIndex( const int &theIndex ){
  PDBatom tmp ;
  nodePDBamino_acid *currentPtr = firstPtr ;
  while( currentPtr ){
    if(  currentPtr->theAA.hasAtomWithIndex( theIndex ) ){ 
      return currentPtr->theAA.getAtomWithIndex( theIndex ) ; 
    }
    currentPtr = currentPtr->next ;
  }
  return tmp ;
}
/*================================================================*/
bool PDBchain::changeCoordOfAtomWithIndex(  const int &theIndex, 
					    double * newCoord ){
  PDBatom tmp ;
  nodePDBamino_acid *currentPtr = firstPtr ;
  while( currentPtr ){
    if(  currentPtr->theAA.hasAtomWithIndex( theIndex ) ){ 
      return currentPtr->theAA.changeCoordOfAtomWithIndex( theIndex, newCoord ) ;
    }
    currentPtr = currentPtr->next ;
  }
  return false ;
}
/*================================================================*/
bool PDBchain::changeCoordOfAtomWithIndex( const int &theIndex,
					     PDBvector &theVector ){
  double *newCoord = new double[ 3 ] ;
  for( int i=0; i<3; i++ ){ newCoord[ i ] = theVector[ i ] ; }
  return this->changeCoordOfAtomWithIndex( theIndex, newCoord ) ;
}
/*================================================================*/
PDBatom PDBchain::getAtomWithNameAndResIndex(const char *name,const int &resIdx ){
  PDBamino_acid *aa;
  aa=this->pointToAminoFromIndex(resIdx) ;
  PDBatom *at;
  at=aa->getAtomFromName(name);
  return *at;
}
/*================================================================*/
PDBatom *PDBchain::pointToAtomWithNameAndResIndex(const char *name,const int &resIdx ){
  PDBamino_acid *aa=NULL;
  aa=this->pointToAminoFromIndex(resIdx) ; /*cout<<*aa<<endl;*/
  PDBatom *at=NULL;
  at=aa->pointToAtomFromName(name);
  return at;
}
/*================================================================*/
bool PDBchain::getNameOfAtomWithAtomIndex( const int &index, 
					   string &atom_name ){
  if( this->isEmpty( ) ){ return false ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    if( currentPtr->theAA.getNameOfAtomWithAtomIndex( index, atom_name ) ){
      return true ;
    }
    currentPtr = currentPtr->next;
  }    
  return false ;  
}
/*================================================================*/
int PDBchain::numberOfAtoms( ){
  int numAt = 0 ;
  nodePDBamino_acid *currentPtr = firstPtr ;
  while( currentPtr ){
    numAt += currentPtr->theAA.numberOfAtoms( ) ;
    currentPtr = currentPtr->next ;
  }
  return numAt ;
}
/*================================================================*/
void PDBchain::createCAchain( listPDBatom &theChain ){
  nodePDBamino_acid *currentPtr = firstPtr ;
  PDBatom *theCA ;

  theChain.empty() ;
  while( currentPtr ){
    theCA = ( (*currentPtr).theAA ).getAtomFromName( " CA " ) ;
    theChain.insertAtBack( *theCA );
    delete theCA ; 
    currentPtr = (*currentPtr).next ;   
  }    
}
/*================================================================*/
int PDBchain::createCBChain( listPDBatom &theCBChain ){
  nodePDBamino_acid *currentPtr = firstPtr ;
  PDBatom *theCB ;

  theCBChain.empty() ;
  while( currentPtr ){
    if( ( (*currentPtr).theAA ).is("GLY") ){
      theCB = ( (*currentPtr).theAA ).getAtomFromName( " CA " ) ;
    }
    else{ theCB = ( (*currentPtr).theAA ).getAtomFromName( " CB " ) ; }
    theCBChain.insertAtBack( *theCB );
    delete theCB ;
    currentPtr = (*currentPtr).next ;   
  }
  return 0;
}
/*================================================================*/
int PDBchain::createCBcontactMap( ostream &map, double &cutOff ){
  listPDBatom CBlist ;
  (*this).createCBChain( CBlist ) ;
  return CBlist.createContactMap( map, cutOff ) ;
}
/*================================================================*/
int PDBchain::createCBcontactMap( double **map, double &cutOff ){
  /*cout<<"PDBchain::createCBcontactMap\n";*/
  listPDBatom CBlist ;
  (*this).createCBChain( CBlist ) ;
  /*cout << CBlist<<"TER\n";*/
  return CBlist.createContactMap( map, cutOff ) ;
}
/*================================================================*/
int PDBchain::createCBcontactMap( PDBchain &chain2, double **map,
				  double &cutOff ){
  /*cout<<"PDBchain::createCBcontactMap\n";*/
  listPDBatom CBlist, CBlist2 ;
  (*this).createCBChain( CBlist ) ;
  chain2.createCBChain( CBlist2 ) ;
  /*cout<<CBlist<<"TER\n"<<CBlist2<<"END\n";exit(1);*/
  return CBlist.createContactMap( CBlist2, map, cutOff ) ;
}
/*================================================================*/
/*don't care contact for |i-j|<skip*/
/*map need previous allocation map[0..l-1][0..l-1]*/
int PDBchain::createHeavyAtomContactMap( int &skip,double**map, double &cutoff ){
  if( this->length( ) < skip+2 ){ return 0 ; }
  int n=0,i=0,j=1;
  PDBamino_acid aa;
  nodePDBamino_acid *aa1Ptr = firstPtr;
  nodePDBamino_acid *aa2Ptr = firstPtr->next;
  while(aa1Ptr){
    aa=aa1Ptr->theAA;
    while(aa2Ptr && j-i<=skip){
      map[i][j]=map[j][i]=0;
      aa2Ptr=aa2Ptr->next;/*avoid next neighbors*/
      j++;
    }
    if(j-i<=skip) break;
    while(aa2Ptr){
      if(aa.doContact(aa2Ptr->theAA,cutoff)){
	n++;
	map[i][j]=map[j][i]=1;
      }
      else{
	map[i][j]=map[j][i]=0;
      }
      aa2Ptr=aa2Ptr->next; 
      j++;
    }
    aa1Ptr=aa1Ptr->next;
    i++;
    if(!aa1Ptr) break;/*special case when skip==0*/
    aa2Ptr=aa1Ptr->next;
    j=i+1;
  }
  return n;
}
/*================================================================*/
int PDBchain::createDistContMap( double**map, double &dmin, double &dmax ){
  /*cout<<"int PDBchain::createDistContMap( double**, double&, double& )\n";*/
  int l = this->length( ) ;
  if( l < 1 ){ return 0 ; }
  double d, range = dmax - dmin ;
  nodePDBamino_acid *aa1Ptr = firstPtr ;
  nodePDBamino_acid *aa2Ptr = firstPtr->next->next ;/*avoid next neighbors*/
  for(int i=1; i<l-1; i++){
    for(int j=i+2; j<=l; j++){
      d=aa1Ptr->theAA.minD( aa2Ptr->theAA ) ;
      if( d<dmin ){ map[j-1][i-1]++ ; }
      else if( dmin<=d && d<=dmax ){  map[j-1][i-1] += (dmax-d)/range ; }
      aa2Ptr = aa2Ptr->next ;
    }
    aa1Ptr = aa1Ptr->next ;
    aa2Ptr = aa1Ptr->next->next ;
  }
}
/*================================================================*/
int PDBchain::createDistContMap( PDBchain &chain2, double **map,
				 double &dmin, double &dmax ){
  if( this->isEmpty( ) || chain2.isEmpty( ) ){ return 0 ; }
  int l = this->length( ) ;
  if( chain2.length( ) != l ){ return 0 ; }
  double d, range = dmax - dmin ;
  nodePDBamino_acid *aa1Ptr = firstPtr ;
  nodePDBamino_acid *aa2Ptr ;
  for(int i=0; i<l; i++){
    aa2Ptr = chain2.firstPtr ;
    for(int j=0; j<l; j++){
      d=aa1Ptr->theAA.minD( aa2Ptr->theAA ) ;
      if( d<dmin ){ map[j][i]++ ; }
      else if( dmin<=d && d<=dmax ){  map[j][i] += (dmax-d)/range ; }
      aa2Ptr = aa2Ptr->next ;
    }
    aa1Ptr = aa1Ptr->next ;
  }
}
/*================================================================*/
/*do not count as contact if |j-i|<=skip*/
/*contact definition given by the cont_sel type and the cutOff*/
/*store contact in symmetric **map[0..l-1][0..l-1]. We asssume **map has been previously allocated*/
int PDBchain::createSymmetricContactMap(int &skip, cont_sel type, double **map, double &cutOff){
  int l=this->length(); 
  if(l<skip+2) return 0;
  /*The following contact selection types are applicable only to an
    asymmetric contact map*/
  if(type==_F_HVBK_S_HVSC_||type==_F_HVSC_S_HVBK_) return 0;
  /*select appropriate pointer to doContact member function*/
  bool (PDBamino_acid::*doContact)(PDBamino_acid &, double &)=sel_doCont(type);
  
  PDBamino_acid *aa;
  nodePDBamino_acid *aa1Ptr = firstPtr;
  nodePDBamino_acid *aa2Ptr = firstPtr->next;
  int nc=0;/*number of contacts*/
  int i=0,j=1;
  while(aa1Ptr){
    aa=&(aa1Ptr->theAA);
    while(aa2Ptr && j-i<=skip){
      map[i][j]=map[j][i]=0.; /*printf("%3d %3d\n",i,j);*/
      aa2Ptr=aa2Ptr->next;/*avoid close neighbors*/
      j++; 
    }
    if(j-i<=skip) break;/*we reached end of chain*/
    while(aa2Ptr){
      if( (aa->*doContact)(aa2Ptr->theAA,cutOff)){
	map[i][j]=map[j][i]=1.; /*printf("%3d %3d\n",i,j);*/
	nc++; 
      }
      else{
	map[i][j]=map[j][i]=0.; /*printf("%3d %3d\n",i,j);*/
      }
      aa2Ptr=aa2Ptr->next; 
      j++; 
    }
    aa1Ptr=aa1Ptr->next;
    i++;
    if(!aa1Ptr) break;/*special case when skip==0*/
    aa2Ptr=aa1Ptr->next;
    j=i+1;
  }
  /*dump_array2(map,l,l); exit(0);*/
  return nc;
}
/*================================================================*/
/*do not count as contact if |j-i|<=skip
  store contact in symmetric **map[0..l-1][0..l-1].
  We asssume **map has been previously allocated
*/
int PDBchain::createSymmetricContactMap(int &skip, gen_atom_type type,
					double **map, double &cutOff){
  int l=this->length(),n; 
  if(l<skip+2) return 0;
  PDBamino_acid *aa;
  nodePDBamino_acid *aa1Ptr = firstPtr;
  nodePDBamino_acid *aa2Ptr = firstPtr->next;
  int nc=0;/*number of contacts*/
  int i=0,j=1;
  while(aa1Ptr){
    aa=&(aa1Ptr->theAA);
    while(aa2Ptr && j-i<=skip){
      map[i][j]=map[j][i]=0.; /*printf("%3d %3d\n",i,j);*/
      aa2Ptr=aa2Ptr->next;/*avoid close neighbors*/
      j++; 
    }
    if(j-i<=skip) break;/*we reached end of chain*/
    while(aa2Ptr){
      n=aa->doContact(aa2Ptr->theAA,type,type,cutOff);/*number of contacts*/
      if(n){
	map[i][j]=map[j][i]=n; /*printf("%3d %3d\n",i,j);*/
	nc++;/*increase the number of contacts between amino acids*/ 
      }
      else{
	map[i][j]=map[j][i]=0.; /*printf("%3d %3d\n",i,j);*/
      }
      aa2Ptr=aa2Ptr->next; 
      j++; 
    }
    aa1Ptr=aa1Ptr->next;
    i++;
    if(!aa1Ptr) break;/*special case when skip==0*/
    aa2Ptr=aa1Ptr->next;
    j=i+1;
  }
  /*dump_array2(map,l,l); exit(0);*/
  return nc;
}
/*================================================================*/
/*map[i][j] is different than map[j][i]*/
int PDBchain::createAsymmetricContactMap(int &skip, cont_sel type, 
					 double **map, double &cutOff){
  /*printf("createAsymmetricContactMap(..)\n");*/
  int l=this->length(); 
  if(l<skip+2) return 0;
  /*The following contact selection types are applicable only to a
    symmetric contact map*/
  if(type==_CA_CA_||type==_CB_CB_||type==_HV_HV_||type==_HVSC_HVSC_||
     type==_HVBK_HVBK_||type==_HVBK_HVSC_) return 0;
  bool (PDBamino_acid::*doContact)(PDBamino_acid &, double &)=sel_doCont(type);
  int i1=1,i2,nc=0;
  PDBamino_acid *aa;
  nodePDBamino_acid *aa1Ptr=firstPtr,*aa2Ptr=NULL;
  while(aa1Ptr){ /*printf("i1=%d\n",i1);if(aa1Ptr->next==NULL)cout<<"last\n";*/
    aa=&(aa1Ptr->theAA);
    aa2Ptr=firstPtr; i2=1;
    while(aa2Ptr){
      while(abs(i2-i1)<=skip && aa2Ptr){ /*skip [i1-skip,i1+skip] neighbors*/
	map[i1-1][i2-1]=0;
	i2++;
	aa2Ptr=aa2Ptr->next;
      }
      if(aa2Ptr){
	if( (aa->*doContact)(aa2Ptr->theAA,cutOff)){
	  map[i1-1][i2-1]=1.; /*printf("%3d %3d\n",i1,i2);*/
	  nc++; 
	}
	else map[i1-1][i2-1]=0.; /*printf("%3d %3d\n",i1,i2);*/
	aa2Ptr=aa2Ptr->next; i2++;
      }
    }
    aa1Ptr=aa1Ptr->next; i1++;
  }
  return nc;
}/*Matches int createAsymmetricContactMap(..)*/
/*================================================================*/
/*map[i][j] is different than map[j][i]*/
int PDBchain::createAsymmetricContactMap(int &skip, gen_atom_type type1,
					 gen_atom_type type2, double **map,
					 double &cutOff){
  /*printf("createAsymmetricContactMap(..)\n");*/
  /*cout<<"type1=\""<<gat2str(type1)<<"\" type2=\""<<gat2str(type2)<<"\"\n";*/
  int l=this->length(); 
  if(l<skip+2) return 0;
  int i1=1,i2,nc=0,n;
  PDBamino_acid *aa;
  nodePDBamino_acid *aa1Ptr=firstPtr,*aa2Ptr=NULL;
  while(aa1Ptr){ /*printf("i1=%d\n",i1);if(aa1Ptr->next==NULL)cout<<"last\n";*/
    aa=&(aa1Ptr->theAA);
    aa2Ptr=firstPtr; i2=1;
    while(aa2Ptr){
      while(abs(i2-i1)<=skip && aa2Ptr){ /*skip [i1-skip,i1+skip] neighbors*/
	map[i1-1][i2-1]=0;
	i2++;
	aa2Ptr=aa2Ptr->next;
      }
      if(aa2Ptr){
	n=aa->doContact(aa2Ptr->theAA,type1,type2,cutOff);
	if(n){	  /*cout<<*aa<<endl<<aa2Ptr->theAA<<"\n\n";*/
	  map[i1-1][i2-1]=n; /*printf("%3d %3d\n",i1,i2);*/
	  nc++; 
	}
	else map[i1-1][i2-1]=0.; /*printf("%3d %3d\n",i1,i2);*/
	aa2Ptr=aa2Ptr->next; i2++;
      }
    }
    aa1Ptr=aa1Ptr->next; i1++;
  }
  return nc;
}/*Matches int createAsymmetricContactMap(..)*/
/*================================================================*/
int PDBchain::createContactMap(int &skip, cont_sel type, double **map,
			       double &cutOff){
  if(symmetric_map(type))
    return this->createSymmetricContactMap(skip,type,map,cutOff);
  return this->createAsymmetricContactMap(skip,type,map,cutOff);
}
/*================================================================*/
int PDBchain::createContactMap(int &skip, gen_atom_type type1,
			       gen_atom_type type2, double **map,
			       double &cutOff){
  if(type1==type2)
    return this->createSymmetricContactMap(skip,type1,map,cutOff);
  /*cout<<"createAsymmetricContactMap\n";*/
  return this->createAsymmetricContactMap(skip,type1,type2,map,cutOff);
}
/*================================================================*/
/*don't care contact for amino acids i<=skip*/
/*definition of contact in type (CA,CB,HV)*/
/*contact cutoff in Angstroms*/
double PDBchain::co(int &skip, cont_sel type, double &cutoff){
  int l=this->length();
  double cont_ord=0.;
  if(l<skip+2) return 0;
  double **map=alloc_array(l,l);
  int nc=this->createContactMap(skip,type,map,cutoff);
  for(int i=0;i<l-skip;i++){
    for(int j=i+1+skip;j<l;j++){
      if(map[i][j]>0.)
	cont_ord+=j-i;	  
    }
  }
  return cont_ord/nc;
}
/*================================================================*/
int PDBchain::createCACBchain( listPDBatom &theCACBChain ){
  nodePDBamino_acid *currentPtr = firstPtr ;
  PDBatom *theCACB ;

  theCACBChain.empty() ;
  while( currentPtr ){
    theCACB = ( (*currentPtr).theAA ).getAtomFromName( " CA " ) ;
    theCACBChain.insertAtBack( *theCACB );
    delete theCACB ; 
    if( !( (*currentPtr).theAA ).is("GLY") ){
      theCACB = ( (*currentPtr).theAA ).getAtomFromName( " CB " ) ;
      theCACBChain.insertAtBack( *theCACB );
      delete theCACB ;
    }
    currentPtr = (*currentPtr).next ;   
  }    
  return 0;
}

void PDBchain::insert_H( ){
  if( ( *this ).isEmpty( ) ){  return ; }
  int inserted=0;
  //  nodePDBamino_acid *prevPtr=firstPtr,*currentPtr=prevPtr->next;//avoid N-terminal
  nodePDBamino_acid *prevPtr,*currentPtr;
  prevPtr=firstPtr;
  currentPtr=prevPtr->next;

  PDBatom *C, *CA, *N, H ;
  PDBvector Cv, CAv, Nv, Hv, NCv, NCAv, NHv ;
  PDBamino_acid *aaprev, *aacurr ;
  string const c=" C  ", n=" N  ", ca=" CA ", hn=" H  ";
  double a = -1.0296697692, b = -1.0166677430, bond_length=1.03 ;
  /*unit_vec(N->H) = a * unit_vec(N->C) + b * unit_vec(N->CA) */
  aaprev=&(prevPtr->theAA);
  while(currentPtr){
    aacurr = &(currentPtr->theAA) ;
    if(!aacurr->hasAtomWithName(hn)){
      C = aaprev->getAtomFromName( c ) ;
      Cv = C->getCoord( );
      N = aacurr->getAtomFromName( n ) ;
      Nv =N->getCoord( );
      CA = aacurr->getAtomFromName( ca ) ;
      CAv = CA->getCoord( );
      NCv = Cv ; NCv -= Nv ; NCv.normalize( ) ; NCv *= a ;
      NCAv = CAv ; NCAv -= Nv ; NCAv.normalize( ) ; NCAv *= b ;
      NHv = NCv ; NHv += NCAv ; NHv *= bond_length ;
      Hv = Nv ; Hv += NHv ;
      H = *N ; /*create atom, initialized to be like the nitrogen*/
      H.assignName((std::string&)(hn)) ;  /*change name to hydrogen*/
      H.assignCoord(Hv) ; /*change coordinates to hydrogen*/
      aacurr->addAtom( H ) ;
      delete C ;       delete CA ;       delete N ;
      inserted++;
    }
    prevPtr=currentPtr;
    aaprev=aacurr;
    currentPtr=prevPtr->next;
  }
  if(inserted){ int first=1; (*this).renumberFullyAtomSerialNumber(first) ; }
}
/*================================================================*/
void PDBchain::remove_all_hydrogens( ){
  nodePDBamino_acid *currentPtr = this->firstPtr ;
  while(currentPtr){
    currentPtr->theAA.remove_all_hydrogens( ) ;
    currentPtr = currentPtr->next ;
  }
}
/*================================================================*/
void PDBchain::translate( const PDBvector &shift ){
  if( ( *this ).isEmpty( ) ){  return ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    currentPtr->theAA.translate( shift ) ;
    currentPtr = currentPtr->next ;   
  }
}
/*================================================================*/
void PDBchain::rotate( const PDBvector &axis, const double &angle ){
  if( ( *this ).isEmpty( ) ){  return ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    currentPtr->theAA.rotate( axis, angle) ;
    currentPtr = currentPtr->next ;   
  }
}
/*================================================================*/
void PDBchain::rotate( const PDBvector &new_origin, const PDBvector &axis,
		       const double &angle ){
  if( ( *this ).isEmpty( ) ){  return ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    currentPtr->theAA.rotate( new_origin, axis, angle) ;
    currentPtr = currentPtr->next ;   
  }
}
/*================================================================*/
void PDBchain::rotate( double ** rotation_matrix ){
  if( ( *this ).isEmpty( ) ){  return ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    currentPtr->theAA.rotate( rotation_matrix ) ;
    currentPtr = currentPtr->next ;   
  }
}
/*================================================================*/
void PDBchain::rotatePhiBy( const double &phi, const int &index ){
  /*cout<<"void PDBchain::rotatePhiBy( const double &, const int & )\n";*/
  if( (phi==0) || this->isEmpty( ) ){ return ; }
  int l = this->length( ) ;
  if( index > l ){ return ; }
  PDBamino_acid *center_aa = this->getAminoFromIndex( index ) ;
  PDBatom *CA = center_aa->getAtomFromName( " CA " ) ;
  PDBatom *N = center_aa->getAtomFromName( " N  " ) ;
  PDBvector origin, axis ;
  /*Looking from the CA atom to the N atom, a positive rotation is 
    counter-clock wise*/
  origin = CA->getCoord( ) ;
  axis = N->getCoord( ) - origin ;
  /*rotate central amino and rest of chain after central amino acid */
  int i=1;
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    if( i>= index ){
      currentPtr->theAA.rotate( origin, axis, phi ) ;
    }
    currentPtr = currentPtr->next ;  i++ ;
  }
}
/*================================================================*/
void PDBchain::rotatePsiBy( const double &psi, const int &index ){
  /*cout<<"void PDBchain::rotatePsiBy( const double &, const int & )\n";*/
  if( (psi==0) || this->isEmpty( ) ){ return ; }
  int l = this->length( ) ;
  if( index > l ){ return ; }
  PDBamino_acid *center_aa = this->getAminoFromIndex( index ) ;
  PDBatom *CA = center_aa->getAtomFromName( " CA " ) ;
  PDBatom *C = center_aa->getAtomFromName( " C  " ) ;
  PDBatom *O = center_aa->getAtomFromName( " O  " ) ;
  PDBvector origin, axis ;
  /*Looking from the C' atom to the CA atom, a positive rotation is 
    counter-clock wise*/
  origin = C->getCoord( ) ;
  axis = CA->getCoord( ) - origin ;
  /*we only have to rotate atom " O  " of the central amino acid*/
  O->rotate( origin, axis, psi ) ;
  double *new_xyz = ( O->getCoord( ) ).getComponents( ) ;
  this->changeCoordOfAtomNameAtIndex( " O  ", index, new_xyz ) ;
  /*now rotate rest of chain after the central amino acid*/
  int i=1;
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    if( i > index ){
      currentPtr->theAA.rotate( origin, axis, psi ) ;
    }
    currentPtr = currentPtr->next ;  i++ ;
  } 
}
/*================================================================
Phi dihedral angle between planes defined by C'-N-CA and N-CA-C''
  cos(Phi) = A * B, where A unit vector parallel to C'N x NA and
                    B unit vector parallel to NCA x CAC''
         C''
        /
  N - CA    <---- In this extended orientation, Phi = -180
 /     
C'                                                          */
double PDBchain::getPhiAtIndex( const int &index ){
  if( this->length( ) < 2 ){
    cout << "chain too short to output Phi at position "<<index<<endl;
    exit(1) ;
  }
  if(index == 1){ return 0.0 ; }
  double Phi ;
  PDBamino_acid *prev = this->getAminoFromIndex( index-1 ) ;
  PDBamino_acid *curr = this->getAminoFromIndex( index ) ;
  PDBatom *prevC  = prev->getAtomFromName( " C  " ) ;
  PDBatom *currN  = curr->getAtomFromName( " N  " ) ;
  PDBatom *currCA = curr->getAtomFromName( " CA " ) ;
  PDBatom *currC  = curr->getAtomFromName( " C  " ) ;
  PDBvector CN  = currN->getCoord( )  -  prevC->getCoord( )  ;
  PDBvector NCA = currCA->getCoord( ) -  currN->getCoord( )  ;
  PDBvector CAC = currC->getCoord( )  -  currCA->getCoord( ) ;
  PDBvector A = CN ^ NCA ;
  PDBvector B = NCA ^ CAC ;
  PDBvector Z = NCA ^ A ; /*provides sign convention*/
  A.normalize( ) ; B.normalize( ) ;
  Phi = 180*acos( A*B )/PI ;
  if( Z*B < 0 ){ Phi *= -1 ; }
  return Phi ;
}
/*================================================================
Psi dihedral angle between the two planes defined by N-CA-C' and CA-C'-N
  cos(Psi) = A * B, where A unit vector parallel to NCA x CAC and
                    B unit vector parallel to CAC x C'N
         N
        /
  CA - C'    <---- In this extended orientation, Phi = -180
 /     
N                                                          */
double PDBchain::getPsiAtIndex( const int &index ){
  if( this->length( ) < 2 ){ cout << "chain too short!\n"; exit(1) ;}
  if(index == this->length( ) ){ return 0.0 ; }
  double Psi ;
  PDBamino_acid *curr = this->getAminoFromIndex( index ) ;
  PDBamino_acid *next = this->getAminoFromIndex( index+1 ) ;
  PDBatom *currN  = curr->getAtomFromName( " N  " ) ;
  PDBatom *currCA = curr->getAtomFromName( " CA " ) ;
  PDBatom *currC  = curr->getAtomFromName( " C  " ) ;
  PDBatom *nextN  = next->getAtomFromName( " N  " ) ;
  PDBvector NCA = currCA->getCoord( ) - currN->getCoord( )  ;
  PDBvector CAC = currC->getCoord( )  - currCA->getCoord( ) ;
  PDBvector CN =  nextN->getCoord( )  - currC->getCoord( )  ;
  PDBvector A = NCA ^ CAC ;
  PDBvector B = CAC ^ CN ;
  PDBvector Z = CAC ^ A ; /*provides sign convention*/
  A.normalize( ) ; B.normalize( ) ;
  Psi = 180*acos( A*B )/PI ;
  if( Z*B < 0 ){ Psi *= -1 ; }
  return Psi ;
}
/*================================================================*/
double PDBchain::drmsPhiPsi( PDBchain &chain2 ){
  double l = this->length( ), drms=0.0, d ;
  if( chain2.length( ) != l ){ return -1 ; }
  if( l < 2 ){ cout << "chain too short!\n"; exit(1) ; }
  for( int i=1; i<=l; i++ ){
    d = this->getPhiAtIndex( i ) - chain2.getPhiAtIndex( i ) ; drms += d*d ;
    /*cout<<"Phi1-Phi2="<<d;*/
    d = this->getPsiAtIndex( i ) - chain2.getPsiAtIndex( i ) ; drms += d*d ;
    /*cout<<" Psi1-Psi2="<<d<<endl;*/
  }
  /*cout<<"drmsPhiPsi="<<sqrt( drms/(2*l))<<endl;*/
  return sqrt( drms/(l-1) ) ;
}
/*================================================================*/
PDBvector PDBchain::get_CA_CM( ){
  listPDBatom ca_chain;
  this->createCAchain( ca_chain ) ; 
  return ca_chain.get_CM( ) ;
}
/*================================================================*/
PDBvector PDBchain::get_GC( ){
  listPDBatom single_list ;
  this->fuseToSingleList( single_list ) ;
  return single_list.get_CM( ) ;
}
/*================================================================*/
double PDBchain::get_CA_surrounding_sphere( ){
  listPDBatom ca_chain;
  this->createCAchain( ca_chain ) ; 
  return ca_chain.get_surrounding_sphere( ) ;
}
/*================================================================*/
void PDBchain::filterOnly( const string *names, const int &n ){
  if( this->isEmpty( ) || n==0 ){  return ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    currentPtr->theAA.filterOnly( names, n ) ;
    currentPtr = currentPtr->next ;   
  } 
}
/*================================================================*/
bool  PDBchain::is_there_resName_name_at_resIndexII( const string &aaname_and_atname, const int &resIndex                                 ){
  if( resIndex<1 || resIndex>this->length( ) ){ return false ; }
  string aaname = aaname_and_atname.substr( 0, 3 ) ;
  string atname = aaname_and_atname.substr( 4, 4 ) ;
  /*cout<<"aaname=|"<<aaname<<"|, atname=|"<<atname<<"|, rexIndex="<<resIndex<<endl;*/
  PDBamino_acid *current_aa = this->getAminoFromIndex( resIndex ) ;
  if( aaname=="BKB" ){/*atom of the backbone*/
    /*exceptions: PRO-N  , first N and last O do not hydrogen bond*/
    if( (current_aa->getName( )=="PRO") && (atname==" N  ") ){ return false ; }
    if( (resIndex==1) && (atname==" N  ") ){ return false ; }
    if( (resIndex==this->length( )) && (atname==" O  ") ){ return false ; }
    /*now check presence of atom*/
    if( current_aa->hasAtomWithName( atname ) ){ return true ; }
  }    
  else{/*atom of the sidechain*/
    if( aaname != current_aa->getName( ) ){ return false ; }
    if( current_aa->hasAtomWithName( atname ) ){ return true ; }
  }
  return false ;
}
/*================================================================*/
bool PDBchain::isAtomNameAtAminoAcidIndex(const string &name,
					       const int &idx){
  if(this->length()<idx) return false;
  nodePDBamino_acid *currentPtr=firstPtr;
  int i=1;
  while(i<idx){ currentPtr=currentPtr->next; i++; }
  if(currentPtr->theAA.hasAtomWithName(name)) return true;
  return false;
}
/*================================================================*/
string PDBchain::output_one_letter_sequence_from_indexes( const int &a,
							  const int &b ){
  string tmp, name, nm ;
  for( int i=a; i<=b; i++ ){
    name = this->getNameFromIndex( i ) ;
    nm =  three_letter_code_to_one_letter_code( name ) ;
    tmp += nm ;
  }
  return tmp ;
}
/*================================================================*/
string PDBchain::output_one_letter_sequence( ){
  int a=1, b=this->length( ) ;
  return this->output_one_letter_sequence_from_indexes( a, b ) ;
}
/*================================================================*/
void PDBchain::output_three_letter_sequence_from_indexes( const int &a,
					    const int &b, ostream &OUT ){
  for( int i=a; i<=b; i++ ){ OUT<< this->getNameFromIndex( i )<<endl ; }
}
/*================================================================*/
void PDBchain::output_three_letter_sequence( ostream &OUT ){
  int a=1, b=this->length( ) ;
  this->output_three_letter_sequence_from_indexes( a, b, OUT ) ;
}
/*================================================================*/
void PDBchain::fuseToSingleList( listPDBatom &theList ){
  listPDBatom aa_list ;
  if( !theList.isEmpty( ) ){ theList.empty( ) ; }
  nodePDBamino_acid *currentPtr = firstPtr ;
  while(currentPtr){
    currentPtr->theAA.fuseToSingleList( aa_list ) ;
    theList.insertAtBack( aa_list ) ;
    currentPtr = currentPtr->next ;   
  } 
}
/*================================================================*/
/*n is lenght of names array*/
double PDBchain::alignByRmsd( PDBchain &native, const string *names, 
			      const int &n ){
  double rmsd, **R=alloc_array( 3, 3 ) ; /*R rotation matrix*/
  if( this == &native ){ return 0.0 ; }/*self-assignment*/
  /*put all atoms of both PDBchains into two listPDBatom lists*/
  listPDBatom natL, confL ;
  native.fuseToSingleList( natL ) ;  this->fuseToSingleList( confL ) ;
  /*leave only the selected atoms*/
  if( names[0] != " ALL" ){
    natL.filterOnly( names, n ) ;  confL.filterOnly( names, n ) ;
  }
  /*superimpose geometric centers on the origin*/
  PDBvector natCM, confCM ;
  natCM  =  natL.get_CM( ) ;   natL.translate( natCM  * (-1) ) ;
  confCM = confL.get_CM( ) ;  confL.translate( confCM * (-1) ) ;
  /*obtain the rotation matrix*/
  rmsd = confL.rot_matrix_to_adapt_by_rmsd_to( natL, R ) ;
  /*cout<<rmsd<<"  "<<rmsd*det_3x3( R )<<endl;*/
  /*Align now the whole protein*/
  this->translate( confCM * (-1) ) ;
  this->rotate( R ) ;
  this->translate( natCM ) ;

  return rmsd ;
}
/*================================================================*/
/*One contact between any two amino acids is enough for the chains to be
  in contact*/
bool PDBchain::do_contact(PDBchain *chain2, bool (*do_aa_contact)(PDBamino_acid*,PDBamino_acid*)){
  nodePDBamino_acid *pt=firstPtr, *pt2;
  PDBamino_acid *aa1,*aa2;
  while(pt){
    pt2=chain2->firstPtr;
    aa1=&(pt->theAA);
    while(pt2){
      aa2=&(pt->theAA);
      if(do_aa_contact(aa1,aa2)){ return true; }
      pt2=pt2->next;
    }
    pt=pt->next;
  }
  return false;
}
/*=====================================================*/
double PDBchain::average_graph_length(int &skip, cont_sel type, double &cutOff){
  int l=this->length();
  if(l<skip+2) return 0.;
  double **map=alloc_array(l,l);
  int ck=0;/*no skip when creating the contact map*/
  this->createContactMap(ck,type,map,cutOff);/* dump_array2(map,l,l);*/
  double **lengths=alloc_array(l,l);/*lenthg of the path from i to j*/
  init_array(0.,lengths,l);
  double **B=alloc_array(l,l);/*stores powers of map array*/
  assign_array(B,map,l);
  for(int k=1;k<l;k++){
    for(int i=0;i<l;i++)
      for(int j=i+1+ck;j<l;j++)
	if(B[i][j]>0. && lengths[i][j]<=0.){
	  lengths[i][j]=k;/*number of minimum links from i to j*/
	}
    mult_arrayBA(B,map,l);/*increase the order of B by one*/
  }
  int nc=0;
  double agl=0.;
  for(int i=0;i<l;i++)
    for(int j=i+1+skip;j<l;j++){
      nc++;
      agl+=lengths[i][j];
    }
  return agl/nc;
}
/*=====================================================*/
int  PDBchain::number_of_contacts(int &skip, cont_sel type, double &cutOff){
  int l=this->length();
  if(l<skip+2) return 0;
  double **map=alloc_array(l,l);
  return this->createContactMap(skip,type,map,cutOff);
}
/*=====================================================*/
void PDBchain::output_fasta_seq( ostream &out, string &header ){
  out<<">"<<header<<endl;
  string all=this->output_one_letter_sequence( ) ;
  int l=this->length(),x=80;
  for(int i=0;i<l;i+=x) out<<all.substr(i,x)<<endl;
}
/*=====================================================*/

/*------------ class nodePDBchain  --------------*/

nodePDBchain::nodePDBchain( PDBchain &chain, const string &theChainID){
  theChain = chain ; chainID = theChainID ; next = NULL ;
}
/*================================================================*/



/*------------ class PDBchains ------------------------*/

ostream &operator<<(ostream &output, const PDBchains &theList ){
  nodePDBchain *currentPtr = theList.firstPtr ;
  while( currentPtr != NULL){
    output << currentPtr->theChain ;
    currentPtr = currentPtr->next ;
    if(currentPtr){ output << "TER\n"; }
    else{ output << "END\n"; }
  }
  return output ;
}
/*==============================================================*/
ifstream &operator>>(ifstream &input, PDBchains &the_chains ){
  /*cout<<"ifstream &operator>>(ifstream &,PDBchains &)\n";*/
  the_chains.importChains( input ) ;
  return input ;
}
/*==============================================================*/
nodePDBchain *PDBchains::getNewNode(PDBchain &chain,const string &theChainID){
  nodePDBchain *ptr = new nodePDBchain( chain, theChainID );
  assert( ptr != 0 ) ;
  return ptr ; 
}
/*================================================================*/
PDBchains::PDBchains( ){ firstPtr = lastPtr = NULL ; }
/*================================================================*/
PDBchains::PDBchains(char *file){
  char *uzfile=new char[512];
  /*in case file is zipped*/
  zip_format fmt=returnZipFormat(file);  unzip(uzfile,file);
  /*cout<<"uzfile=|"<<uzfile<<"|\n";*/
  ifstream IN(uzfile);
  firstPtr = lastPtr = NULL ;
  this->importChains(IN);
  IN.close();  
  zip(uzfile,fmt);
}
/*================================================================*/
PDBchains::~PDBchains( ){ this->empty( ) ; }
/*================================================================*/
void PDBchains::empty( ){
 if(!this->isEmpty( ))
   while(firstPtr != NULL)
     this->rmFromFront( ); 
}
/*================================================================*/
bool PDBchains::isEmpty( ){ return firstPtr == NULL ; }
/*================================================================*/
void PDBchains::rmFromFront( ){
  nodePDBchain *currentPtr = firstPtr ;
  firstPtr = currentPtr->next ; delete currentPtr ;
}
/*================================================================*/
void PDBchains::insertAtBack( PDBchain &the_chain , const string &thechainID ){
  nodePDBchain *newPtr = getNewNode( the_chain , thechainID ) ;
  if( this->isEmpty() ){ firstPtr = lastPtr = newPtr ; }
  else{ lastPtr->next = newPtr ; lastPtr = newPtr ; }
}
/*================================================================*/
int PDBchains::importChains( ifstream &PDB ){
  /*cout<<"int PDBchains::importChains(ifstream &)\n";*/
  string line,chainID ;
  int position;
  PDBchain the_chain ;

  this->empty();
  if( PDB.eof() ){ return 1 ; }
  /*while(getline(PDB,line)) cout<<line;*/
  do{
    PDB >> the_chain ; /*cout<<the_chain<<endl;*/
    if(!the_chain.isEmpty()){
      chainID=the_chain.printChainID();
      /*cout<<"chainID=|"<<chainID<<"|\n";exit(0);*/
      this->insertAtBack( the_chain, chainID ) ;
    }
    position=PDB.tellg(); getline(PDB,line); /*cout<<line<<endl;*/
    if(is_ATOM(line)) PDB.seekg(position);
  }while( !is_END(line) && !PDB.eof()) ;
  /*cout<<"number of chains="<<this->length()<<endl;*/
  return 0 ;
}
/*=====================================================*/
int PDBchains::length( ){
  int count = 0; nodePDBchain  *currentPtr = firstPtr ;
  while( currentPtr != NULL ){ 
    count++ ;  currentPtr = (*currentPtr).next ;
  }
  return count ;

}
/*=====================================================*/
PDBchain* PDBchains::pointToPDBchainFromIndex( const int &index ){
  if( (index > this->length( )) || (index<1) ){ return NULL; }
  nodePDBchain *currentPtr = firstPtr ; 
  int count = 1 ;
  while( count < index ){ currentPtr = currentPtr->next ;  count++ ; }
  return &(currentPtr->theChain) ;
}
/*=====================================================*/
PDBchain* PDBchains::getPDBchainFromIndex( const int &index ){
  PDBchain *tmp=NULL;
  if( (index > this->length( )) || (index<1) ){ return tmp; }
  nodePDBchain *currentPtr = firstPtr ; 
  int count = 1 ;
  while( count < index ){ currentPtr = currentPtr->next ;  count++ ; }
  tmp=new PDBchain; *tmp = currentPtr->theChain ;
  return tmp ;
}
/*=====================================================*/ 
bool PDBchains::getPDBchainFromIndex(PDBchain &out, const int &idx ){
  /*cout<<"idx="<<idx<<" l="<<this->length()<<endl;*/
  if(idx>this->length() || idx<1){return false;}
  nodePDBchain *currentPtr = firstPtr ; 
  int count = 1 ;
  while( count < idx ){ currentPtr = currentPtr->next ;  count++ ; }
  out=currentPtr->theChain; /*cout<<out;*/
  return !out.isEmpty();
}
/*=====================================================*/
PDBchain* PDBchains::pointToPDBchainFromChainID(string &ID){
  nodePDBchain *currentPtr = firstPtr ;
  if(ID=="_"){
    return this->pointToPDBchainFromIndex(1);
  }
  while(currentPtr){
    if(currentPtr->chainID==ID) return &(currentPtr->theChain);
    currentPtr=currentPtr->next;
  }
  return NULL;
}
/*=====================================================*/
 bool PDBchains::getPDBchainFromChainID(PDBchain &out, const string &ID ){
   nodePDBchain *currentPtr = firstPtr ; 
   /*cout<<"l="<<this->length()<<endl; cout<<"ID="<<ID<<endl;*/
   while(currentPtr){
     /*cout<<"current ID="<<currentPtr->chainID<<endl;*/
     if(currentPtr->chainID==ID){
       out=currentPtr->theChain; /*cout<<out;exit(0);*/
       return !out.isEmpty();
     }
     currentPtr=currentPtr->next;
   }
   return false;  
 }
/*=====================================================*/

/*
int main(void){



  return 0 ;
}
*/
