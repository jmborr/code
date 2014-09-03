#include <stdio.h>       
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SEG 15
#define NLAYER1 330
#define NLAYER2 150
#define NLAYER3 2
#define MLAYER1 60
#define MLAYER2 45
#define MLAYER3 2
#define MAX_LEN 2000
#define LenScale 1000

void Test(char sprffile[20]);
int Readmat(char file[20]);
void FeedNet(int I);
void Netout();
void FeedNet1(int I);
void Netout1();
void result1(char *file);
int Max(float *p,int N);

char StrA1[MAX_LEN],StrA2[MAX_LEN],prdA2[MAX_LEN];
float Score[MAX_LEN][21];
float nout[MAX_LEN][NLAYER3],pout[MAX_LEN][NLAYER3];
int Len;
float w1[NLAYER1][NLAYER2],Dw1[NLAYER1][NLAYER2];
float w2[NLAYER2][NLAYER3],Dw2[NLAYER2][NLAYER3];
float in1[NLAYER1],in2[NLAYER2],out2[NLAYER2],in3[NLAYER3],out[NLAYER3];
float tgtout[NLAYER3];
float v1[MLAYER1][MLAYER2],Dv1[MLAYER1][MLAYER2];
float v2[MLAYER2][MLAYER3],Dv2[MLAYER2][MLAYER3];
float im1[MLAYER1],im2[MLAYER2], mout2[MLAYER2],im3[MLAYER3],mout[MLAYER3];

float dw1[NLAYER1][NLAYER2],dw2[NLAYER2][NLAYER3];
float dv1[MLAYER1][MLAYER2],dv2[MLAYER2][MLAYER3]; // Derivatives
 
int n[NLAYER3][NLAYER3], m[MLAYER3][MLAYER3];
int Ntotal, Ngood, Mgood;
int Mmax=0;
 
int main(int argu, char *argv[]){
  int II,i,j,k,count;

  char wghtfile[30];
  FILE *fpout;
   

  if(argu<3) {
    fprintf(stderr,"USAGE: solve wgtfile pdb\n");
    return(-1);
  }

   
   strcpy(wghtfile, argv[1]);
    if((fpout=fopen(wghtfile,"r"))==NULL) {
	    printf("cannot open wghtfile %s\n", wghtfile);
	    exit(-1);
    }
    
    for(i=0;i<NLAYER1;i++){
      for(j=0;j<NLAYER2;j++){
        fscanf(fpout,"%f\n",&w1[i][j]);
        Dw1[i][j]=0.0;
      }
    }
    for(i=0;i<NLAYER2;i++){
      for(j=0;j<NLAYER3;j++){
        fscanf(fpout,"%f\n",&w2[i][j]);
        Dw2[i][j]=0.0;
      }
    }
    for(i=0;i<MLAYER1;i++){
      for(j=0;j<MLAYER2;j++){
        fscanf(fpout,"%f\n",&v1[i][j]);
        Dv1[i][j]=0.0;
      }
    }
    for(i=0;i<MLAYER2;i++){
      for(j=0;j<MLAYER3;j++){
        fscanf(fpout,"%f\n",&v2[i][j]);
        Dv2[i][j]=0.0;
      }
    }
    fclose(fpout);

    Test(argv[2]);
    
    return(0);

}

/******************************************************/
void Test(char sprffile[20])
{ int i,j;

  
	Readmat(sprffile);
	
      for(i=1;i<=Len;i++){
      FeedNet(i);
      Netout(); 
      for(j=0;j<NLAYER3;j++)
	nout[i][j]=out[j];
    }

    for(i=1;i<=Len;i++)  {
      FeedNet1(i);
      Netout1();
      for(j=0;j<MLAYER3;j++)
         pout[i][j]=mout[j];
      } 
      result1(sprffile); 

     Ntotal+=Len;

 } 
  
/*************************************************************/
int Readmat(char file[20]){
     int i,j,k;
     FILE *fp;
     int tmp[20];
     char ss[1000], filenm[45]="./";
     unsigned char c;

     strcat(filenm, file);
     strcat(filenm,".mat3");

     if((fp=fopen(filenm,"r"))==NULL){
           fprintf(stderr,"cannot open %s\n",filenm);
	   exit(-1);
      }

     

     for(i=0;i<3;i++){
       j=fgetc(fp);
       while(j!='\n'){
	 if(j==EOF){
	   fprintf(stderr,"ERROR from solve.c while reading %s\n",filenm);
	   exit(1); /*error, file is empty*/
	 }
	 j=fgetc(fp);
       }
     }

     Len=0;
     while(1){
	   fgets(ss,1000,fp);
	   if(strlen(ss)<50) break;
	    Len++;
	   sscanf(ss+6,"%c",&StrA1[Len]);
	         for(i=0;i<20;i++){
	        sscanf(ss+10+3*i,"%3f",&Score[Len][i]);
	      Score[Len][i]=1.0/(1.0+exp(-1.0*Score[Len][i]));
	       }
    }

     fclose(fp);

     return(0);
 }

/****************************************************************/
void FeedNet(int I){
  int i,j,k;
  
  k=0;
  for(i=(I-SEG/2);i<(I-SEG/2)+SEG;i++){
    if(i<1 | i>Len){
      for(j=0;j<20;j++){
	in1[k]=0.0;
	k++;
      }
      in1[k]=1.0;
      k++;
    in1[k]=(float)Len/LenScale;
    k++;
    }
    else{
      for(j=0;j<20;j++){
	in1[k]=Score[i][j];
	k++;
      }
      in1[k]=0.0;
      k++;
    in1[k]=(float)Len/LenScale;
    k++;
   }   
  }
  for(i=0;i<NLAYER3;i++)
    tgtout[i]=0.0;
  tgtout[(StrA2[I]-48)%2]=1.0; 

}


/***************************************************************/
void Netout(){
  int i,j;

  for(i=0;i<NLAYER2;i++){
    in2[i]=0.0;
    for(j=0;j<NLAYER1;j++){
      in2[i]+=w1[j][i]*in1[j];
    }
    out2[i]=1.0/(1.0+exp(-1.0*in2[i]));
  }

  for(i=0;i<NLAYER3;i++){
    in3[i]=0.0;
    for(j=0;j<NLAYER2;j++){
      in3[i]+=w2[j][i]*out2[j];
    }
    out[i]=1.0/(1.0+exp(-1.0*in3[i]));
  }
}

/***************************************************************/
void FeedNet1(int I){
  int i,j,k;
  
  k=0;
  for(i=(I-SEG/2);i<(I-SEG/2)+SEG;i++){
    if(i<1 | i>Len){
      for(j=0;j<NLAYER3;j++){
	im1[k]=0.0;
	k++;
      }
      im1[k]=1.0;
      k++;
      im1[k]=(float)Len/LenScale;
      k++;

    }
    else{
      for(j=0;j<NLAYER3;j++){
	im1[k]=nout[i][j];
	k++;
      }
      im1[k]=0.0;
      k++;
      im1[k]=(float)Len/LenScale;
      k++;
    }
  } 
  for(i=0;i<MLAYER3;i++)
    tgtout[i]=0.0;
  
  tgtout[(StrA2[I]-48)%2]=1.0; 
}

/***************************************************************/
void Netout1(){
  int i,j;

  for(i=0;i<MLAYER2;i++){
    im2[i]=0.0;
    for(j=0;j<MLAYER1;j++){
      im2[i]+=v1[j][i]*im1[j];
    }
    mout2[i]=1.0/(1.0+exp(-1.0*im2[i]));
  }

  for(i=0;i<MLAYER3;i++){
    im3[i]=0.0;
    for(j=0;j<MLAYER2;j++){
      im3[i]+=v2[j][i]*mout2[j];
    }
    mout[i]=1.0/(1.0+exp(-1.0*im3[i]));
  }
}

/***************************************************************/

void result1(char *file){
  int i,j;
  FILE *fp;
  char filenm[25]="./", suffix[6]=".neu";
  char true_exp;

  strcat(filenm,file);
  strcat(filenm,suffix);
  if((fp=fopen(filenm,"w"))==NULL) {
            printf("cannot open file %s\n",filenm);
             exit(-1);
         }
  fprintf(fp,"%5d\n",Len);

  j=0; 
  for(i=1;i<=Len;i++){

    prdA2[i]=Max(pout[i],MLAYER3)+48;

   fprintf(fp,"%5d %c  %5.2f  %5.2f  %c\n",i,StrA1[i],pout[i][0],
                               pout[i][1],prdA2[i]);
 
  }

  fclose(fp);
}

/***************************************************************/

int Max(float *p, int N){
  int i,n;
  float tmp=-9999.9;

  for(i=0;i<N;i++)
    if(*(p+i)>tmp){
      tmp=*(p+i);
      n=i;
    }
  return(n);
}





















