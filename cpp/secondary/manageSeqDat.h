#include<iostream>

#define INF_MSD 100000 /*veeeeery long sequence*/
#define alphaCO 4  /*minimum helix size*/
#define betaCO  3  /*minimum beta size*/

#define alphaCode 2
#define betaCode 4

double alpha_score=1.0/alphaCO+0.0001; /*extra 0.0001 to avoid rounding errors*/
double beta_score =1.0/betaCO +0.0001;

/*From a seq.dat-like file, return where do the secondary structure
  segments begin (either alpha or beta) and the length of each one. We
  need beginSec and resperSec be already allocated and big enough to
  contain all infor for the nSec secondary structure elements. If we
  are interested in a chunk of seq.dat, instead of all the file, then
  we pass a segment [nter,cter]. Residue indexes stored in beginSec follow
  the Fortran convention, i.e., begin with 1, not 0*/
int initSecSeg(char *seqdatfile, int *beginSec, int *resperSec,int nter=0, int cter=INF_MSD){
  FILE *in=fopen(seqdatfile,"r");
  int i,k;
  int confidence;
  double score=0.;
  char *aa=new char[4];
  int prev;       /*previous secondary structure code*/
  int curr=1;     /*secondary structure code of current residue*/
  int nSec=0;
  int l=0;
  int ls=0; /*number of residues in secondary structure segments*/
  double minp=1.0/3;
  int chunk=10;

  for(k=0;k<nter;k++) fscanf(in,"%d %s %d %d",&i,aa,&curr,&confidence);/*skip residues below nter*/
  k=0;/*will be number of residues in between nter and cter*/
  do{
    prev=curr;
    if(k==cter+1) break; /*skip residues above cter*/
    if(fscanf(in,"%d %s %d %d",&i,aa,&curr,&confidence)==EOF) break; /*will count number of lines*/
    k++;
    /*update score*/
    if(curr==prev){
      if(curr==alphaCode){ score+=alpha_score; }
      else if(curr==betaCode){ score+=beta_score; }
      else{ score=0.; }
    }
    else{
      if(curr==alphaCode){ score=alpha_score; }
      else if(curr==betaCode){ score=beta_score;}
      else{ score=0.;}
    }
    /*did we find a segment big enough?*/
    if(score>1){
      if(curr==alphaCode){ l=alphaCO; beginSec[nSec]=1+i-alphaCO-nter;}
      else{ l=betaCO;  beginSec[nSec]=1+i-betaCO-nter; }
      do{
	prev=curr;
	if(k==cter+1) break; /*skip residues above cter*/
	if(fscanf(in,"%d %s %d %d",&i,aa,&curr,&confidence)==EOF) break;
	k++;
	if(curr==prev){ l++;}
	else{
	  if     (curr==alphaCode) score=alpha_score;
	  else if(curr==betaCode)  score=beta_score;
	  else                  score=0.;	  
	}
      }while(score>1);
      resperSec[nSec]=l;
      ls+=l;
      l=0;
      nSec++;
    }
  }while(1);
  fclose(in);
  /*if residues in secondary structure elements account for less than
    minp, or if there's only two elements of secondary structure, then
    just divide the sequence in chunks of 11-residue segments, and use
    them secondary structure fragments*/
  /*printf("ls=%d k=%d ratio=%lf\n",ls,k,(1.0*ls)/k);exit(1);*/
  if( (1.0*ls)/k<minp || nSec<=2 ){
    nSec=0;
    i=0;
    while(k>chunk){
      resperSec[nSec]=chunk;
      beginSec[nSec++]=1+i+nter;
      i+=chunk;
      k-=chunk;
    }
  }
  /*printf("nSec=%d\n",nSec);*/
  return nSec;
}
