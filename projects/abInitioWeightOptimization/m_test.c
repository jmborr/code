
/* Simple optimization example. See README for details. */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define N_WEIGHTS 10
#define MAX_PROTEINS 1000
#define MAX_DECOYS 1500

int N_DECOYS[MAX_PROTEINS];

//#define N_DECOYS 100

int N_PROTEINS;

// data:
// 0 - tmscore
// 1.. - energy

double ***data;
int *native;

// 0 - eta
// 1 - bk
// 2.. - weights

double *parameters;

char filelist[1000][100];
int nfiles;

float initvalues[1000][2];

int read_data(char *lname)
{
  FILE *inp;
  char buf[1000];
	int i, j, k, nfiles, num, nat;
	double min;
		  
    nfiles = 0;
    inp = fopen(lname,"r");
    while (!feof(inp)) {
      if (fgets(buf,1000,inp)==buf) {
        sscanf(buf,"%s",filelist[nfiles++]);
      }
    }
    fclose(inp);
    
    printf("Reading %d files...\n", nfiles);

		native = (int*)calloc(sizeof(int)*nfiles,1);
		
		data = (double***)calloc(sizeof(double**)*nfiles,1);
		for (i=0;i<nfiles;i++) {
			data[i] = (double**)calloc(sizeof(double*)*MAX_DECOYS,1);
			for (j=0;j<MAX_DECOYS;j++) {
			  data[i][j] = (double*)calloc(sizeof(double)*15,1);
			}
		}

		N_PROTEINS = 0;
		    
		printf("\n");

    for (i=0;i<nfiles;i++) {
      inp = fopen(filelist[i],"r");
      if (inp) {
        N_DECOYS[N_PROTEINS] = 0;
        j = 0;
      	do {
  		    if (fgets(buf,1000,inp)==buf) {
    		    num = sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
    		    	&data[N_PROTEINS][j][0],
    		    	&data[N_PROTEINS][j][1],
    		    	&data[N_PROTEINS][j][2],
    		    	&data[N_PROTEINS][j][3],
    		    	&data[N_PROTEINS][j][4],
    		    	&data[N_PROTEINS][j][5],
    		    	&data[N_PROTEINS][j][6],
    		    	&data[N_PROTEINS][j][7],
    		    	&data[N_PROTEINS][j][8],
    		    	&data[N_PROTEINS][j][9],
    		    	&data[N_PROTEINS][j][10],
    		    	&data[N_PROTEINS][j][11],    		   
    		    	&data[N_PROTEINS][j][12],
    		    	&data[N_PROTEINS][j][13]);    		   
    		    if (num>10) { N_DECOYS[N_PROTEINS]++; j++; }
	      	} 
		    } while (!feof(inp));
		    printf("FILE %2d: %s, %d decoys\n", N_PROTEINS, filelist[i], N_DECOYS[i]);		    
		    N_PROTEINS++;
		  }
      fclose(inp);
    }
		printf("\n");
		
		for (i=0;i<N_PROTEINS;i++) {
			min = data[i][0][1];
			nat = 0;
			for (j=0;j<N_DECOYS[i];j++) {
				if (data[i][j][1]>min) { min=data[i][j][1]; nat=j; }
			}
			native[i]=nat;
			printf("PROTEIN: %3d  NATIVE: %3d %7.2f\n", i, nat, min);
		}
		
		
		printf("\n");

		printf("%d files read\n", N_PROTEINS);
		
	  parameters = (double*)calloc(sizeof(double)*(N_WEIGHTS+5),1);

		for (i=0;i<N_WEIGHTS+5;i++)
		  parameters[i] = 1.0;		
}

// calculates an energy for k-th protein j-th decoy
double E(int k, int j)
{
  double sum;
  int l;
  
    sum = 0.0;

//printf("ENE %d %d\n", k, j);
      for (l=0;l<7;l++) {
        sum += parameters[l]*data[k][j][6+l]; // start from column # 6
//printf("%5d %lf sum %lf   \t", l, parameters[l]*data[k][j][6+l], sum);
      }
//printf("sum: %lf\n", sum);
      
  return sum;
    
}

double RMS(int k, int j) 
{
  return 1.0-data[k][j][1];
//  return data[k][j][12];
}


void dump(int k)
{
  int j;
  
    printf("\n\n");

//    for (k=0;k<n;k++) {
      for (j=0;j<N_DECOYS[k];j++) {
        printf("%10.4lf %10.4lf\n", RMS(k,j), E(k,j));
      }
//    }
    
    printf("\n\n");
}

// calculates a correlation coefficient
double CC1(int k) 
{
  double sumx, sumy, sumx2, sumy2, sumxy, sumdx2, sumdy2, sumdxdy, id, cc;
  int j;

    id = 1./(double)N_DECOYS[k];
    

      sumx = 0.0;
      for (j=0;j<N_DECOYS[k];j++) {
        sumx += RMS(k,j);
      }
      
      sumx2 = 0.0;
      for (j=0;j<N_DECOYS[k];j++) {
        sumx2 += RMS(k,j)*RMS(k,j);
      }

      sumy = 0.0;
      for (j=0;j<N_DECOYS[k];j++) {
        sumy += E(k,j);
      }
      
      sumy2 = 0.0;
      for (j=0;j<N_DECOYS[k];j++) {
        sumy2 += E(k,j)*E(k,j);
      }

      sumxy = 0.0;
      for (j=0;j<N_DECOYS[k];j++) {
        sumxy += E(k,j)*RMS(k,j);
      }

      sumdx2 = sumx2 - (sumx*sumx)*id;
      sumdy2 = sumy2 - (sumy*sumy)*id;
      sumdxdy = sumxy - (sumx*sumy)*id;
      
      cc = sumdxdy/sqrt(sumdx2*sumdy2);
    
  return cc;    
}

// calculates a correlation coefficient
double CC(void) 
{
  double tot;
  int k;

    tot = 0.0;
    
    for (k=0;k<N_PROTEINS;k++) {
      tot += CC1(k);
    }
    
    tot *= (1./N_PROTEINS);
    
  return tot;    
}

void printCC(int init)
{
  int k;
  double tot;
  
    printf("\n");
    
    tot = 0.0;
    for (k=0;k<N_PROTEINS;k++) {
      printf("CC %d: %8.3lf, R2: %8.3lf  ", k, CC1(k), CC1(k)*CC1(k));
      if (init) initvalues[k][0] = CC1(k);
      else printf("before: %8.3f, diff: %8.3f", initvalues[k][0], CC1(k)-initvalues[k][0]);
      printf("\n");
      tot += CC1(k);
    }
    tot /= (double)N_PROTEINS;
    
    printf("\nAVERAGE CC: %8.3lf, R2: %8.3lf\n\n", tot, tot*tot);
    if (init) initvalues[k][0] = tot;    	
    else printf("before: %8.3f, diff: %8.3f", initvalues[k][0], tot-initvalues[k][0]);    
}



double G1(void)
{
  return (1./(1.+CC()));
}

double G2(void)
{
  double tot, sum, sumw, tmp, rms;
  int j, k, l;
  
    tot = 0.0;
    for (k=0;k<N_PROTEINS;k++) {
      sum = 0.0;
      for (j=0;j<N_DECOYS[k];j++) {
        tmp = RMS(k,j)-parameters[7]*E(k,j)+parameters[8];
        tmp *= tmp;
        tmp /= RMS(k,j);
        sum += tmp;
      }
      sum /= (double)N_DECOYS[k];
      tot += sum;
    }
    
  return tot;
}

// calculate z-score
double Z1(int k)
{
  double sume, sume2, enat, id, sig;
  int j;

    id = 1./(double)N_DECOYS[k];
    

      sume = 0.0;
      sume2 = 0.0;
      
      for (j=0;j<N_DECOYS[k];j++) {
        sume += E(k,j);
        sume2 += E(k,j)*E(k,j);
      }
      sume *= id;
      sume2 *= id;
      
			enat = E(k,native[k]);

			sig = sume2-sume*sume;
		
//printf("enat: %f sig: %f\n", enat, sig);

/*		if (sig>0.)*/ return (sume-enat)/sqrt(sig);    
			
//  return 0.0;    
	
}

double Z(void) 
{
  double tot;
  int k;

    tot = 0.0;
    
    for (k=0;k<N_PROTEINS;k++) {
      tot += Z1(k);
    }
    
    tot *= (1./N_PROTEINS);
    
  return tot;    
}


double G3(void)
{
  return (1./(1.+Z()));
}


void printZ(int init)
{
  int k;
  double tot;
  
    printf("\n");
    
    tot = 0.0;
    for (k=0;k<N_PROTEINS;k++) {
      printf("Z-score %d: %8.3lf ", k, Z1(k));
      if (init) initvalues[k][1] = Z1(k);
	    else printf("before: %8.3f, diff: %8.3f", initvalues[k][1], Z1(k)-initvalues[k][1]);    
      tot += Z1(k);
      printf("\n");
    }
    tot /= (double)N_PROTEINS;
    
    printf("\nAVERAGE Z-score: %8.3lf\n\n", tot);
    if (init) initvalues[k][1] = tot;
	  else printf("before: %8.3f, diff: %8.3f", initvalues[k][1], tot-initvalues[k][1]);    
}


#include "cfortran.h"
#include "minuitfcn.h"
#include "minuit.h"

#define MAXLINE 256

typedef struct {
  char *name;
  double value;
  double error;
  double min;
  double max;
} par_t;

par_t pars[] = {

  { "W4",   1.0, 1e-2, -100., 100. },
  { "W5",   1.0, 1e-2, -100., 100. },
  { "W6",   1.0, 1e-2, -100., 100. },
  { "W7",   1.0, 1e-2, -100., 100. },
  { "W8",   1.0, 1e-2, -100., 100. },
  { "W9",   1.0, 1e-2, -100., 100. },
  { "W10",   1.0, 1e-2, -100., 100. },
  { "W11",   1.0, 1e-2, -100., 100. },
  { "W12",   1.0, 1e-2, -100., 100. },
  { "eta",  1.0, 1e-2, -1000., 1000. },
  { "bk", 	1.0, 1e-2, -1000., 1000. }  

  
/*
  { "w1",   0.0, 1e-2, -0.1, 100. },
  { "w2",   1.0, 1e-2, -100., 100. },
  { "w3",   1.0, 1e-2, -100., 100. },
  { "w4",   1.0, 1e-2, -100., 100. },
  { "w5",   1.0, 1e-2, -100., 100. },
  { "w6",   0.0, 1e-2, -0.1, 100. },
  { "w7",   0.0, 1e-2, -0.1, 100. },
  { "eta",  0.0, 1e-2, -0.1, 100. },
  { "bk", 	0.0, 1e-2, -0.1, 100. }
*/  
};

void fcn (int npar, double* grad, double* fcnval, double* xval, int iflag, void* futil)
{
	int i;
	
	for (i=0;i<7;i++)
	  parameters[i] = xval[i];
  
  *fcnval = G1()*G3();
  
  /*(xval[0] - 2.0) * (xval[0] - 2.0) +
	    (xval[1] - 3.0) * (xval[1] - 3.0) +
	    (xval[2] + 5.0) * (xval[2] + 5.0); */
}

int main (int argc, char *argv[])
{
  int i, dummy;
  char *command = NULL;

  setbuf(stdout,0);
  printf("start\n");
  
  read_data(argv[1]);
	
  printf("\n");
	
  command = (char *) malloc (MAXLINE * sizeof (char));

  MNINIT (5,6,7);

  for (i = 0; i < 7; i++) {
    MNPARM (i+1, pars[i].name, pars[i].value, pars[i].error, pars[i].min, pars[i].max, dummy);
  }
	
  printCC(1);
  printf("\n\n");
  printZ(1);

  dump(32);
  
//  snprintf (command, MAXLINE, "SET STRATEGY 2");
//  i = 2;
//  MNCOMD (minuitfcn, command, i, NULL);

  snprintf (command, MAXLINE, "MIGRAD");
//  snprintf (command, MAXLINE, "SCAN");
  MNCOMD (minuitfcn, command, dummy, NULL);

//  snprintf (command, MAXLINE, "IMPROVE");
//  MNCOMD (minuitfcn, command, dummy, NULL);

  printCC(0);
  printf("\n\n");
  printZ(0);
  
  dump(32);

  exit (0);
}






