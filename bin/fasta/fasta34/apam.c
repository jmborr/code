/*	pam.c	19-June-86
	copyright (c) 1987 William R. Pearson
	read in the alphabet and pam matrix data
	designed for universal matcher

	This version reads BLAST format (square) PAM files
*/

/* $Name: fa34t20b3 $ - $Id: apam.c,v 1.17 2002/10/09 13:58:12 wrp Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

int
initpam (char *mfname, struct pstruct *ppst)
{
   char    line[512], *lp;
   int     i, iaa, ipam;
   int *hsq, nsq;
   int *sascii;
   char *sq;
   int ess_tmp, max_val;
   int have_es = 0;
   FILE   *fmat;

   if ((fmat = fopen (mfname, "r")) == NULL)
   {
      printf (" cannot open scoring matrix file %s\n", mfname);
      return 0;
   }

   hsq = ppst->hsq;
   nsq = ppst->nsq;
   sq = ppst->sq;

   while (fgets (line, sizeof(line), fmat) != NULL && line[0]=='#');

   /* decide whether this is a protein or DNA matrix */
   if (ppst->dnaseq == 1) sascii = &nascii[0];
   else sascii = &aascii[0];

   ess_tmp = sascii[','];

/*	clear out sascii	*/
   for (i = 0; i <= AAMASK; i++) sascii[i] = NA;

/*	set end of line stop	*/
   sascii[0] = sascii['\r'] = sascii['\n'] = EL;

   sascii[','] = ess_tmp;

/* read the alphabet */
   sq[0] = '\0';
   for (i = 0, nsq = 1; line[i]; i++) {
     if (line[i] == '*') have_es = 1;
     if (line[i] > ' ') sq[nsq++] = toupper (line[i]);
   }
   sq[nsq]='\0';
   nsq--;

/* set end of sequence stop */
   fprintf(stderr,"sq[%d]: %s\n",nsq,sq+1);

/* initialize sascii */
   for (iaa = 1; iaa <= nsq; iaa++) {
     sascii[sq[iaa]] = iaa;
   }
   if (ppst->dnaseq==1) {
     sascii['U'] = sascii['T'];
     sascii['u'] = sascii['t'];
   }

/* setup hnt values */
   hsq[0] = 0;
   for (iaa = 1; iaa <= nsq; iaa++) {
     hsq[iaa]=iaa;
   }
   if (ppst->dnaseq==1) {
     hsq[sascii['R']]=hsq[sascii['M']]=hsq[sascii['W']]=hsq[sascii['A']];
     hsq[sascii['D']]=hsq[sascii['H']]=hsq[sascii['V']]=hsq[sascii['A']];
     hsq[sascii['N']]=hsq[sascii['X']]=hsq[sascii['A']];
     hsq[sascii['Y']]=hsq[sascii['S']]=hsq[sascii['B']]=hsq[sascii['C']];
     hsq[sascii['K']]=hsq[sascii['G']];
   }
   else if (ppst->dnaseq == -1 || (ppst->nsq >= 20 && ppst->nsq <= 23)) {
     hsq[sascii['B']] = hsq[sascii['N']];
     hsq[sascii['Z']] = hsq[sascii['E']];
     hsq[sascii['X']] = hsq[sascii['A']];
   }
   /* here if non-DNA, non-protein sequence */
   else ppst->dnaseq = 2;

   max_val = -1;
   for (iaa = 1, ipam = 0; iaa <= nsq; iaa++) {
     if (fgets(line,sizeof(line),fmat)==NULL) {
       fprintf (stderr," error reading pam line: %s\n",line);
       exit (1);
     }
     /*     fprintf(stderr,"%d/%d %s",iaa,nsq,line); */
     strtok(line," \t\n");
     for (i = 0; i < iaa; i++,ipam++) {
       lp=strtok(NULL," \t\n");
       pam[ipam]=atoi(lp);
       if (pam[ipam] > max_val) max_val = pam[ipam];
     }
   }

   /*
   for (iaa = 1, ipam = 0; iaa <= nsq; iaa++) {
     fprintf(stderr,"%c",sq[iaa]);
     for (i=0; i< iaa; i++) 
       fprintf(stderr," %2d",pam[ipam++]);
     fprintf(stderr,"\n");
   }
   */

   if (have_es==0) {
     sascii['*']=nsq;
     for (i=1; i<=nsq; i++) pam[ipam++]= -1;
     pam[ipam]= max_val/2;
     nsq++;
     sq[nsq]='*';
     sq[nsq+1]='\0';
   }

   ppst->nsq = nsq;
   ppst->sqx[0]='\0';
   for (i=1; i<= nsq; i++) {
     ppst->sqx[i] = sq[i];
     ppst->sqx[i+nsq] = tolower(sq[i]);
     if (sascii[aa[i]] < NA && sq[i] >= 'A' && sq[i] <= 'Z')
       sascii[aa[i] - 'A' + 'a'] = sascii[aa[i]]+nsq;
   }
   ppst->nsqx = nsq*2;
   
   strncpy (ppst->pamfile, mfname, 40);
   ppst->pamfile[39]='\0';
   fclose (fmat);
   return 1;
}

/* make a DNA scoring from +match/-mismatch values */

void mk_n_pam(int *arr,int siz, int mat, int mis)
{
  int i, j, k;
  /* current default match/mismatch values */
  int max_mat = +5;
  int min_mis = -4;
  float f_val, f_scale;
  
  f_scale = (float)(mat - mis)/(float)(max_mat - min_mis);
  
  k = 0;
  for (i = 0; i<nnt-1; i++)
    for (j = 0; j <= i; j++ ) {
      if (arr[k] == max_mat) arr[k] = mat;
      else if (arr[k] == min_mis) arr[k] = mis;
      else if (arr[k] != -1) { 
	f_val = (arr[k] - min_mis)*f_scale + 0.5;
	arr[k] = f_val + mis;
      }
      k++;
    }
}

struct std_pam_str {
  char abbrev[6];
  char name[10];
  int *pam;
  int gdel, ggap;
};

static
struct std_pam_str std_pams[] = {
  {"P120","PAM120",apam120,-20,-3},
  {"P250","PAM250",apam250,-12,-2},
  {"P10","MD_10",a_md10,-27,-4},
  {"M10","MD_10",a_md10,-27,-4},
  {"P20","MD_20",a_md20,-26,-4},
  {"M20","MD_20",a_md20,-26,-4},
  {"P40","MD_40",a_md40,-25,-4},
  {"M40","MD_40",a_md40,-25,-4},
  {"BL50","BL50",abl50,-12,-2},
  {"BL62","BL62",abl62,-8,-1},
  {"BL80","BL80",abl80,-12,-2},
  {"\0","\0",NULL,0,0},
};

int
standard_pam(char *smstr, struct pstruct *ppst) {

  struct std_pam_str *std_pam_p;

  for (std_pam_p = std_pams; std_pam_p->abbrev[0]; std_pam_p++ ) {
    if (strcmp(smstr,std_pam_p->abbrev)==0) {
      pam = std_pam_p->pam;
      strncpy(ppst->pamfile,std_pam_p->name,40);
      ppst->pamfile[39]='\0';
#ifndef GAP_OPEN
      ppst->gdelval = std_pam_p->gdel;
#else
      ppst->gdelval = std_pam_p->gdel-std_pam_p->ggap;
#endif
      ppst->ggapval = std_pam_p->ggap;
      return 1;
    }
  }
  return 0;
}

void
build_xascii(int *qascii) {
  int i;
  int comma_char;

  /* preserve comma character */
  comma_char = qascii[','];

  for (i=1; i<128; i++) {
    qascii[i]=NA;
  }
  /* range of values in aax, ntx is from 1..naax,nntx - 
     do not zero-out qascii[0] - 9 Oct 2002 */

  for (i=1; i<naax; i++) {
    qascii[aax[i]]=aax[i];
  }

  for (i=1; i<nntx; i++) {
    qascii[ntx[i]]=ntx[i];
  }

  qascii['\n']=qascii['\r']=qascii[0] = EL;
  qascii[','] = comma_char;
}




/* 
   checks for lower case letters in *sq array;
   if not present, map lowercase to upper
*/
void
init_ascii(int is_ext, int *sascii, int is_dna) {

  int isq, have_lc;
  char *sq;
  int nsq;
  
  if (is_dna==SEQTYPE_UNK) return;

  if (is_dna==SEQTYPE_DNA) {
    if (is_ext) { 
      sq = &ntx[0];
      nsq = nntx;
    }
    else {sq = &nt[0]; nsq = nnt;}
  }
  else {
    if (is_ext) { sq = &aax[0]; nsq = naax; }
    else {sq = &aa[0]; nsq = naa;}
  }


/* initialize sascii from sq[], checking for lower-case letters */
  have_lc = 0;
  for (isq = 1; isq <= nsq; isq++) {
     sascii[sq[isq]] = isq;
     if (sq[isq] >= 'a' && sq[isq] <= 'z') have_lc = 1;
  }

  /* no lower case letters in alphabet, map lower case to upper */
  if (have_lc != 1) { 
    for (isq = 1; isq <= nsq; isq++) {
      if (sq[isq] >= 'A' && sq[isq] <= 'Z') sascii[sq[isq]-'A'+'a'] = isq;
    }
    if (is_dna==1) sascii['u'] = sascii['t'];
  }
}
