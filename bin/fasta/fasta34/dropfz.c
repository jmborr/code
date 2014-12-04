
/* copyright (c) 1998, 1999 William R. Pearson and the U. of Virginia */

/* $Name: fa34t20b3 $ - $Id: dropfz.c,v 1.19 2002/08/21 14:20:22 wrp Exp $ */

/* implements an improved version of the fasty algorithm, see:

   W. R. Pearson, T. Wood, Z. Zhang, A W. Miller (1997) "Comparison of
   DNA sequences with protein sequences" Genomics 46:24-36

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"
#define XTERNAL
#include "upam.h"
#include "uascii.h"

/* globals for fasta */
#define MAXWINDOW 64

#ifndef MAXSAV
#define MAXSAV 10
#endif

#ifndef ALLOCN0
static char *verstr="3.36 June 2000";
#else
static char *verstr="3.36an0 June 2000";
#endif

struct dstruct		/* diagonal structure for saving current run */
{			
   int     score;	/* hash score of current match */
   int     start;	/* start of current match */
   int     stop;	/* end of current match */
   struct savestr *dmax;   /* location in vmax[] where best score data saved */
};

struct savestr
{
   int     score;		/* pam score with segment optimization */
   int     score0;		/* pam score of best single segment */
   int     gscore;		/* score from global match */
   int     dp;			/* diagonal of match */
   int     start;		/* start of match in lib seq */
   int     stop;		/* end of match in lib seq */
};

void savemax();
void kpsort();

struct sx_s {int C1, C2, C3, I1, I2, I3, flag; };

struct wgt { int  iii, ii, iv;};
struct wgtc {char c2, c3, c4, c5;};

#if defined(FASTX) || defined(TFASTX)
extern int aatran(unsigned char *ntseq, unsigned char *aaseq, int maxs, int frame);
extern int saatran(unsigned char *ntseq, unsigned char *aaseq, int maxs, int frame);

struct f_struct {
  struct dstruct *diag;
  struct savestr vmax[MAXSAV];	/* best matches saved for one sequence */
  struct savestr *vptr[MAXSAV];
  struct savestr *lowmax;
  int ndo;
  int noff;
  int hmask;			/* hash constants */
  int *pamh1;			/* pam based array */
  int *pamh2;			/* pam based kfact array */
  int *link, *harr;		/* hash arrays */
  int kshft;			/* shift width */
  int nsav, lowscor;		/* number of saved runs, worst saved run */
#ifdef FASTX
  unsigned char *aa0x, *aa0v;	/* aa0x - 111122223333 */
#endif				/* aa0v - computed codons */
#ifdef TFASTX
  unsigned char *aa1x, *aa1v;	/* aa1x - 111122223333 */
#endif				/* aa1v - computed codons */
#if defined(FASTX) || defined(TFASTX)
  struct sx_s *cur;
  struct wgt **weight0;
  struct wgt **weight1;
  struct wgtc **weight_c;
#endif
  int *waa;
  int *res;
  int max_res;
};

static int dmatchx(unsigned char *aa0, int n0,
		   unsigned char *aa1, int n1,
		   int hoff, int window, 
		   int **pam2, int gdelval, int ggapval, int gshift,
		   struct f_struct *f_str);
#endif

/* initialize for fasta */
/* modified 30-August-1999 by Zheng Zhang to work with an extended alphabet */
/* Assume naa=47, and wgts[47][23] matches both upper and lower case
amoino acids with another amino acid.  And also assume the DNA letter
does not have upper/lower case difference.  If you also allow DNA
sequence to be upper/lower case letters, more need be changed. Not
only here, but also in the alignment code, the way that pack a codon
into a number between 0-63 need be changed. */

/* modified so that if **weightci==NULL, do not fiddle with characters */

void
init_weights(struct wgt ***weighti, struct wgtc ***weightci,
	     int **wgts, int gshift, int gsubs, int naa)
{
  int i, j, do_wgtc=0;
  int aa, b, a, x, y, z;
  int *wwt, e;
  struct wgt **weight;
  struct wgtc **weightc;
  char aacmap[64];
  int temp[47][64]; /*change*/
  char le[47][64];


  if ((*weighti=(struct wgt **)calloc((size_t)(naa+1),sizeof(struct wgt *)))
      ==NULL) {
    fprintf(stderr," cannot allocate weights array: %d\n",naa);
    exit(1);
  }

  weight = *weighti;
  for (aa=0; aa <= naa; aa++) {
    if ((weight[aa]=(struct wgt *)calloc((size_t)256,sizeof(struct wgt)))
	==NULL) {
      fprintf(stderr," cannot allocate weight[]: %d/%d\n",aa,naa);
      exit(1);
    }
  }

  if (weightci !=NULL) {
    if ((*weightci=(struct wgtc **)calloc((size_t)(naa+1),
					  sizeof(struct wgtc *)))==NULL) {
      fprintf(stderr," cannot allocate weight_c array: %d\n",naa);
      exit(1);
    }
    weightc = *weightci;

    for (aa=0; aa <= naa; aa++) {
      if ((weightc[aa]=(struct wgtc *)calloc((size_t)256,sizeof(struct wgtc)))
	  ==NULL) {
	fprintf(stderr," cannot allocate weightc[]: %d/%d\n",aa,naa);
	exit(1);
      }
    }
    do_wgtc = 1;
  }
  else do_wgtc = 0;

  aagetmap(aacmap,64);

  for (aa = 0; aa <= naa; aa++) { /* change*/
      wwt = wgts[aa];
      for (i = 0; i < 64; i++) {	/* j iterates through the codons */
	  x = -1000;
	  y = i;
	  for (j = 0; j < 64; j++) {	/* j iterates through the codons */
	      z = ((~i & j) | (i & ~j));
	      b = 0;		/* score = 0 */
	      if (z % 4) b-= gsubs;
	      if (z /16) b-= gsubs;
	      if ((z /4) % 4) b -= gsubs;   
	      b += wwt[aascii[aacmap[j]]];  /* add the match score for char j*/
	      if (b > x) {
		x = b;		/* x has the score */
		y = j;		/* y has the character */
	      }
	  }
	  /*	  if (y < 0 || y > 63) printf("%d %d %d %d ",aa, i, x, y); */
	  temp[aa][i] = x;
	  le[aa][i] = y;
      }
      /*            printf("\n"); */
  }

  for (aa= 0; aa <= naa; aa++) { 
      wwt = temp[aa];
      for (i = 0; i < 256; i++) {
          for (x=-100,b = 0; b < 4; b++) {
              z = (i/ (1 << ((b+1)*2)))*(1<<(b*2))+(i%(1<<(b*2)));
	      if (x < (e=wwt[z])) {
		  x = e;
		  if (do_wgtc) weightc[aa][i].c4 = aacmap[le[aa][z]];
	      }
          }
          weight[aa][i].iv=x-gshift;
          weight[aa][i].iii = wwt[i%64];

	  if (do_wgtc) {
	    weightc[aa][i].c5 = aacmap[le[aa][i%64]];
	    weightc[aa][i].c3 = aacmap[i%64];
	  }
          x = i %16;
          for (y = -100, b = 0; b < 3; b++) {
              z = ((x >> (b*2)) << (b*2+2)) + (x % (1 << (b*2))); 
              for (a = 0; a < 4; a++) {
		  if ((e =wwt[z+(a<<(b*2))]) > y) {
		      y = e;
		      if (do_wgtc) 
			weightc[aa][i].c2 = aacmap[le[aa][z+(a<<(b*2))]];
		  }
              }
          }
          weight[aa][i].ii = y-gshift;
      }
  }
}

void
free_weights(struct wgt ***weighti0, struct wgt ***weighti1, 
	     struct wgtc ***weightci, int naa)
{
  int aa;
  struct wgt **weight0;
  struct wgt **weight1;
  struct wgtc **weightc;

  weight0 = *weighti0;
  weight1 = *weighti1;
  weightc = *weightci;

  for (aa=0; aa <= naa; aa++) {free(weight0[aa]);}
  for (aa=0; aa <= naa; aa++) {free(weight1[aa]);}
  for (aa=0; aa <= naa; aa++) {free(weightc[aa]);}

  free(weight0);
  free(weight1);
  free(weightc);
}

void
init_work (unsigned char *aa0, int n0, 
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
   int mhv, phv;
   int hmax;
   int i0, hv;
   int pamfact;
   int btemp;
   struct f_struct *f_str;
   struct bdstr *bss;
   /* these used to be globals, but do not need to be */
   int ktup, fact, kt1, lkt;

   int maxn0;
   int *pwaa;
   int i, j, q;
   struct swstr *ss, *r_ss;
   int *waa;
   int *res;
   int nsq, ip, *hsq, naat;
#if defined(FASTX)
   int last_n0, itemp, dnav;
   unsigned char *fd, *fs, *aa0x, *aa0v;
   int n0x, n0x3;
#endif

  if (ppst->ext_sq_set) {
    nsq = ppst->nsqx; ip = 1;
    hsq = ppst->hsqx;
  }
  else {
    nsq = ppst->nsq; ip = 0;
    hsq = ppst->hsq;
  }

   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   btemp = 2 * ppst->param_u.fa.bestoff / 3 +
      n0 / ppst->param_u.fa.bestscale +
      ppst->param_u.fa.bkfact *
      (ppst->param_u.fa.bktup - ppst->param_u.fa.ktup);
   btemp = min (btemp, ppst->param_u.fa.bestmax);
   if (btemp > 3 * n0) btemp = 3 * shscore(aa0,n0,ppst->pam2[0]) / 5;

   ppst->param_u.fa.cgap = btemp + ppst->param_u.fa.bestoff / 3;
   if (ppst->param_u.fa.optcut_set != 1)
#if !defined(TFASTX)
      ppst->param_u.fa.optcut = (btemp*5)/4;
#else
      ppst->param_u.fa.optcut = (btemp*4)/3;
#endif

   ppst->param_u.fa.pgap = ppst->gdelval + ppst->ggapval;
   pamfact = ppst->param_u.fa.pamfact;
   ktup = ppst->param_u.fa.ktup;
   fact = ppst->param_u.fa.scfact * ktup;

#ifdef FASTX
   /* before hashing, we must set up some space and translate the sequence */

   maxn0 = n0 + 2;
   if ((aa0x =(unsigned char *)calloc((size_t)maxn0,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "cannot allocate aa0x array %d\n", maxn0);
     exit (1);
   }
   aa0x++;
   f_str->aa0x = aa0x;


   if ((aa0v =(unsigned char *)calloc((size_t)maxn0,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "cannot allocate aa0v array %d\n", maxn0);
     exit (1);
   }
   aa0v++;
   f_str->aa0v = aa0v;

   /* make a precomputed codon number series */
   dnav = (hnt[aa0[0]]<<2) + hnt[aa0[1]];
   for (i=2; i<n0; i++) {
     dnav = ((dnav<<2)+hnt[aa0[i]])&255;
     aa0v[i-2]=dnav;
   }

   last_n0 = 0;
   for (itemp=0; itemp<3; itemp++) {
     n0x=saatran(aa0,&aa0x[last_n0],n0,itemp);
     /*         for (i=0; i<n0x; i++) {
	   fprintf(stderr,"%c",aa[aa0x[last_n0+i]]);
	   if ((i%60)==59) fprintf(stderr,"\n");
	   }
	   fprintf(stderr,"\n");
	   */
     last_n0 += n0x+1;
   }

   /*     fprintf(stderr,"\n"); */
   n0x = n0;
   n0x3 = n0x/3;

   /* now switch aa0 and aa0x for hashing functions */
   fs = aa0;
   aa0 = aa0x;
   aa0x = fs;
#endif
#if defined(FASTX) || defined(TFASTX)
   if (ppst->ext_sq_set) naat = 46;
   else naat = 23;

   init_weights(&f_str->weight0, NULL,
		ppst->pam2[ip],-ppst->gshift,-ppst->gsubs,naat);
   init_weights(&f_str->weight1, &f_str->weight_c,
		ppst->pam2[0],-ppst->gshift,-ppst->gsubs,naat);
#endif

   if (pamfact == -1)
      pamfact = 0;
   else if (pamfact == -2)
      pamfact = 1;

   for (i0 = 1, mhv = -1; i0 < ppst->nsq; i0++)
     if (hsq[i0] < NMAP && hsq[i0] > mhv)
       mhv = ppst->hsq[i0];

   if (mhv <= 0)
   {
      fprintf (stderr, " maximum hsq <=0 %d\n", mhv);
      exit (1);
   }

   for (f_str->kshft = 0; mhv > 0; mhv /= 2) f_str->kshft++;

/*      kshft = 2;	*/
   kt1 = ktup - 1;
   hv = 1;
   for (i0 = 0; i0 < ktup; i0++)
      hv = hv << f_str->kshft;
   hmax = hv;
   f_str->hmask = (hmax >> f_str->kshft) - 1;


   if ((f_str->harr = (int *) calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate hash array\n");
     exit (1);
   }
   if ((f_str->pamh1 = (int *) calloc (ppst->nsq+1, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate pamh1 array\n");
     exit (1);
   }
   if ((f_str->pamh2 = (int *) calloc (hmax, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate pamh2 array\n");
     exit (1);
   }
   if ((f_str->link = (int *) calloc (n0, sizeof (int))) == NULL) {
     fprintf (stderr, " cannot allocate hash link array");
     exit (1);
   }

   for (i0 = 0; i0 < hmax; i0++)
      f_str->harr[i0] = -1;
   for (i0 = 0; i0 < n0; i0++)
      f_str->link[i0] = -1;

   /* encode the aa0 array */

   phv = hv = 0;
   lkt = kt1;
   for (i0 = 0; i0 < min(n0,lkt); i0++) {
     if (hsq[aa0[i0]] >= NMAP) {hv=phv=0; lkt=i0+ktup; continue;}
     hv = (hv << f_str->kshft) + ppst->hsq[aa0[i0]];
     phv += ppst->pam2[ip][aa0[i0]][aa0[i0]] * ktup;
   }

   for (; i0 < n0; i0++) {
     if (hsq[aa0[i0]] >= NMAP) {
       hv=phv=0;
       /* restart hv, phv calculation */
       for (lkt = i0+kt1; i0 < min(lkt,n0); i0++) {
	 if (hsq[aa0[i0]] >= NMAP) {hv=phv=0; lkt=i0+ktup; continue;}
	 hv = (hv << f_str->kshft) + hsq[aa0[i0]];
	 phv += ppst->pam2[ip][aa0[i0]][aa0[i0]]*ktup;
       }
       continue;
     }
     hv = ((hv & f_str->hmask) << f_str->kshft) + ppst->hsq[aa0[i0]];
     f_str->link[i0] = f_str->harr[hv];
     f_str->harr[hv] = i0;
     if (pamfact) {
       f_str->pamh2[hv] = (phv += ppst->pam2[ip][aa0[i0]][aa0[i0]] * ktup);
       if (hsq[aa0[i0-kt1]]<NMAP)
	 phv -= ppst->pam2[ip][aa0[i0 - kt1]][aa0[i0 - kt1]] * ktup;
     }
     else f_str->pamh2[hv] = fact * ktup;
   }

/* this has been modified from 0..<ppst->nsq to 1..<=ppst->nsq because the
   pam2[0][0] is now undefined for consistency with blast
*/

   if (pamfact)
      for (i0 = 1; i0 <= ppst->nsq; i0++)
	 f_str->pamh1[i0] = ppst->pam2[ip][i0][i0] * ktup;
   else
      for (i0 = 1; i0 <= ppst->nsq; i0++)
	 f_str->pamh1[i0] = fact;

   f_str->ndo = 0;	/* used to save time on diagonals with long queries */


#ifndef ALLOCN0
   if ((f_str->diag = (struct dstruct *) calloc ((size_t)MAXDIAG,
						 sizeof (struct dstruct)))==NULL) {
      fprintf (stderr," cannot allocate diagonal arrays: %ld\n",
	      (long) MAXDIAG *sizeof (struct dstruct));
      exit (1);
     };
#else
   if ((f_str->diag = (struct dstruct *) calloc ((size_t)n0,
					      sizeof (struct dstruct)))==NULL) {
      fprintf (stderr," cannot allocate diagonal arrays: %ld\n",
	      (long)n0*sizeof (struct dstruct));
      exit (1);
     };
#endif

#ifdef FASTX
   /* done hashing, now switch aa0, aa0x back */
   fs = aa0;
   aa0 = aa0x;
   aa0x = fs;
#endif

#if !defined(FASTX) && !defined(TFASTX)
   /* allocate space for the scoring arrays */
   maxn0 = n0 + 2;
   if ((ss = (struct swstr *) calloc (maxn0, sizeof (struct swstr)))
	 == NULL) {
     fprintf (stderr, "cannot allocate ss array %3d\n", n0);
     exit (1);
   }
   ss++;
   f_str->ss = ss;

   if ((r_ss = (struct swstr *) calloc (maxn0, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, "cannot allocate r_ss array %3d\n", n0);
     exit (1);
   }
   r_ss++;
   f_str->r_ss = r_ss;
#endif
#ifdef TFASTX
   if ((f_str->aa1x =(unsigned char *)calloc((size_t)ppst->maxlen+4,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "cannot allocate aa1x array %d\n", ppst->maxlen+4);
     exit (1);
   }
   f_str->aa1x++;

   if ((f_str->aa1v =(unsigned char *)calloc((size_t)ppst->maxlen+4,
					     sizeof(unsigned char))) == NULL) {
     fprintf (stderr, "cannot allocate aa1v array %d\n", ppst->maxlen+4);
     exit (1);
   }
   f_str->aa1v++;

#endif

   if ((waa= (int *)malloc (sizeof(int)*(ppst->nsq+1)*n0)) == NULL) {
     fprintf(stderr,"cannot allocate waa struct %3d\n",ppst->nsq*n0);
     exit(1);
   }

   pwaa = waa;
   for (i=0; i<=ppst->nsq; i++) {
     for (j=0;j<n0; j++) {
       *pwaa = ppst->pam2[ip][i][aa0[j]];
       pwaa++;
     }
   }
   f_str->waa = waa;

#ifndef TFASTX
   maxn0 = max(2*n0,MIN_RES);
#else
   maxn0 = max(4*n0,MIN_RES);
#endif
   if ((res = (int *)calloc((size_t)maxn0,sizeof(int)))==NULL) {
     fprintf(stderr,"cannot allocate alignment results array %d\n",maxn0);
     exit(1);
   }
   f_str->res = res;
   f_str->max_res = maxn0;

   *f_arg = f_str;
}

/* pstring1 is a message to the manager, currently 512 */
/* pstring2 is the same information, but in a markx==10 format */
void
get_param (struct pstruct *pstr, char *pstring1, char *pstring2)
{
#ifdef FASTX
  char *pg_str="FASTY";
#else
#ifdef TFASTX
  char *pg_str="TFASTY";
#endif
#endif

#if !defined(FASTX) && !defined(TFASTX)
   if (!pstr->param_u.fa.optflag)
      sprintf (pstring1, "%s (%s) function [%s matrix (%d:%d)%s] ktup: %d\n join: %d, gap-pen: %d/%d, width: %3d",pg_str,verstr,
	       pstr->pamfile, pstr->pam_h,pstr->pam_l, 
	       (pstr->ext_sq_set) ? "xS":"\0",
	       pstr->param_u.fa.ktup, pstr->param_u.fa.cgap,
	       pstr->gdelval, pstr->ggapval, pstr->param_u.fa.optwid);
   else
      sprintf (pstring1, "%s (%s) function [optimized, %s matrix (%d:%d)%s] ktup: %d\n join: %d, opt: %d, gap-pen: %3d/%3d, width: %3d",pg_str,verstr,
	       pstr->pamfile, pstr->pam_h,pstr->pam_l, 
	       (pstr->ext_sq_set) ? "xS":"\0",
	       pstr->param_u.fa.ktup, pstr->param_u.fa.cgap,
	       pstr->param_u.fa.optcut, pstr->gdelval, pstr->ggapval,
	       pstr->param_u.fa.optwid);
#else  /* FASTX or TFASTX */
   if (!pstr->param_u.fa.optflag)
      sprintf (pstring1, "%s (%s) function [%s matrix (%d:%d)%s] ktup: %d\n join: %d, gap-pen: %d/%d, shift: %d subs: %d width: %3d",pg_str,verstr,
	       pstr->pamfile, pstr->pam_h,pstr->pam_l, 
	       (pstr->ext_sq_set) ? "xS":"\0",
	       pstr->param_u.fa.ktup, pstr->param_u.fa.cgap,
	       pstr->gdelval, pstr->ggapval, pstr->gshift, pstr->gsubs,
	       pstr->param_u.fa.optwid);
   else
      sprintf (pstring1, "%s (%s) function [optimized, %s matrix (%d:%d)%s] ktup: %d\n join: %d, opt: %d, gap-pen: %3d/%3d shift: %3d, subs: %3d width: %3d",pg_str,verstr,
	       pstr->pamfile, pstr->pam_h,pstr->pam_l,
	       (pstr->ext_sq_set) ? "xS":"\0",
	       pstr->param_u.fa.ktup, pstr->param_u.fa.cgap,
	       pstr->param_u.fa.optcut, pstr->gdelval, pstr->ggapval,
	       pstr->gshift,pstr->gsubs,pstr->param_u.fa.optwid);
#endif
   if (pstr->param_u.fa.iniflag) strcat(pstring1," init1");
   /*
   if (pstr->zsflag==0) strcat(pstring1," not-scaled");
   else if (pstr->zsflag==1) strcat(pstring1," reg.-scaled");
   */

   if (pstring2 != NULL) {
     sprintf (pstring2, "; pg_name: %s\n; pg_ver: %s\n; pg_matrix: %s (%d:%d)%s\n\
; pg_gap-pen: %d %d\n; pg_ktup: %d\n; pg_optcut: %d\n; pg_cgap: %d\n",
	      pg_str,verstr,pstr->pamfile, pstr->pam_h,pstr->pam_l,
	      (pstr->ext_sq_set) ? "xS":"\0",  pstr->gdelval,
              pstr->ggapval,pstr->param_u.fa.ktup,pstr->param_u.fa.optcut,
	      pstr->param_u.fa.cgap);
   }
}

void
close_work (unsigned char *aa0, int n0,
	    struct pstruct *ppst,
	    struct f_struct **f_arg)
{
  struct f_struct *f_str;
  int naat;

  f_str = *f_arg;

  if (f_str != NULL) {
#if !defined(FASTX) && !defined(TFASTX)
    f_str->ss--;
    free(f_str->ss);
    f_str->r_ss--;
    free(f_str->r_ss);
#else
   if (ppst->ext_sq_set) naat = 46;
   else naat = 23;
    free_weights(&f_str->weight0,&f_str->weight1,&f_str->weight_c,naat);
    free(f_str->cur);
#ifdef FASTX
    f_str->aa0v--;
    free(f_str->aa0v);
    f_str->aa0x--;
    free(f_str->aa0x);
#else	/* TFASTX */
    f_str->aa1x--;
    free(f_str->aa1x);
    f_str->aa1v--;
    free(f_str->aa1v);
#endif
#endif
    free(f_str->res);
    free(f_str->waa);
    free(f_str->diag);
    free(f_str->link);
    free(f_str->pamh2); 
    free(f_str->pamh1);
    free(f_str->harr);
    free(f_str);
    *f_arg = NULL;
  }
}

void do_fasta (unsigned char *aa0, int n0,
	       unsigned char *aa1, int n1,
	       struct pstruct *ppst, struct f_struct *f_str,
	       struct rstruct *rst, int *hoff)
{
   int     nd;		/* diagonal array size */
   int     lhval;
   int     kfact;
   int i;
   register struct dstruct *dptr;
   register int tscor;

#ifndef ALLOCN0
   register struct dstruct *diagp;
#else
   register int dpos;
   int     lposn0;
#endif
   struct dstruct *dpmax;
   register int lpos;
   int     tpos;
   struct savestr *vmptr;
   int     scor, tmp;
   int     im, ib, nsave;
   int     cmps ();		/* comparison routine for ksort */
   int ktup, kt1, *hsq, ip, lkt;
#ifdef FASTX
   int n0x31, n0x32;
   n0x31 = (n0-2)/3;
   n0x32 = n0x31+1+(n0-n0x31-1)/2;
#endif
#ifdef TFASTX
   unsigned char *fs, *fd;
   int n1x31, n1x32, last_n1, itemp;
   n1x31 = (n1-2)/3;
   n1x32 = n1x31+1+(n1-n1x31-1)/2;
#endif

  if (ppst->ext_sq_set) {
    ip = 1;
    hsq = ppst->hsqx;
  }
  else {
    ip = 0;
    hsq = ppst->hsq;
  }

   ktup = ppst->param_u.fa.ktup;
   kt1 = ktup-1;

   if (n1 < ktup) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     return;
   }

   if (n0+n1+1 >= MAXDIAG) {
     fprintf(stderr,"n0,n1 too large: %d, %d\n",n0,n1);
     rst->score[0] = rst->score[1] = rst->score[2] = -1;
     return;
   }

   f_str->noff = n0 - 1;

#ifdef ALLOCN0
   nd = n0;
#endif

#ifndef ALLOCN0
   nd = n0 + n1;
#endif

   dpmax = &f_str->diag[nd];
   for (dptr = &f_str->diag[f_str->ndo]; dptr < dpmax;)
   {
      dptr->stop = -1;
      dptr->dmax = NULL;
      dptr++->score = 0;
   }

   for (vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++)
      vmptr->score = 0;
   f_str->lowmax = f_str->vmax;
   f_str->lowscor = 0;


   /* start hashing */
   lhval = 0;
   lkt = kt1;
   for (lpos = 0; lpos < min(lkt,n1);) {
     if (hsq[aa1[lpos]]>=NMAP) {
       lhval = 0; lpos++; lkt=lpos+ktup; continue;
     }
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos++]];
   }

#ifndef ALLOCN0
   diagp = &f_str->diag[f_str->noff + lkt];
   for (; lpos < n1; lpos++, diagp++) {
     if (hsq[aa1[lpos]]>=NMAP) {
       lpos++ ; diagp++;
       while (lpos < n1 && hsq[aa1[lpos]]>=NMAP) {lpos++; diagp++;}
       lhval = 0;
     }
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval]; tpos >= 0; tpos = f_str->link[tpos]) {
       if ((tscor = (dptr = &diagp[-tpos])->stop) >= 0) {
#else
   lposn0 = f_str->noff + lpos;
   for (; lpos < n1; lpos++, lposn0++) {
     if (hsq[aa1[lpos]]>=NMAP) {lhval = 0; goto loopl;}
     lhval = ((lhval & f_str->hmask) << f_str->kshft) + ppst->hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval]; tpos >= 0; tpos = f_str->link[tpos]) {
       dpos = lposn0 - tpos;
       if ((tscor = (dptr = &f_str->diag[dpos % nd])->stop) >= 0) {
#endif
	 tscor += ktup;
	 if ((tscor -= lpos) <= 0) {
	   scor = dptr->score;
	   if ((tscor += (kfact = f_str->pamh2[lhval])) < 0 && f_str->lowscor < scor)
#ifdef ALLOCN0
	     savemax (dptr, dpos, f_str);
#else
	     savemax (dptr, f_str);
#endif
	     if ((tscor += scor) >= kfact) {
	       dptr->score = tscor;
	       dptr->stop = lpos;
	     }
	     else {
	       dptr->score = kfact;
	       dptr->start = (dptr->stop = lpos) - kt1;
	     }
	 }
	 else {
	   dptr->score += f_str->pamh1[aa0[tpos]];
	   dptr->stop = lpos;
	 }
       }
       else {
	 dptr->score = f_str->pamh2[lhval];
	 dptr->start = (dptr->stop = lpos) - kt1;
       }
     }				/* end tpos */

#ifdef ALLOCN0
      /* reinitialize diag structure */
   loopl:
     if ((dptr = &f_str->diag[lpos % nd])->score > f_str->lowscor)
       savemax (dptr, lpos, f_str);
     dptr->stop = -1;
     dptr->dmax = NULL;
     dptr->score = 0;
#endif
   }				/* end lpos */

#ifdef ALLOCN0
   for (tpos = 0, dpos = f_str->noff + n1 - 1; tpos < n0; tpos++, dpos--) {
     if ((dptr = &f_str->diag[dpos % nd])->score > f_str->lowscor)
       savemax (dptr, dpos, f_str);
   }
#else
   for (dptr = f_str->diag; dptr < dpmax;) {
     if (dptr->score > f_str->lowscor) savemax (dptr, f_str);
     dptr->stop = -1;
     dptr->dmax = NULL;
     dptr++->score = 0;
   }
   f_str->ndo = nd;
#endif

/*
        at this point all of the elements of aa1[lpos]
        have been searched for elements of aa0[tpos]
        with the results in diag[dpos]
*/
   for (nsave = 0, vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++)
   {
      if (vmptr->score > 0)
      {
	 vmptr->score = spam (aa0, aa1, vmptr, ppst->pam2[0], f_str);
	 f_str->vptr[nsave++] = vmptr;
      }
   }

   if (nsave <= 0) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     return;
   }
       
#ifdef FASTX
   /* FASTX code here to modify the start, stop points for 
      the three phases of the translated protein sequence
      */

   /*
     fprintf(stderr,"n0x: %d; n0x31:%d; n0x32: %d\n",n0,n0x31,n0x32);
     for (ib=0; ib<nsave; ib++) {
       fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	       f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
	       f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
	       f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
	       f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
     }

     fprintf(stderr,"---\n");
     */

   for (ib=0; ib<nsave; ib++) {
     if (f_str->noff-f_str->vptr[ib]->dp+f_str->vptr[ib]->start >= n0x32)
       f_str->vptr[ib]->dp += n0x32;
     if (f_str->noff-f_str->vptr[ib]->dp +f_str->vptr[ib]->start >= n0x31)
       f_str->vptr[ib]->dp += n0x31;
   }
	    
   /*
     for (ib=0; ib<nsave; ib++) {
       fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	       f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
	       f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
	       f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
	       f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
     }
     */
#endif /* FASTX */
#ifdef TFASTX
   /* TFASTX code here to modify the start, stop points for 
	     the three phases of the translated protein sequence
	     TFASTX modifies library start points, rather than 
	     query start points
	     */

     /*
   fprintf(stderr,"n0: %d; noff: %d; n1: %d; n1x31: %d n1x32 %d\n",n0, f_str->noff,n1,n1x31,n1x32);
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	     f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
	     f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
	     f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
	     f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
   }

   fprintf(stderr,"---\n");
   */

   for (ib=0; ib<nsave; ib++) {
     if (f_str->vptr[ib]->start >= n1x32) {
       f_str->vptr[ib]->start -= n1x32;
       f_str->vptr[ib]->stop -= n1x32;
       f_str->vptr[ib]->dp -= n1x32;
     }
     if (f_str->vptr[ib]->start >= n1x31) {
       f_str->vptr[ib]->start -= n1x31;
       f_str->vptr[ib]->stop -= n1x31;
       f_str->vptr[ib]->dp -= n1x31;
     }
   }
	    
   /*
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	     f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
	     f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
	     f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
	     f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
   }
   */

#endif /* TFASTX */

   scor = sconn (f_str->vptr, nsave, ppst->param_u.fa.cgap, 
		 ppst->param_u.fa.pgap, f_str);

   for (vmptr=f_str->vptr[0],ib=1; ib<nsave; ib++)
     if (f_str->vptr[ib]->score > vmptr->score) vmptr=f_str->vptr[ib];

/*  kssort (f_str->vptr, nsave); */

   rst->score[1] = vmptr->score;
   rst->score[0] = max (scor, vmptr->score);
   rst->score[2] = rst->score[0];		/* initn */

   if (ppst->param_u.fa.optflag) {
     if (rst->score[0] > ppst->param_u.fa.optcut) {
#ifdef FASTX
       rst->score[2] = dmatchx(aa0, n0,aa1,n1,*hoff=f_str->noff - vmptr->dp,
			     ppst->param_u.fa.optwid, ppst->pam2[0],
			     ppst->gdelval,ppst->ggapval,ppst->gshift,f_str);
#else /* TFASTX */
     /* generate f_str->aa1x */
/*
     for (i=0; i<n1; i++) {
       fputc(ppst->sq[aa1[i]],stderr);
       if (i%60==59) fputc('\n',stderr);
     }
     fprintf(stderr,"\n-----\n");
*/
/*
     fprintf(stderr,"n1: %d, aa1x[n1]: %d; EOSEQ: %d\n",
	     n1,f_str->aa1x[n1],EOSEQ);
     for (fs=aa1,itemp=0; itemp <3; itemp++,fs++) {
       for (fd= &f_str->aa1x[itemp]; *fs!=EOSEQ; fd += 3, fs++) *fd = *fs;
       fprintf(stderr,"fs stopped at: %d\n",(int)(fs-f_str->aa1x));
       *fd=EOSEQ;
     }
*/
/*
     for (i=0; i<n1; i++) {
       fputc(ppst->sq[f_str->aa1x[i]],stderr);
       if (i%60==59) fputc('\n',stderr);
     }
*/
     rst->score[2] = dmatchx(aa0, n0, aa1, n1, *hoff=vmptr->dp-f_str->noff,
			     ppst->param_u.fa.optwid, ppst->pam2[0],
			     ppst->gdelval,ppst->ggapval,ppst->gshift,f_str);
#endif /* TFASTX */
     }
   }
}

void do_work (unsigned char *aa0, int n0,
	      unsigned char *aa1, int n1,
	      int frame,
	      struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, struct rstruct *rst)
{
  int hoff;
  int last_n1, itx, dnav, n10, i, ir;
  unsigned char *aa1x;

  rst->escore = 1.0;
  rst->segnum = rst->seglen = 1;

  if (n1 < ppst->param_u.fa.ktup) {
    rst->score[0] = rst->score[1] = rst->score[2] = 0;
    return;
  }

#ifdef FASTX
  do_fasta (f_str->aa0x, n0, aa1, n1, ppst, f_str, rst, &hoff);
#else
   /* make a precomputed codon number series */

  if (frame == 0) {
    dnav = (hnt[aa1[0]]<<2) + hnt[aa1[1]];
    for (i=2; i<n1; i++) {
      dnav = ((dnav<<2)+hnt[aa1[i]])&255;
      f_str->aa1v[i-2]=dnav;
    }
  }
  else {
    dnav = (3-hnt[aa1[n1-1]]<<2) + 3-hnt[aa1[n1-2]];
    for (i=2, ir=n1-3; i<n1; i++,ir--) {
      dnav = ((dnav<<2)+3-hnt[aa1[ir]])&255;
      f_str->aa1v[i-2]=dnav;
    }
  }

  /* make translated sequence */
  last_n1 = 0;
  aa1x = f_str->aa1x;
  for (itx= frame*3; itx< frame*3+3; itx++) {
    n10  = saatran(aa1,&aa1x[last_n1],n1,itx);
    /*
    fprintf(stderr," itt %d frame: %d\n",itx,frame);
    for (i=0; i<n10; i++) {
      fprintf(stderr,"%c",aa[f_str->aa1x[last_n1+i]]);
      if ((i%60)==59) fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"n10: %d aa1x[] %d last_n1: %d\n",n10,aa1x[last_n1+n10],
	    last_n1);
    */
    last_n1 += n10+1;
  }
  n10 = last_n1-1;

  do_fasta (aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff);
#endif
}

void do_opt (unsigned char *aa0, int n0,
	     unsigned char *aa1, int n1,
	     int frame,
	     struct pstruct *ppst,
	     struct f_struct *f_str,
	     struct rstruct *rst)
{
  int optflag, tscore, hoff;

  optflag = ppst->param_u.fa.optflag;
  ppst->param_u.fa.optflag = 1;

#ifdef FASTX
  do_fasta (f_str->aa0x, n0, aa1, n1, ppst, f_str, rst, &hoff);
#else
  do_fasta (aa0, n0, aa1, n1, ppst, f_str, rst, &hoff);
#endif

  ppst->param_u.fa.optflag = optflag;
}

#ifdef ALLOCN0
void
savemax (dptr, dpos, f_str)
  register struct dstruct *dptr;
  int  dpos;
  struct f_struct *f_str;
{
   register struct savestr *vmptr;
   register int i;

#else
void
savemax (dptr, f_str)
  register struct dstruct *dptr;
  struct f_struct *f_str;
{
   register int dpos;
   register struct savestr *vmptr;
   register int i;

   dpos = (int) (dptr - f_str->diag);

#endif

/* check to see if this is the continuation of a run that is already saved */

   if ((vmptr = dptr->dmax) != NULL && vmptr->dp == dpos &&
	 vmptr->start == dptr->start)
   {
      vmptr->stop = dptr->stop;
      if ((i = dptr->score) <= vmptr->score)
	 return;
      vmptr->score = i;
      if (vmptr != f_str->lowmax)
	 return;
   }
   else
   {
      i = f_str->lowmax->score = dptr->score;
      f_str->lowmax->dp = dpos;
      f_str->lowmax->start = dptr->start;
      f_str->lowmax->stop = dptr->stop;
      dptr->dmax = f_str->lowmax;
   }

   for (vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++)
      if (vmptr->score < i)
      {
	 i = vmptr->score;
	 f_str->lowmax = vmptr;
      }
   f_str->lowscor = i;
}

int spam (unsigned char *aa0, unsigned char *aa1,
	  struct savestr *dmax, int **pam2,
	  struct f_struct *f_str)
{
   int     lpos;
   int     tot, mtot;
   struct {
     int     start, stop, score;
   } curv, maxv;
   register unsigned char *aa0p, *aa1p;

   aa1p = &aa1[lpos = dmax->start];
   aa0p = &aa0[lpos - dmax->dp + f_str->noff];
   curv.start = lpos;

   tot = curv.score = maxv.score = 0;
   for (; lpos <= dmax->stop; lpos++) {
     tot += pam2[*aa0p++][*aa1p++];
     if (tot > curv.score) {
       curv.stop = lpos;
       curv.score = tot;
      }
      else if (tot < 0) {
	if (curv.score > maxv.score) {
	  maxv.start = curv.start;
	  maxv.stop = curv.stop;
	  maxv.score = curv.score;
	}
	tot = curv.score = 0;
	curv.start = lpos+1;
      }
   }

   if (curv.score > maxv.score) {
     maxv.start = curv.start;
     maxv.stop = curv.stop;
     maxv.score = curv.score;
   }

/*	if (maxv.start != dmax->start || maxv.stop != dmax->stop)
		printf(" new region: %3d %3d %3d %3d\n",maxv.start,
			dmax->start,maxv.stop,dmax->stop);
*/
   dmax->start = maxv.start;
   dmax->stop = maxv.stop;

   return maxv.score;
}

#if defined(FASTX) || defined(TFASTX)
#define XFACT 10
#else
#define XFACT 0
#endif

int sconn (struct savestr **v, int n, 
       int cgap, int pgap, struct f_struct *f_str)
{
   int     i, si, cmpp ();
   struct slink
   {
      int     score;
      struct savestr *vp;
      struct slink *next;
   }      *start, *sl, *sj, *so, sarr[MAXSAV];
   int     lstart, tstart, plstop, ptstop;

/*	sort the score left to right in lib pos */

   kpsort (v, n);

   start = NULL;

/*	for the remaining runs, see if they fit */

   for (i = 0, si = 0; i < n; i++)
   {

/*	if the score is less than the gap penalty, it never helps */
      if (v[i]->score < cgap)
	 continue;
      lstart = v[i]->start;
      tstart = lstart - v[i]->dp + f_str->noff;

/*	put the run in the group */
      sarr[si].vp = v[i];
      sarr[si].score = v[i]->score;
      sarr[si].next = NULL;

/* 	if it fits, then increase the score */
      for (sl = start; sl != NULL; sl = sl->next)
      {
	 plstop = sl->vp->stop;
	 ptstop = plstop - sl->vp->dp + f_str->noff;
	 if (plstop < lstart+XFACT && ptstop < tstart+XFACT) {
	   sarr[si].score = sl->score + v[i]->score + pgap;
	   break;
	 }
      }

/*	now recalculate where the score fits */
      if (start == NULL)
	 start = &sarr[si];
      else
	 for (sj = start, so = NULL; sj != NULL; sj = sj->next)
	 {
	    if (sarr[si].score > sj->score)
	    {
	       sarr[si].next = sj;
	       if (so != NULL)
		  so->next = &sarr[si];
	       else
		  start = &sarr[si];
	       break;
	    }
	    so = sj;
	 }
      si++;
   }

   if (start != NULL)
      return (start->score);
   else
      return (0);
}

void
kssort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->score >= v[j + gap]->score)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

void
kpsort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->start <= v[j + gap]->start)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

static int dmatchx(unsigned char *aa0, int n0,
		   unsigned char *aa1, int n1,
		   int hoff, int window, 
		   int **pam2, int gdelval, int ggapval, int gshift,
		   struct f_struct *f_str)
{

   hoff -= window/2;

#ifdef FASTX
   return lx_band(aa1,n1,f_str->aa0v,n0-2, 
		  pam2,-(gdelval-ggapval),-ggapval,-gshift,
		  hoff,window,f_str);
#endif
#ifdef TFASTX
   return lx_band(aa0,n0,f_str->aa1v,n1-2, 
		  pam2,-(gdelval-ggapval),-ggapval,-gshift,
		  hoff,window,f_str);
#endif
}

static void init_row(row, sp)
    struct sx_s *row;
    int sp;
{
  int i;
  for (i = 0; i < sp; i++) {
      row[i].C1 = row[i].I1 = 0;
      row[i].C2 = row[i].I2 = 0;
      row[i].C3 = row[i].I3 = 0;
      row[i].flag = 0;
  }
}

int lx_band(unsigned char *prot_seq,  /* array with protein sequence numbers*/
	    int len_prot,    /* length of prot. seq */
	    unsigned char *dna_prot_seq, /* translated DNA sequence numbers*/
	    int len_dna_prot,   /* length trans. seq. */
	    int **pam_matrix,   /* scoring matrix */
	    int gopen, int gext, /* gap open, gap extend penalties */
	    int gshift,         /* frame-shift penalty */
	    int start_diag,     /* start diagonal of band */
	    int width,         /* width for band alignment */
	    struct f_struct *f_str)
{
  void *ckalloc();
  int i, j, bd, bd1, x1, x2, sp, p1=0, p2=0;
  struct sx_s *last, *tmp;
  int sc, del, best = 0, cd,ci, e1, e2, e3, cd1, cd2, cd3, f, gg;
  register unsigned char *dp;
  register struct sx_s *ap, *aq;
  struct wgt *wt, *ww;
  int aa, b, a,x,y,z;

  sp = width+7;
  gg = gopen+gext;
  /* sp = sp/3+1; */

  if (f_str->cur == NULL) {
    f_str->cur = (struct sx_s *) ckalloc(sizeof(struct sx_s)*sp);
  }

  init_row(f_str->cur, sp);

  /*
  if (start_diag %3 !=0) start_diag = start_diag/3-1;
  else start_diag = start_diag/3;
  if (width % 3 != 0) width = width/3+1;
  else width = width /3;
  */

  x1 = start_diag; 		/* x1 = lower bound of DNA */
  x2 = 1;               /* the amount of position shift from last row*/

  /* i counts through protein sequence, x1 through DNAp */

  for (i = max(0, -width-start_diag), x1+=i; i < len_prot; i++, x1++) {
      bd = min(x1+width, (len_dna_prot+2)/3);	/* upper bound of band */
      bd1 = max(0,x1);	                /* lower bound of band */
      wt = f_str->weight0[prot_seq[i]];
      del = 1-x1;   /*adjustment*/
      bd += del; 
      bd1 +=del;

      ap = &f_str->cur[bd1]; aq = ap+1;
      e1 = f_str->cur[bd1-1].C3; e2 = ap->C1; cd1 = cd2= cd3= 0;
      for (dp = &dna_prot_seq[(bd1-del)*3]; ap < &f_str->cur[bd]; ap++) {
	  ww = &wt[(unsigned char) *dp++];
	  sc = max(max(e1+ww->iv, (e3=ap->C2)+ww->ii), e2+ww->iii);
	  if (cd1 > sc) sc = cd1;
	  cd1 -= gext;
	  if ((ci = aq->I1) > 0) {
	      if (sc < ci) { ap->C1 = ci; ap->I1 = ci-gext;}
	      else {
		  ap->C1 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd1 < sc) cd1 = sc;
		      ap->I1 = max(ci-gext, sc);
		  } else ap->I1 = ci-gext;
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I1 = ap->C1 = 0;
	      } else {
		  ap->C1 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd1 < sc) cd1 = sc;
		      ap->I1 = sc;
		  } else ap->I1 = 0;
	      }
	  }
	  ww = &wt[(unsigned char) *dp++];
	  sc = max(max(e2+ww->iv, (e1=ap->C3)+ww->ii), e3+ww->iii);
	  if (cd2 > sc) sc = cd2;
	  cd2 -= gext;
	  if ((ci = aq->I2) > 0) {
	      if (sc < ci) { ap->C2 = ci; ap->I2 = ci-gext;}
	      else {
		  ap->C2 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd2 < sc) cd2 = sc;
		      ap->I2 = max(ci-gext, sc);
		  }
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I2 = ap->C2 = 0;
	      } else {
		  ap->C2 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd2 < sc) cd2 = sc;
		      ap->I2 = sc;
		  } else ap->I2 = 0;
	      }
	  }
	  ww = &wt[(unsigned char)*dp++];
	  sc = max(max(e3+ww->iv, (e2=aq->C1)+ww->ii), e1+ww->iii);
	  if (cd3 > sc) sc = cd3;
	  cd3 -= gext;
	  if ((ci = aq++->I3) > 0) {
	      if (sc < ci) { ap->C3 = ci; ap->I3 = ci-gext;}
	      else {
		  ap->C3 = sc;
		  sc -= gg;
		  if (sc > 0) {
		      if (sc > best) best =sc;
		      if (cd3 < sc) cd3 = sc;
		      ap->I3 = max(ci-gext, sc);
		  }
	      }
	  } else {
	      if (sc <= 0) {
		  ap->I3 = ap->C3 = 0;
	      } else {
		  ap->C3 = sc; sc-=gg;
		  if (sc >0) {
		      if (sc > best) best =sc;
		      if (cd3 < sc) cd3 = sc;
		      ap->I3 = sc;
		  } else ap->I3 = 0;
	      }
	  }
      }
  }
  /*  printf("The best score is %d\n", best); */
  return best+gg;
}

/* ckalloc - allocate space; check for success */
void *ckalloc(size_t amount)
{
	void *p;

	if ((p = (void *)malloc( (size_t)amount)) == NULL)
		w_abort("Ran out of memory.","");
	return(p);
}

/* calculate the 100% identical score */
shscore(unsigned char *aa0, int n0, int **pam2)
{
  int i, sum;
  for (i=0,sum=0; i<n0; i++)
    sum += pam2[aa0[i]][aa0[i]];
  return sum;
}

#define SGW1 100
#define SGW2 300
#define WIDTH 60

typedef struct mat *match_ptr;

typedef struct mat {
	int i, j, l;
	match_ptr next;
} match_node;

typedef struct {
	int i,j;
} state;

typedef state *state_ptr;

typedef struct st_s { int C, I, D;} *st_ptr;

static st_ptr up=NULL, down, tp;

static int *st_up;

static int gop, gext, shift;

void *ckalloc();
static match_ptr small_global(), global();
static int local_align(), find_best();
static void init_row2(), init_ROW();

extern int pro_dna(unsigned char *prot_seq,  /* array with prot. seq. numbers*/
		   int len_prot,    /* length of prot. seq */
		   unsigned char *dna_prot_seq, /* trans. DNA seq. numbers*/
		   int len_dna_prot,   /* length trans. seq. */
		   int **pam_matrix,   /* scoring matrix */
		   int gopen, int gex, /* gap open, gap extend penalties */
		   int gshift,         /* frame-shift penalty */
		   struct f_struct *f_str,
		   int *alignment,  /*store the alignment*/
		   int max_res,
		   int *nres,		/* length of alignment */
		   struct a_struct *aln) /* alignment info */
{
	match_ptr align, ap, aq;
	int x, y, ex, ey, i, score;

	gext = gex; gop = gopen; shift = gshift;

	if (up==NULL) {
	  up = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+10));
	  down = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+10));
	  tp = (st_ptr) ckalloc(sizeof(struct st_s)*(len_dna_prot+10));
	  st_up = (int *) ckalloc(sizeof(int)*(len_dna_prot+10));
	}

	/*local alignment find the best local alignment x and y
	  is the starting position of the best local alignment
	  and ex ey is the ending position */

	score= local_align(&x, &y, &ex, &ey, pam_matrix,
			   dna_prot_seq, len_dna_prot,
			   prot_seq, len_prot, f_str);

	up += 3; down += 3; tp += 3;
	/* aln->min_n = y; aln->max_n = ey; */

	align = global(x, y, ex, ey, pam_matrix, dna_prot_seq, prot_seq,
		       0, 0, f_str);

	/* x, y - start in prot, dna_prot */
	alignment[0] = x; /* start of alignment in DNA */
	alignment[1] = y; /* start of alignment in prot */
	for (ap = align, i= 2; ap; i++) {
	    if (i < max_res) alignment[i] = ap->l;
	    aq = ap->next; free(ap); ap = aq;
	}
	if (i >= max_res)
	  fprintf(stderr,"***alignment truncated: %d/%d***\n", max_res,i);

	up = &up[-3]; down = &down[-3]; tp = &tp[-3];
	free(up); free(tp); free(down); free(st_up);
	up = NULL;

	*nres = i;
	return score;
}

static void swap(void **a, void **b)
{
    void *t = *a;
    *a = *b;   *b = t;
}

/*
   local alignment find the best local alignment x and y
   is the starting position of the best local alignment
   and ex ey is the ending position 
*/
static local_align(x, y, ex, ey, wgts, dnap, ld, pro, lp, f_str)
     int *x, *y, *ex, *ey, ld, **wgts, lp;
     unsigned char *dnap;
     unsigned char *pro;
     struct f_struct *f_str;
{
	int i, j,  score, x1,x2,x3,x4, e1 = 0, e2 = 0, e3,
	    sc, del,  e, best = 0,  cd, ci, c;
	struct wgt *wt, *ww;
	state_ptr cur_st, last_st, cur_i_st;
	st_ptr cur, last;
	unsigned char *dp;
	int *cur_d_st;

/*      
   Array rowiC stores the best scores of alignment ending at a position
   Arrays rowiD and rowiI store the best scores of alignment ending
                 at a position with a deletion or insrtion
   Arrays sti stores the starting position of the best alignment whose
              score stored in the corresponding row array.
   The program stores two rows to complete the computation, same is
        for the global alignment routine.
*/
	ld += 2;
	init_ROW(up, ld+1);
	init_ROW(down, ld+1);
	init_row2(st_up, ld+3);
	cur = up+1;
	last = down+1; 
	cur_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
	last_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
	cur_i_st = (state_ptr) ckalloc(sizeof(state)*(ld+1));
	cur_d_st = st_up; 
	dp = dnap-2;
	for (i = 0; i < lp; i++) {
	        wt = f_str->weight1[pro[i]];  e2 =0; e1 = last[0].C;
		for (j = 0; j < 2; j++) {
		    cur_st[j].i = i+1;
		    cur_st[j].j = j+1;
		}
		for (j = 2; j < ld; j++) {
			ww = &wt[(unsigned char) dp[j]];
			del = -1;
			if (j >= 3) {
			    sc = 0;
			    e3 = e2; e2 = e1;
			    e1 = last[j-2].C; 
			    if ((e=e2+ww->iii) > sc) {sc = e; del = 3;}
			    if ((e=e1+ww->ii) > sc) {sc = e; del = 2;}
			    if ((e = e3+ww->iv) > sc) {sc = e; del = 4;} 
			} else {
			    sc = e2  = 0;
			    if (ww->iii > 0) {sc = ww->iii; del = 3;}
			}
			if (sc < (ci=last[j].I)) {
			    sc = ci; del = 0;
			}
			if (sc < (cd=cur[j].D)) {
			    sc = cd; del = 5;
			}
			cur[j].C = sc;
			e = sc  - gop;
			if (e > cd) {
			    cur[j+3].D = e-gext;
			    cur_d_st[j+3] = 3;
			} else {
			    cur[j+3].D = cd-gext;
			    cur_d_st[j+3] = cur_d_st[j]+3;
			}
			switch(del) {
			case 5:
			    c = cur_d_st[j];
			    cur_st[j].i = cur_st[j-c].i;
			    cur_st[j].j = cur_st[j-c].j;
			    break;
			case 0:
			    cur_st[j].i = cur_i_st[j].i;
			    cur_st[j].j = cur_i_st[j].j;
			    break;
			case 2:
			case 3:
			case 4:
			    if (i) {
				if (j-del >= 0) {
				    cur_st[j].i = last_st[j-del].i;
				    cur_st[j].j = last_st[j-del].j;
				} else {
				    cur_st[j].i = i;
				    cur_st[j].j = 0;
				}
			    } else {
				cur_st[j].i = 0;
				cur_st[j].j = max(0, j-del+1);
			    }
			    break;
			case -1:
			    cur_st[j].i = i+1;
			    cur_st[j].j = j+1;
			    break;
			}
			if (e > ci) {
			    cur[j].I  = e -gext;
			    cur_i_st[j].i = cur_st[j].i;
			    cur_i_st[j].j = cur_st[j].j;
			} else {
			    cur[j].I  = ci- gext;
			}
			if (sc > best) {
				x1 = cur_st[j].i;
				x2 = cur_st[j].j;
				best =sc;
				x3 = i;
				x4 = j;
			}
		}
		swap((void *)&last, (void *)&cur);
		swap((void *)&cur_st, (void *)&last_st);
	}
	/*	printf("The best score is %d\n", best);*/
	*x = x1; *y = x2; *ex = x3; *ey = x4;
	free(cur_st); free(last_st); free(cur_i_st); 
	return best;
}

/* 
   Both global_up and global_down do linear space score only global 
   alignments on subsequence pro[x]...pro[ex], and dna[y]...dna[ey].
   global_up do the algorithm upwards, from row x towards row y.
   global_down do the algorithm downwards, from row y towards x.
*/

static void global_up(row1, row2, x, y, ex, ey, wgts, dnap, pro, N, f_str)
     st_ptr *row1, *row2;
     int  x, y, ex, ey, **wgts;
     unsigned char *dnap;
     unsigned char  *pro;
     int N;
     struct f_struct *f_str;
{
	int i, j, k, sc, e, e1, e2, e3, t, ci, cd, score;
	struct wgt *wt, *ww;
	st_ptr cur, last;

	cur = *row1; last = *row2;
	sc = -gop;
	for (j = 0; j <= ey-y+1; j++) {
	    if (j % 3 == 0) {last[j].C = sc; sc -= gext; last[j].I = sc-gop;}
	    else { last[j].I = last[j].C = -10000;}
	}  
	last[0].C = 0; cur[0].D = cur[1].D = cur[2].D = -10000;
	last[0].D = last[1].D = last[2].D = -10000;
	if (N) last[0].I = -gext;
	for (i = 1; i <= ex-x+1; i++) {
	        wt = f_str->weight1[pro[i+x-1]]; e1 = -10000; e2 = last[0].C;
		for (j = 0; j <= ey-y+1; j++) {
		    t = j+y;
		    sc = -10000; 
		    ww = &wt[(unsigned char) dnap[t-3]]; 
		    if (j < 4) {
			if (j == 3) {
			    sc = e2+ww->iii;
			} else if (j == 2) {
			    sc = e2 + ww->ii;
			}
		    } else {
			e3 = e2; e2 = e1;
			e1 = last[j-2].C;
			sc = max(e2+ww->iii, max(e1+ww->ii, e3+ww->iv));
		    }
		    sc = max(sc, max(ci=last[j].I, cd = cur[j].D));
		    cur[j].C = sc;
		    cur[j+3].D = max(cd, sc-gop)-gext;
		    cur[j].I = max(ci, sc-gop)-gext;
		}
		swap((void *)&last, (void *)&cur);
	}
	/*printf("global up score =%d\n", last[ey-y+1].C);*/
	for (i = 0; i <= ey-y+1; i++) last[i].I = cur[i].I;	
	if (*row1 != last) swap((void *)row1, (void *)row2);
}

static void global_down(row1, row2,  x, y, ex, ey, wgts, dnap, pro, N, f_str)
     st_ptr *row1, *row2;
     int x, y, ex, ey, **wgts;
     unsigned char *dnap;
     char  *pro;
     int N;
     struct f_struct *f_str;
{
	int i, j, k, sc, del, *tmp, e,  t, e1,e2,e3, ci,cd, score;
	struct wgt *wt, *w1, *w2, *w3;
	st_ptr cur, last;

	cur = (*row1); last = *row2;
	sc = -gop;
	for (j = ey-y+1; j >= 0; j--) {
	    if ((ey-y+1-j) % 3) {last[j].C = sc; sc-=gext; last[j].I = sc-gop;}
	    else  last[j].I =  last[j].C = -10000;
	    cur[j].I = -10000;
	} 
	last[ey-y+1].C = 0;
	if (N) last[ey-y+1].I = -gext;
	cur[ey-y+1].D = cur[ey-y].D = cur[ey-y-1].D = -10000;
	last[ey-y+1].D = last[ey-y].D = last[ey-y-1].D = -10000;
	for (i = ex-x; i >= 0; i--) {
	        wt = f_str->weight1[pro[i+x]]; e2 = last[ey-y+1].C; 
		e1 = -10000;
		w3 = &wt[(unsigned char) dnap[ey]]; 
		w2 = &wt[(unsigned char) dnap[ey-1]];
		for (j = ey-y+1; j >= 0; j--) {
		    t = j+y;
		    w1 = &wt[(unsigned char) dnap[t-1]];
		    sc = -10000;
		    if (t+3 > ey) {
			if (t+2 == ey) {
			    sc = e2+w2->iii; 
			} else if (t+1 == ey) {
			    sc = e2+w1->ii;
			}
		    } else {
			e3 = e2; e2 = e1;
			e1 = last[j+2].C;
			sc = max(e2+w2->iii, max(e1+w1->ii,e3+w3->iv)) ;
		    }
		    if (sc < (cd= cur[j].D)) {
			sc = cd; 
			cur[j-3].D = cd-gext;
		    } else cur[j-3].D =max(cd, sc-gop)-gext;
		    if (sc < (ci= last[j].I)) {
			sc = ci;
			cur[j].I = ci - gext;
		    } else cur[j].I = max(sc-gop,ci)-gext;
		    cur[j].C = sc;
		    w3 = w2; w2 = w1;
		}
		swap((void *)&last, (void *)&cur);
	}
	for (i = 0; i <= ey-y+1; i++) last[i].I = cur[i].I;
	if (*row1 != last) swap((void *)row1, (void *)row2);
}

static void init_row2(row, ld)
int *row, ld;
{
	int i;
	for (i = 0; i < ld; i++) row[i] = 0;
}

static void init_ROW(row, ld)
st_ptr row;
int ld;
{
    int i;
    for (i = 0; i < ld; i++) row[i].I = row[i].D = row[i].C = 0;
}

static match_ptr combine(x1, x2, st)
match_ptr x1, x2;
int st;
{
	match_ptr x;

	if (x1 == NULL) return x2;
	for (x = x1; x->next; x = x->next);
	x->next = x2;
	if (st) {
	    for (x = x2; x; x = x->next) {
		x->j++;
		if (x->l == 3 || x->l == 4) break;
	    }
	    x->l--;
	}
	return x1;
}

/*
   global use the two upwards and downwards score only linear
   space global alignment subroutine to recursively build the
   alignment.
*/

match_ptr global(x,y, ex, ey, wgts, dnap, pro, N1, N2, f_str)
     int x, y, ex, ey, **wgts;
     unsigned char *dnap;
     unsigned char *pro;
     int N1, N2;
     struct f_struct *f_str;
{
	int m;
	int m1, m2;
	match_ptr x1, x2, mm1, mm2;
	/*printf("%d %d %d %d %d %d\n", x,y, ex, ey, N1, N2);*/
/*
   if the space required is limited, we can do a quadratic space
   algorithm to find the alignment.
*/
	if (ex <= x) {
	    mm1  = NULL;
	    for (m = y+3; m <= ey; m+=3) {
		x1 = (match_ptr) ckalloc(sizeof(match_node));
		x1->l = 5; x1->next = mm1; 
		if (mm1== NULL) mm2 = x1;
		mm1 = x1;
	    }
	    if (ex == x) {
		if ((ey-y) % 3 != 0) {
		    x1  = (match_ptr) ckalloc(sizeof(match_node));
		    x1->l = ((ey-y) % 3) +1; x1->next = NULL;
		    if (mm1) mm2->next = x1; else mm1 = x1;
		} else mm2->l = 4;
	    }
	    return mm1;
	}
	if (ey <= y) {
	    mm1  = NULL;
	    for (m = x; m <= ex; m++) {
		x1 = (match_ptr) ckalloc(sizeof(match_node));
		x1->l = 0; x1->next = mm1; mm1 = x1;
	    }
	    return mm1;
	}
	if (ex -x < SGW1 && ey-y < SGW2) 
	    return small_global(x,y,ex,ey,wgts, dnap, pro, N1, N2,f_str);
	m = (x+ex)/2;
/*     
   Do the score only global alignment from row x to row m, m is
   the middle row of x and ex. Store the information of row m in
   upC, upD, and upI.
*/
	global_up(&up, &tp,  x, y, m, ey, wgts, dnap, pro, N1, f_str);
/* 
   Do the score only global alignment downwards from row ex
   to row m+1, store information of row m+1 in downC downI and downD
*/
	global_down(&down, &tp, m+1, y, ex, ey, wgts, dnap, pro, N2, f_str);
/*
   Use these information of row m and m+1, to find the crossing
   point of the best alignment with the middle row. The crossing
   point is given by m1 and m2. Then we recursively call global
   itself to compute alignments in two smaller regions found by
   the crossing point and combine the two alignments to form a
   whole alignment. Return that alignment.
*/
	if (find_best(up, down, &m1, &m2, ey-y+1, y)) {
	    x1 = global(x, y, m, m1, wgts, dnap, pro, N1, 0, f_str);
	    x2 = global(m+1, m2, ex, ey, wgts, dnap, pro, 0, N2, f_str);
	    if (m1 == m2) x1 = combine(x1,x2,1);
	    else x1 = combine(x1, x2,0);
	} else {
	    x1 = global(x, y, m-1, m1, wgts, dnap, pro, N1, 1, f_str);
	    x2 = global(m+2, m2, ex, ey, wgts, dnap, pro, 1, N2, f_str);
	    mm1 = (match_ptr) ckalloc(sizeof(match_node));
	    mm1->i = m; mm1->l = 0; mm1->j = m1;
	    mm2 = (match_ptr) ckalloc(sizeof(match_node));
	    mm2->i = m+1; mm2->l = 0; mm2->j = m1;
	    mm1->next = mm2; mm2->next = x2;
	    x1 = combine(x1, mm1, 0);
	}
	return x1;
}

static find_best(up, down, m1, m2, ld, y)
     st_ptr up, down; 
     int ld, y;
     int *m1, *m2;
{
	int i, best = -1000, j = 0, s1, s2, s3, s4, st;

	for (i = 1; i < ld; i++) {
	    s2 = up[i].C + down[i].C;
	    s4 = up[i].I + down[i].I + gop;
	    if (best < s2) {
		best = s2; j = i; st = 1;
	    }
	    if (best < s4) {
		best = s4; j = i; st = 0;
	    }
	}
	*m1 = j-1+y;
	*m2 = j+y;
	/*printf("score=%d\n", best);*/
	return st;
} 

/*
   An alignment is represented as a linked list whose element
   is of type match_node. Each element represent an edge in the
   path of the alignment graph. The fields of match_node are
   l ---  gives the type of the edge.
   i, j --- give the end position.
*/

static match_ptr small_global(x, y, ex, ey, wgts, dnap, pro, N1, N2, f_str)
     int x, y, ex, ey, **wgts;
     unsigned char *dnap;
     unsigned char *pro;
     int N1, N2;
     struct f_struct *f_str;
{
        static int C[SGW1+1][SGW2+1], st[SGW1+1][SGW2+1], D[SGW2+7], I[SGW2+1];
        int i, j, e, sc, score, del, k, t,  ci, cd;
	int *cI, *cD, *cC, *lC, *cst, e2, e3, e4;
	match_ptr mp, first;
	struct wgt *wt, *ww;

	/*printf("small_global %d %d %d %d\n", x, y, ex, ey);*/
	sc = -gop-gext; C[0][0] = 0; 
	if (N1) I[0] = -gext; else I[0] = sc;
	for (j = 1; j <= ey-y+1; j++) {
	    if (j % 3== 0) {
		C[0][j] = sc; sc -= gext; I[j] = sc-gop;
	    } else I[j] = C[0][j] = -10000;
	    st[0][j] = 5;
	}
	lC = &C[0][0]; cD = D; D[0] = D[1] = D[2] = -10000;
	cI = I;
	for (i = 1; i <= ex-x+1; i++) {
	    cC = &C[i][0];	
	    wt = f_str->weight1[pro[i+x-1]]; cst = &st[i][0];
	    for (j = 0; j <=ey-y+1; j++) {
		ci = cI[j];
		cd= cD[j];
		t = j+y;
		ww = &wt[(unsigned char) dnap[t-3]];
		if (j >= 4) {
		    sc = lC[j-3]+ww->iii; e2 = lC[j-2]+ww->ii;  
		    e4 = lC[j-4]+ww->iv; del = 3;
		    if (e2 > sc) { sc = e2; del = 2;}
		    if (e4 >= sc) { sc = e4; del = 4;}
		} else {
		    if (j == 3) {
			sc = lC[0]+ww->iii; del =3;
		    } else if (j == 2) {
			sc = lC[0]+ww->ii; del = 2;
		    } else {sc = -10000; del = 0;}
		}
		if (sc < ci) {
		    sc = ci; del = 0; 
		}
		if (sc <= cd) {
		    sc = cd;
		    del = 5;
		}
		cC[j] = sc;
		sc -= gop;
		if (sc <= cd) {
		    del += 10;
		    cD[j+3] = cd - gext;
		} else cD[j+3] = sc -gext;
		if (sc < ci) {
		    del += 20;
		    cI[j] = ci-gext;
		} else cI[j] = sc-gext;
		*(cst++) = del;
	    }
	    lC = cC;
	}
	/*printf("small global score =%d\n", C[ex-x+1][ey-y+1]);*/
	if (N2 && cC[ey-y+1] <  ci+gop) st[ex-x+1][ey-y+1] =0;
	first = NULL; e = 1;
	for (i = ex+1, j = ey+1; i > x || j > y; i--) {
		mp = (match_ptr) ckalloc(sizeof(match_node));
		mp->i = i-1;
		k  = (t=st[i-x][j-y])%10;
		mp->j = j-1;
		if (e == 5 && (t/10)%2 == 1) k = 5;
		if (e == 0 && (t/20)== 1) k = 0;
		if (k == 5) { j -= 3; i++; e=5;}
		else {j -= k;if (k==0) e= 0; else e = 1;}
		mp->l = k;
		mp->next = first;
		first = mp;
	}

/*	for (i = 0; i <= ex-x; i++) {
		for (j = 0; j <= ey-y; j++) 
			printf("%d ", C[i][j]);
		printf("\n");
	}
*/
	return first;	
}

#define XTERNAL
#include "upam.h"

void
display_alig(a, dna, pro,length, ld, f_str)
     int *a;
     unsigned char *dna, *pro;
     int length, ld;
     struct f_struct *f_str;
{
	int len = 0, i, j, x, y, lines, k, iaa;
	static char line1[100], line2[100], line3[100],
		 tmp[10] = "         ", *st;
	char *dna1, c1, c2, c3;

	line1[0] = line2[0] = line3[0] = '\0'; x= a[0]; y = a[1]-3;

	printf("\n%5d\n%5d", y+3, x);
	for (len = 0, j = 2, lines = 0; j < length; j++) {
		i = a[j];
		line3[len] = ' ';
		switch (i) {
		case 3: 
		    y += 3;
		    line2[len] = aa[iaa=pro[x++]];
		    line1[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c5;
		    if (line1[len] != f_str->weight_c[iaa][(unsigned char) dna[y]].c3)
			line3[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c3;
		    break;
		case 2:
		    y += 2;
		    line1[len] = '\\';
		    line2[len++] = ' ';
		    line2[len] = aa[iaa=pro[x++]];
		    line1[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c2;
		    line3[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c3;
		    break;
		case 4:
		    y += 4;
		    line1[len] = '/';
		    line2[len++] = ' ';
		    line2[len] = aa[iaa=pro[x++]];
		    line1[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c4;
		    line3[len] = f_str->weight_c[iaa][(unsigned char) dna[y]].c3;
		    break;
		case 5:
		    y += 3;
		    line1[len] = f_str->weight_c[0][(unsigned char) dna[y]].c3;
		    line2[len] = '-';
		    break;
		case 0:
		    line1[len] = '-';
		    line2[len] = aa[pro[x++]];
		    break;
		}
		len++;
		line1[len] = line2[len]  = line3[len]  = '\0'; 
		if (len >= WIDTH) {
		    for (k = 10; k <= WIDTH; k+=10) 
			printf("    .    :");
		    if (k-5 < WIDTH) printf("    .");
		    c1 = line1[WIDTH]; c2 = line2[WIDTH]; c3 = line3[WIDTH];
		    line1[WIDTH] = line2[WIDTH] = line3[WIDTH] = '\0';
		    printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
		    line1[WIDTH] = c1; line2[WIDTH] = c2;
		    strcpy(line1, &line1[WIDTH]);
		    strcpy(line2, &line2[WIDTH]);
		    strcpy(line3, &line3[WIDTH]);
		    len = len - WIDTH;
		    printf("\n%5d\n%5d", y+3, x);
		}
        }
	for (k = 10; k < len; k+=10) 
	    printf("    .    :");
	if (k-5 < len) printf("    .");
	printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
}


/* alignment store the operation that align the protein and dna sequence.
   The code of the number in the array is as follows:
   0:     delete of an amino acid.
   2:     frame shift, 2 nucleotides match with an amino acid
   3:     match an  amino acid with a codon
   4:     the other type of frame shift
   5:     delete of a codon
   

   Also the first two element of the array stores the starting point 
   in the protein and dna sequences in the local alignment.

   Display looks like where WIDTH is assumed to be divisible by 10.

    0    .    :    .    :    .    :    .    :    .    :    .    :
     AACE/N\PLK\G\HK\Y/LWA\S\C\E/P\PRIRZ/G\HK\Y/LWA\S\C\E/P\PRIRZ
          I S   G S  V F   N R Q L A     G S  V F   N R Q L A    
     AACE P P-- G HK Y TWA A C E P P---- G HK Y TWA A C E P P----

   60    .    :    .    :    .    :    .    :    .    :    .    :
     /G\HK\Y/LWA\S\C\E/P\PRIRZ/G\HK\Y/LWA\S\C\E/P\PRIRZ/G\HK\Y/LW
      G S  V F   N R Q L A     G S  V F   N R Q L A     G S  V F 
      G HK Y TWA A C E P P---- G HK Y TWA A C E P P---- G HK Y TW

For frame shift, the middle row show the letter in the original sequence,
and the letter in the top row is the amino acid that is chose by the 
alignment (translated codon chosen from 4 nucleotides, or 2+1).
*/

/* fatal - print message and die */
void
fatal(msg)
char *msg;
{
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

int do_walign (unsigned char *aa0, int n0,
	       unsigned char *aa1, int n1,
	       int frame,
	       struct pstruct *ppst, 
	       struct f_struct *f_str, 
	       int **ares, int *nres, struct a_struct *aln)
{
  int score;
  int i, ir, last_n1, itemp, n10, itx, dnav;
  unsigned char *aa1x;
  
#ifdef FASTX
  score = pro_dna(aa1, n1, f_str->aa0v, n0-2, ppst->pam2[0],
		 -(ppst->gdelval - ppst->ggapval), -ppst->ggapval,
		 -ppst->gshift,
		 f_str, f_str->res, f_str->max_res, nres, aln);
  /*    display_alig(f_str->res,f_str->aa0v+2,aa1,*nres,n0-2,f_str); */
  aln->llrev = 0;
  aln->llfact = 1;
  aln->llmult = 1;
  aln->qlfact = 3;
  aln->frame = 0;
  if (frame > 0) aln->qlrev = 1;
  else aln->qlrev = 0;
#endif
#ifdef TFASTX
   /* make a precomputed codon number series */
  if (frame==0) {
    dnav = (hnt[aa1[0]]<<2) + hnt[aa1[1]];
    for (i=2; i<n1; i++) {
      dnav = ((dnav<<2)+hnt[aa1[i]])&255;
      f_str->aa1v[i-2]=dnav;
    }
  }
  else { /* must do things backwards */
    dnav = (3-hnt[aa1[n1-1]]<<2) + 3-hnt[aa1[n1-2]];
    for (i=2,ir=n1-3; i<n1; i++,ir--) {
      dnav = ((dnav<<2)+3-hnt[aa1[ir]])&255;
      f_str->aa1v[i-2]=dnav;
    }
  }

  /* make translated sequence */
  last_n1 = 0;
  aa1x = f_str->aa1x;
  for (itx= frame*3; itx< frame*3+3; itx++) {
    n10  = saatran(aa1,&aa1x[last_n1],n1,itx);
    /*
    fprintf(stderr," itt %d itx: %d\n",itt,itx);
    for (i=0; i<n10; i++) {
      fprintf(stderr,"%c",aa[f_str->aa1x[last_n1+i]]);
      if ((i%60)==59) fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
    */
    last_n1 += n10+1;
  }
  n10 = last_n1-1;

  score = pro_dna(aa0, n0, f_str->aa1v, n1-2, ppst->pam2[0],
		 -(ppst->gdelval - ppst->ggapval), -ppst->ggapval,
		 -ppst->gshift,
		 f_str, f_str->res, f_str->max_res, nres, aln);
  /*   display_alig(f_str->res,f_str->aa0y,aa1,*nres,n0,f_str); */
  aln->qlfact = 1;
  aln->qlrev = 0;
  aln->llfact = 3;
  aln->llmult = 1;
  aln->frame = 0;
  if (frame > 0) aln->llrev = 1;
  else aln->llrev = 0;
#endif


  *ares = f_str->res;

  return score;
}


#include "structs.h"

int calcons(unsigned char *aa0, int n0,
	    unsigned char *aa1, int n1,
	    int *res, int nres, int *nc,
	    struct a_struct *aln, struct pstruct pst,
	    char *seqc0, char *seqc1,
	    struct f_struct *f_str)
{
  int i0, i1;
  int lenc, not_c, itmp;
  char *sp0, *sp1, *sq;
  unsigned char aap, *ap0, *ap1;
  int *rp, *rpmax;
  int ngap_p, ngap_d, nfs;
  
  /* don't fill in the ends */

  if (pst.ext_sq_set) {sq = pst.sqx;}
  else {sq = pst.sq;}

  /* res[0] has start of protein sequence */
  /* res[1] has start of translated DNA sequence */

#ifdef FASTX
  ap0 = f_str->aa0v;		/* computed codons -> ap0*/
  ap1 = aa1;			/* protein sequence -> ap1 */
  aln->smin1 = aln->min1= *res++;	/* start in protein sequence */
  aln->smin0= *res++;		/* start in DNA/codon sequence */
  aln->min0 = aln->smin0-3;	/* codon_start - 3 */
#else	/* TFASTXYZ */
  ap0 = f_str->aa1v;		/* computed codons -> ap0*/
  ap1 = aa0;			/* protein sequence */
  aln->smin0 = aln->min1 = *res++;	/* start in protein sequence */
  aln->smin1 = *res++;		/* start in codon sequence */
  aln->min0 = aln->smin1-3;	/* codon_start - 3 */
#endif

  rp = res;			/* start of alignment info */
  rpmax = &res[nres-2];		/* end of alignment info */

  aln->smins = aln->mins = 0;

/* now get the middle */
#ifdef FASTX
  sp0 = seqc0;		/* sp0/seqc0 is codon sequence */
  sp1 = seqc1;		/* sp1/seqc1 is protein sequence */
#endif
#ifdef TFASTX
  sp1 = seqc0;		/* sp1/seqc0 is protein sequence */
  sp0 = seqc1;		/* sp0/seqc1 is codon sequence */
#endif

  lenc = not_c = aln->nident = ngap_p = ngap_q = nfs = 0;
  i0 = aln->min0;	/* start of codon sequence */
  i1 = aln->min1;	/* start of protein sequence */

  while (rp < rpmax ) {
    switch (*rp++) {
    case 3:
      i0 += 3;
      *sp1 = sq[aap=ap1[i1++]];
      *sp0 = f_str->weight_c[aap][ap0[i0]].c5;
      if (toupper(*sp0) == toupper(*sp1)) aln->nident++;
      sp0++; sp1++;
      lenc++;
      break;
    case 2:
      nfs++;
      i0 += 2;
      *sp0++ = '/';
      *sp1++ = '-';
      not_c++;
      *sp1 = sq[aap=ap1[i1++]];
      *sp0 = f_str->weight_c[aap][ap0[i0]].c2;
      if (toupper(*sp0) == toupper(*sp1)) aln->nident++;
      sp0++; sp1++;
      lenc++;
      break;
    case 4:
      nfs++;
      i0 += 4;
      *sp0++ = '\\';
      *sp1++ = '-';
      not_c++;
      *sp1 = sq[aap=ap1[i1++]];
      *sp0 = f_str->weight_c[aap][ap0[i0]].c4;
      if (toupper(*sp0) == toupper(*sp1)) aln->nident++;
      sp0++; sp1++;
      lenc++;
      break;
    case 5:
      i0 += 3;
      *sp0++ = f_str->weight_c[0][ap0[i0]].c3;
      *sp1++ = '-';
      lenc++;
      ngap_p++;
      break;
    case 0:
      *sp0++ = '-';
      *sp1++ = sq[ap1[i1++]];
      lenc++;
      ngap_d++;
      break;
    }
  }

#ifdef FASTX
  aln->max0 = i0+3;	/* end of codon sequence */
  aln->max1 = i1;	/* end of protein sequence */
  aln->ngap_q = ngap_d;
  aln->ngap_l = ngap_p;
  aln->nfs = nfs;
#endif
#ifdef TFASTX
  aln->max1 = i0+3;	/* end of codon sequence */
  aln->max0 = i1;	/* end of protein sequence */
  aln->ngap_q = ngap_p;
  aln->ngap_l = ngap_d;
  aln->nfs = nfs;
#endif
  aln->min0 = aln->smin0;
  aln->min1 = aln->smin1;

  if (lenc < 0) lenc = 1;

  *nc = lenc;
/*	now we have the middle, get the right end */

  return lenc+not_c;
}

int calc_id(const unsigned char *aa0, const int n0,
	    const unsigned char *aa1, const int n1,
	    int *res, int nres,
	    struct a_struct *aln, struct pstruct pst,
	    struct f_struct *f_str)
{
  int i0, i1;
  int lenc, not_c, itmp;
  char sp0, sp1;
  unsigned char aap;
  const unsigned char *ap0, *ap1;
  int *rp, *rpmax;
  
  /* don't fill in the ends */

  /* res[0] has start of protein sequence */
  /* res[1] has start of translated DNA sequence */

#ifdef FASTX
  ap0 = f_str->aa0v;		/* computed codons -> ap0*/
  ap1 = aa1;			/* protein sequence -> ap1 */
  aln->smin1 = aln->min1= *res++;	/* start in protein sequence */
  aln->smin0= *res++;		/* start in DNA/codon sequence */
  aln->min0 = aln->smin0-3;	/* codon_start - 3 */
#else	/* TFASTXYZ */
  ap0 = f_str->aa1v;		/* computed codons -> ap0*/
  ap1 = aa0;			/* protein sequence */
  aln->smin0 = aln->min1 = *res++;	/* start in protein sequence */
  aln->smin1 = *res++;		/* start in codon sequence */
  aln->min0 = aln->smin1-3;	/* codon_start - 3 */
#endif

  rp = res;			/* start of alignment info */
  rpmax = &res[nres-2];		/* end of alignment info */

  aln->smins = aln->mins = 0;

/* now get the middle */

  lenc = not_c = aln->nident = ngap_p = ngap_d = 0;
  i0 = aln->min0;	/* start of codon sequence */
  i1 = aln->min1;	/* start of protein sequence */

  while (rp < rpmax ) {
    switch (*rp++) {
    case 3:
      i0 += 3;
      sp1 = pst.sq[aap=ap1[i1++]];
      sp0 = f_str->weight_c[aap][ap0[i0]].c5;
      if (sp0 == sp1) aln->nident++;
      lenc++;
      break;
    case 2:
      i0 += 2;
      not_c++;
      sp1 = pst.sq[aap=ap1[i1++]];
      sp0 = f_str->weight_c[aap][ap0[i0]].c2;
      if (sp0 == sp1) aln->nident++;
      lenc++;
      break;
    case 4:
      i0 += 4;
      not_c++;
      sp1 = pst.sq[aap=ap1[i1++]];
      sp0 = f_str->weight_c[aap][ap0[i0]].c4;
      if (sp0 == sp1) aln->nident++;
      lenc++;
      break;
    case 5:
      i0 += 3;
      lenc++;
      ngap_p++;
      break;
    case 0:
      lenc++;
      ngap_d++;
      break;
    }
  }

#ifdef FASTX
  aln->max0 = i0+3;	/* end of codon sequence */
  aln->max1 = i1;	/* end of protein sequence */
  aln->ngap_q = ngap_d;
  aln->ngap_l = ngap_p;
#endif
#ifdef TFASTX
  aln->max1 = i0+3;	/* end of codon sequence */
  aln->max0 = i1;	/* end of protein sequence */
  aln->ngap_q = ngap_p;
  aln->ngap_l = ngap_d;
#endif
  aln->min0 = aln->smin0;
  aln->min1 = aln->smin1;

  if (lenc < 0) lenc = 1;

/*	now we have the middle, get the right end */

  return lenc+not_c;
}

#ifdef PCOMPLIB
#include "p_mw.h"
void
update_params(struct qmng_str *qm_msg, struct pstruct *ppst)
{
  ppst->n0 = qm_msg->n0;
}
#endif
