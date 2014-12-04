
/* copyright (c) 1998, 1999 William R. Pearson and the U. of Virginia */

/* $Name: fa34t20b3 $ - $Id: dropffa.c,v 1.12 2002/08/28 21:09:17 wrp Exp $ */

/* this code implements the "fastf" algorithm, which is designed to
   deconvolve mixtures of protein sequences derived from mixed-peptide
   Edman sequencing.  The expected input is:

   >test | 40001 90043 | mgstm1
   MGCEN,
   MIDYP,
   MLLAY,
   MLLGY

   Where the ','s indicate the length/end of the sequencing cycle
   data.  Thus, in this example, the sequence is from a mixture of 4
   peptides, M was found in the first position, G,I, and L(2) at the second,
   C,D, L(2) at the third, etc.

   Because the sequences are derived from mixtures, there need not be
   any partial sequence "MGCEN", the actual deconvolved sequence might be
   "MLDGN".
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"

/* globals for fasta */
#define MAXWINDOW 64

#define EOSEQ 0 
#define ESS 49
#define MAXHASH 32
#define NMAP MAXHASH+1
#define NMAP_X 23	/* re-code NMAP for 'X' */
#define NMAP_Z 24	/* re-code NMAP for '*' */

#ifndef MAXSAV
#define MAXSAV 10
#endif

static char *verstr="3.36 June 2000";

#ifdef TFASTF
extern int aatran(unsigned char *ntseq, unsigned char *aaseq, int maxs, int frame);
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
   int     start0;		/* score from global match */
   int     dp;			/* diagonal of match */
   int     start;		/* start of match in lib seq */
   int     stop;		/* end of match in lib seq */
};

struct swstr { int H, E;};

struct bdstr { int CC, DD, CP, DP;};

struct hlstr { int next, pos;};

#define NM_MAX 10

struct f_struct {
  struct dstruct *diag;
  struct savestr vmax[MAXSAV];	/* best matches saved for one sequence */
  struct savestr *vptr[MAXSAV];
  struct savestr *lowmax;
  int nsave;			/* number of results saved in f_str->vmax */
  int ndo;
  int noff;			/* offset used in diagonal calculation */
  int nm0, nmoff;		/* number of segments, segment length */
  unsigned char *aa0;
  unsigned char *aa0t;
  int aa0ix;			/* counter for current offset in aa0t */
  int hmask;			/* hash constants */
  int *pamh1;			/* pam based array */
  int *pamh2;			/* pam based kfact array */
  struct hlstr *link, *harr;		/* hash arrays */
  int kshft;			/* shift width */
  int nsav, lowscor;		/* number of saved runs, worst saved run */
#ifdef TFASTF
  unsigned char *aa1x;
  int n10;
#endif
  struct bdstr *bss;
  struct swstr *ss;
  struct swstr *r_ss;
  int *waa;
  int *res;
  int max_res;
};

static int dmatch (unsigned char *aa0, int n0,
		   unsigned char *aa1, int n1,
		   int hoff, int window, 
		   int **pam2, int gdelval, int ggapval,
		   struct f_struct *f_str);

#ifdef ALLOCN0
void savemax (struct dstruct *, int, struct f_struct *);
#else
void savemax (struct dstruct *, struct f_struct *);
#endif

int m0_spam(unsigned char *, unsigned char *, int, struct savestr *, int **, struct f_struct *);
int m1_spam(unsigned char *, unsigned char *, int, struct savestr *, int **, int, struct f_struct *);

int sconn(struct savestr **, int nsave, int cgap, int pgap, struct f_struct *);
void kpsort(struct savestr **, int);
void kssort(struct savestr **, int);
int sconn_a(unsigned char *, int, struct f_struct *, int cgap, int pgap, struct a_struct *);
void kpsort(struct savestr **, int);

/* initialize for fasta */

void
init_work (unsigned char *aa0, int n0, 
	   struct pstruct *ppst,
	   struct f_struct **f_arg)
{
   int mhv, phv;
   int hmax;
   int i0, ii0, hv;
   int pamfact;
   int btemp;
   struct f_struct *f_str;
   /* these used to be globals, but do not need to be */
   int ktup, fact, kt1;

   int maxn0;
   int *pwaa;
   int i, j, q;
   struct bdstr *bss;
   struct swstr *ss, *r_ss;
   int *waa;
   int *res;

   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   ppst->sw_flag = 0;
   ppst->param_u.fa.pgap = ppst->gdelval + ppst->ggapval;
   pamfact = ppst->param_u.fa.pamfact;
   ktup = ppst->param_u.fa.ktup = 1;
   fact = ppst->param_u.fa.scfact * ktup;

   if (pamfact == -1) pamfact = 0;
   else if (pamfact == -2) pamfact = 1;

   /* fastf3 cannot work with lowercase symbols as low complexity;
      thus, NMAP must be disabled; this depends on aascii['X']  */
   if (ppst->hsq[NMAP_X] == NMAP ) {ppst->hsq[NMAP_X]=1;}
   if (ppst->hsq[NMAP_Z] == NMAP ) {ppst->hsq[NMAP_Z]=1;}

   /*   this does not work for share ppst structs, as in threads */
   /*else {fprintf(stderr," cannot find 'X'==NMAP\n");} */

   for (i0 = 1, mhv = -1; i0 <= ppst->nsq; i0++)
      if (ppst->hsq[i0] < NMAP && ppst->hsq[i0] > mhv) mhv = ppst->hsq[i0];

   if (mhv <= 0) {
      fprintf (stderr, " maximum hsq <=0 %d\n", mhv);
      exit (1);
   }

   for (f_str->kshft = 0; mhv > 0; mhv /= 2)
      f_str->kshft++;

/*      kshft = 2;	*/
   kt1 = ktup - 1;
   hv = 1;
   for (i0 = 0; i0 < ktup; i0++)
      hv = hv << f_str->kshft;
   hmax = hv;
   f_str->hmask = (hmax >> f_str->kshft) - 1;

   if ((f_str->aa0 = (unsigned char *) calloc(n0+1, sizeof(char))) == NULL) {
     fprintf (stderr, " cannot allocate f_str->aa0 array; %d\n",n0+1);
     exit (1);
   }
   for (i=0; i<n0; i++) f_str->aa0[i] = aa0[i];
   aa0 = f_str->aa0;

   if ((f_str->aa0t = (unsigned char *) calloc(n0+1, sizeof(char))) == NULL) {
     fprintf (stderr, " cannot allocate f_str0->aa0t array; %d\n",n0+1);
     exit (1);
   }
   f_str->aa0ix = 0;

   if ((f_str->harr = (struct hlstr *) calloc (hmax, sizeof (struct hlstr))) == NULL) {
     fprintf (stderr, " cannot allocate hash array; hmax: %d hmask: %d\n",
	      hmax,f_str->hmask);
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
   if ((f_str->link = (struct hlstr *) calloc (n0, sizeof (struct hlstr))) == NULL) {
     fprintf (stderr, " cannot allocate hash link array");
     exit (1);
   }

   for (i0 = 0; i0 < hmax; i0++) {
      f_str->harr[i0].next = -1;
      f_str->harr[i0].pos = -1;
   }

   for (i0 = 0; i0 < n0; i0++) {
      f_str->link[i0].next = -1;
      f_str->link[i0].pos = -1;
   }

   /* encode the aa0 array */
   /*
     this code has been modified to allow for mixed peptide sequences
      aa0[] = 5 8 9 3 4 NULL 5 12 3 7 2 NULL
      the 'NULL' character resets the hash position counter, to indicate that
      any of several residues can be in the same position.
      We also need to keep track of the number of times this has happened, so that
      we can redivide the sequence later

      i0 counts through the sequence
      ii0 counts through the hashed sequence

      */

   f_str->nm0 = 1;
   f_str->nmoff = -1;
   phv = hv = 0;
   for (i0 = ii0 = 0; i0 < kt1; i0++, ii0++) {
     if (aa0[i0] != ESS && aa0[i0] != 0) {
       hv = (hv << f_str->kshft) + ppst->hsq[aa0[i0]];
       phv += ppst->pam2[0][aa0[i0]][aa0[i0]] * ktup;
     }
     else {
       fprintf(stderr," sequence too short: %d\n",i0);
       exit(1);
     }
   }

   for (; i0 < n0; i0++, ii0++)
   {
     /* reset the counter and start hashing again */
     if (aa0[i0] == ESS || aa0[i0] == 0) {
       aa0[i0] = 0;	/* set ESS to 0 */
       /*       fprintf(stderr," converted ',' to 0\n");*/
       i0++;	/* skip over the blank */
       f_str->nm0++;
       if (f_str->nmoff < 0) f_str->nmoff = i0;
       phv = hv = 0;
       for (ii0=0; ii0 < kt1; i0++, ii0++) {
	 hv = (hv << f_str->kshft) + ppst->hsq[aa0[i0]];
	 phv += ppst->pam2[0][aa0[i0]][aa0[i0]] * ktup;
       }
     }
     hv = ((hv & f_str->hmask) << f_str->kshft) + ppst->hsq[aa0[i0]];
     f_str->link[i0].next = f_str->harr[hv].next;
     f_str->link[i0].pos = f_str->harr[hv].pos;
     f_str->harr[hv].next = i0;
     f_str->harr[hv].pos = ii0;
     if (pamfact) {
       f_str->pamh2[hv] = (phv += ppst->pam2[0][aa0[i0]][aa0[i0]] * ktup);
       phv -= ppst->pam2[0][aa0[i0 - kt1]][aa0[i0 - kt1]] * ktup;
     }
     else
       f_str->pamh2[hv] = fact * ktup;
   }
   if (f_str-> nmoff < 0) f_str->nmoff = n0;

/*
   fprintf(stderr," nmoff: %d/%d nm0: %d\n", f_str->nmoff, n0,f_str->nm0);
*/

   /*
   fprintf(stderr," hmax: %d\n",hmax);
   for ( hv=0; hv<hmax; hv++)
       fprintf(stderr,"%2d %c %3d %3d\n",hv,
	       (hv > 0 && hv < ppst->nsq ) ? ppst->sq[ppst->hsq[hv]] : ' ',
	       f_str->harr[hv].pos,f_str->harr[hv].next);
   fprintf(stderr,"----\n");
   for ( hv=0; hv<n0; hv++)
       fprintf(stderr,"%2d: %3d %3d\n",hv,
	       f_str->link[hv].pos,f_str->link[hv].next);

   */

/* this has been modified from 0..<ppst->nsq to 1..<=ppst->nsq because the
   pam2[0][0] is now undefined for consistency with blast
*/

   if (pamfact)
      for (i0 = 1; i0 <= ppst->nsq; i0++)
	 f_str->pamh1[i0] = ppst->pam2[0][i0][i0] * ktup;
   else
      for (i0 = 1; i0 <= ppst->nsq; i0++)
	 f_str->pamh1[i0] = fact;


   ppst->param_u.fa.cgap = shscore(aa0,f_str->nmoff-1,ppst->pam2[0],ppst->nsq)/3;
   if (ppst->param_u.fa.cgap > ppst->param_u.fa.bestmax/4)
     ppst->param_u.fa.cgap = ppst->param_u.fa.bestmax/4;

   f_str->ndo = 0;
#ifndef ALLOCN0
   if (f_str->diag==NULL) 
     f_str->diag = (struct dstruct *) calloc ((size_t)MAXDIAG,
					      sizeof (struct dstruct));
#else
   if (f_str->diag==NULL) 
     f_str->diag = (struct dstruct *) calloc ((size_t)n0,
					      sizeof (struct dstruct));
#endif

   if (f_str->diag == NULL)
   {
      printf (" cannot allocate diagonal arrays: %ld\n",
	      (long) MAXDIAG * (long) (sizeof (struct dstruct)));
      exit (1);
   }

#ifdef TFASTF
   if ((f_str->aa1x =(unsigned char *)calloc((size_t)ppst->maxlen+2,
					     sizeof(unsigned char)))
       == NULL) {
     fprintf (stderr, "cannot allocate aa1x array %d\n", ppst->maxlen+2);
     exit (1);
   }
   f_str->aa1x++;
#endif

   f_str->bss = (struct bdstr *) calloc((size_t)ppst->param_u.fa.optwid+4,
					sizeof(struct bdstr));
   f_str->bss++;
  
   /* allocate space for the scoring arrays */
   maxn0 = n0 + 4;
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

   if ((waa= (int *)calloc ((size_t)(ppst->nsq+1)*n0,sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate waa struct %3d\n",ppst->nsq*n0);
     exit(1);
   }

   pwaa = waa;
   for (i=0; i<=ppst->nsq; i++) {
     for (j=0;j<n0; j++) {
       *pwaa = ppst->pam2[0][i][aa0[j]];
       pwaa++;
     }
   }
   f_str->waa = waa;

   maxn0 = max(3*n0/2,MIN_RES);
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
#ifndef TFASTF
  char *pg_str="FASTF";
#else
  char *pg_str="TFASTF";
#endif

  sprintf (pstring1, "%s (%s) function [%s matrix (%d:%d)] ktup: %d\n join: %d, gap-pen: %d/%d, width: %3d",pg_str,verstr,
	       pstr->pamfile, pstr->pam_h,pstr->pam_l, pstr->param_u.fa.ktup, pstr->param_u.fa.cgap,
	       pstr->gdelval, pstr->ggapval, pstr->param_u.fa.optwid);

  if (pstr->param_u.fa.iniflag) strcat(pstring1," init1");
  /*
  if (pstr->zsflag==0) strcat(pstring1," not-scaled");
  else if (pstr->zsflag==1) strcat(pstring1," reg.-scaled");
  */

  if (pstring2 != NULL) {
     sprintf (pstring2, "; pg_name: %s\n; pg_ver: %s\n; pg_matrix: %s (%d:%d)\n\
; pg_gap-pen: %d %d\n; pg_ktup: %d\n; pg_cgap: %d\n",
	      pg_str,verstr,pstr->pamfile, pstr->pam_h,pstr->pam_l, pstr->gdelval,
              pstr->ggapval,pstr->param_u.fa.ktup,pstr->param_u.fa.cgap);
   }
}

void
close_work (unsigned char *aa0, int n0,
	    struct pstruct *ppst,
	    struct f_struct **f_arg)
{
  struct f_struct *f_str;

  f_str = *f_arg;

  if (f_str != NULL) {
    f_str->ss--;
    f_str->r_ss--;
    f_str->bss--;

    free(f_str->res);
    free(f_str->waa);
    free(f_str->r_ss);
    free(f_str->ss);
    free(f_str->bss);
    free(f_str->diag);
    free(f_str->link);
    free(f_str->pamh2); 
    free(f_str->pamh1);
    free(f_str->harr);
    free(f_str->aa0t);
    free(f_str->aa0);
    free(f_str);
    *f_arg = NULL;
  }
}

int do_fasta (unsigned char *aa0, int n0,
	      unsigned char *aa1, int n1,
	      struct pstruct *ppst, struct f_struct *f_str,
	      struct rstruct *rst, int *hoff)
{
   int     nd;		/* diagonal array size */
   int     lhval;
   int     kfact;
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
   int     tpos, npos;
   struct savestr *vmptr;
   int     scor, tmp;
   int     im, ib, nsave;
   int     cmps ();		/* comparison routine for ksort */
   int ktup, kt1;
   unsigned char   *aa1ptr;
   int  itt, lcont, ocont, loff;	/* lcont is returned by getlib to
					 * indicate there is more sequence
					 * remaining.  ocont is the previous
					 * value of lcont, for going back
					 * later. loff corrects maxn for the
					 * modified size of aa1 for continued
					 * sequences */

   ktup = ppst->param_u.fa.ktup;
   kt1 = ktup-1;

   if (n1 < ktup) {
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     return 1;
   }

   if (n0+n1+1 >= MAXDIAG) {
     fprintf(stderr,"n0,n1 too large: %d, %d\n",n0,n1);
     rst->score[0] = rst->score[1] = rst->score[2] = -1;
     return -1;
   }

   f_str->noff = n0 - 1;

#ifdef ALLOCN0
   nd = n0;
#endif

/*
	these initializations have been added to deal with reading
	sequences in chunks
*/

   aa1ptr = aa1;

   lcont = 0;
   ocont = 0;
   loff = 0;

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

   lpos = 0;
#ifndef ALLOCN0
   diagp = &f_str->diag[f_str->noff + kt1];
   for (; lpos < n1; lpos++, diagp++) {
     lhval = ppst->hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval].pos, npos = f_str->harr[lhval].next;
	  tpos >= 0; tpos = f_str->link[npos].pos, npos = f_str->link[npos].next) {
       if ((tscor = (dptr = &diagp[-tpos])->stop) >= 0) {
#else
   lposn0 = f_str->noff + lpos;
   for (; lpos < n1; lpos++, lposn0++) {
     lhval = ppst->hsq[aa1[lpos]];
     for (tpos = f_str->harr[lhval].pos; tpos >= 0; tpos = f_str->link[tpos].pos) {
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

   /* set up pointers for sorting */

   for (nsave = 0, vmptr = f_str->vmax; vmptr < &f_str->vmax[MAXSAV]; vmptr++) {
     if (vmptr->score > 0) {
       vmptr->score = m0_spam (aa0, aa1, n1, vmptr, ppst->pam2[0], f_str);
       f_str->vptr[nsave++] = vmptr;
     }
   }

   /* sort them */
   kssort (f_str->vptr, nsave);

   /*
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	     f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
	     f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
	     f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
	     f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
     for (im=f_str->vptr[ib]->start; im<=f_str->vptr[ib]->stop; im++)
       fprintf(stderr," %c:%c",ppst->sq[aa0[f_str->noff+im-f_str->vptr[ib]->dp]],
	       ppst->sq[aa1[im]]);
     fputc('\n',stderr);
   }

   fprintf(stderr,"---\n");
    */

   /* now use m_spam to re-evaluate */

   /*
   for (tpos = 0; tpos < n0; tpos++) {
     fprintf(stderr,"%c %2d",ppst->sq[aa0[tpos]],aa0[tpos]);
     if (tpos %10 == 9) fputc(' ',stderr);
   }
   fputc('\n',stderr);
   */

   f_str->aa0ix = 0;
   for (ib=0; ib < nsave; ib++) {
     if ((vmptr=f_str->vptr[ib])->score > 0) {
       vmptr->score = m1_spam (aa0, aa1, n1, vmptr,
			       ppst->pam2[0], ppst->pam_l, f_str);
     }
   }
   /* reset aa0 */
   for (tpos = 0; tpos < n0; tpos++) {
     if (aa0[tpos] >= 32) aa0[tpos] -= 32;
   }

   kssort(f_str->vptr,nsave);

   for ( ; nsave > 0; nsave--) 
     if (f_str->vptr[nsave-1]->score >0) break;

   if (nsave <= 0) {
     f_str->nsave = 0;
     rst->score[0] = rst->score[1] = rst->score[2] = 0;
     return 1;
   }
   else f_str->nsave = nsave;

   /*
   fprintf(stderr,"n0: %d; n1: %d; noff: %d\n",n0,n1,f_str->noff);
   for (ib=0; ib<nsave; ib++) {
     fprintf(stderr,"0: %4d-%4d  1: %4d-%4d  dp: %d score: %d\n",
	     f_str->noff+f_str->vptr[ib]->start-f_str->vptr[ib]->dp,
	     f_str->noff+f_str->vptr[ib]->stop-f_str->vptr[ib]->dp,
	     f_str->vptr[ib]->start,f_str->vptr[ib]->stop,
	     f_str->vptr[ib]->dp,f_str->vptr[ib]->score);
     for (im=f_str->vptr[ib]->start; im<=f_str->vptr[ib]->stop; im++)
       fprintf(stderr," %c:%c",ppst->sq[aa0[f_str->noff+im-f_str->vptr[ib]->dp]],
	       ppst->sq[aa1[im]]);
     fputc('\n',stderr);
   }

   fprintf(stderr,"---\n");
   */

   scor = sconn (f_str->vptr, nsave, ppst->param_u.fa.cgap, 
		 ppst->param_u.fa.pgap, f_str);

   for (vmptr=f_str->vptr[0],ib=1; ib<nsave; ib++)
     if (f_str->vptr[ib]->score > vmptr->score) vmptr=f_str->vptr[ib];

   rst->score[1] = vmptr->score;
   rst->score[0] = rst->score[2] = max (scor, vmptr->score);
   return 1;
}

void do_work (unsigned char *aa0, int n0,
	      unsigned char *aa1, int n1,
	      int frame,
	      struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, struct rstruct *rst)
{
  int hoff, n10, i;

#ifdef TFASTF 
  n10=aatran(aa1,f_str->aa1x,n1,frame);
  if (ppst->debug_lib)
    for (i=0; i<n10; i++)
      if (f_str->aa1x[i]>ppst->nsq) {
	fprintf(stderr,
		"residue[%d/%d] %d range (%d)\n",i,n1,
		f_str->aa1x[i],ppst->nsq);
	f_str->aa1x[i]=0;
	n10=i-1;
      }

  do_fasta (f_str->aa0, n0, f_str->aa1x, n10, ppst, f_str, rst, &hoff);
#else	/* FASTA */
  do_fasta (f_str->aa0, n0, aa1, n1, ppst, f_str, rst, &hoff);
#endif

}

void do_opt (unsigned char *aa0, int n0,
	     unsigned char *aa1, int n1,
	     int frame,
	     struct pstruct *pst,
	     struct f_struct *f_str,
	     struct rstruct *rst)
{
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

/* this version of spam() is designed to work with a collection of
   subfragments, selecting the best amino acid at each position so
   that, from each subfragment, each position is only used once.

   As a result, m_spam needs to know the number of fragments.

   In addition, it now requires a global alignment to the fragment
   and resets the start and stop positions

   */

int m1_spam (unsigned char *aa0, unsigned char *aa1, int n1,
	  struct savestr *dmax, int **pam2, int pam_l, struct f_struct *f_str)
{
  int     tpos, lpos, im, ii, nm, ci;
  int     tot, mtot, ctot, pv;
  int a0stop;
   struct {
     int     start, stop, score;
   } curv, maxv;
   register unsigned char *aa0p, *aa1p;
   unsigned char *aa0pt;

   lpos = dmax->start;			/* position in library sequence */
   tpos = lpos - dmax->dp + f_str->noff; /* position in query sequence */
   /* force global alignment, reset start*/
   if (tpos < lpos) {
     lpos = dmax->start -= tpos;
     tpos = 0;
   }
   else {
     tpos -= lpos;
     lpos = dmax->start = 0;
   }

   dmax->stop = dmax->start + (f_str->nmoff -2 - tpos);
   if (dmax->stop > n1) dmax->stop = n1;

   /*
   if (dmax->start < 0) {
     tpos = -dmax->start;
     lpos = dmax->start=0;
   }
   else tpos = 0;
   */

   aa1p = &aa1[lpos];
   aa0p = &aa0[tpos];

   nm = f_str->nm0;

   tot = curv.score = maxv.score = 0;
   for (; lpos <= dmax->stop; lpos++,aa0p++,aa1p++) {
     ctot = pam_l;
     ci = -1;
     for (im = 0, ii=0; im < nm; im++,ii+=f_str->nmoff) {
       if (aa0p[ii] < 32 && (pv = pam2[aa0p[ii]][*aa1p]) > ctot) {
	 ctot = pv;
	 ci = ii;
       }
     }
     tot += ctot;
     if (ci >= 0 && aa0p[ci] < 32) {
       aa0p[ci] +=  32;
     }
     
   }
   return tot;
}

int ma_spam (unsigned char *aa0, int n0, unsigned char *aa1,
	  struct savestr *dmax, struct pstruct *ppst,
	  struct f_struct *f_str)
{
  int **pam2;
  int     tpos, lpos, im, ii, nm, ci, lp0;
  int     tot, mtot, ctot, pv;
  struct {
    int     start, stop, score;
  } curv, maxv;
   register unsigned char *aa0p, *aa1p;
   unsigned char *aa0pt;
   int aa0t_flg;

   pam2 = ppst->pam2[0];
   aa0t_flg = 0;

   lpos = dmax->start;			/* position in library sequence */
   tpos = lpos - dmax->dp + f_str->noff; /* position in query sequence */
   lp0 = lpos = dmax->start;
   aa1p = &aa1[lpos];
   aa0p = &aa0[tpos];			/* real aa0 sequence */

			/* the destination aa0 sequence (without nulls) */
   aa0pt = &f_str->aa0t[f_str->aa0ix];

   curv.start = lpos;
   nm = f_str->nm0;

   /* sometimes, tpos may be > 0, with lpos = 0 - fill with 'X' */
   if (lpos == 0 && tpos > 0)
     for (ii = 0; ii < tpos; ii++) *aa0pt++ = 31;  /* filler character */

   tot = curv.score = maxv.score = 0;
   for (; lpos <= dmax->stop; lpos++) {
     ctot = ppst->pam_l;
     ci = -1;
     for (im = 0, ii=0; im < nm; im++,ii+=f_str->nmoff) {
       if (aa0p[ii] < 32 && (pv = pam2[aa0p[ii]][*aa1p]) > ctot) {
	 ctot = pv;
	 ci = ii;
       }
     }
     tot += ctot;
     if (ci >= 0) {
       if (ci >= n0) {fprintf(stderr," warning - ci off end %d/%d\n",ci,n0);}
       else {
	 *aa0pt++ = aa0p[ci];
	 aa0p[ci] +=  32;
	 aa0t_flg=1;
       }
     }
     aa0p++; aa1p++;
   }

   if (aa0t_flg) {
     dmax->dp -= f_str->aa0ix;		/* shift ->dp for aa0t */
     if ((ci=(int)(aa0pt-f_str->aa0t)) > n0) {
       fprintf(stderr," warning - aapt off %d/%d end\n",ci,n0);
     }
     else 
       *aa0pt++ = 0;			/* skip over NULL */

     aa0pt = &f_str->aa0t[f_str->aa0ix];
     aa1p = &aa1[lp0];

     /*
     for (im = 0; im < f_str->nmoff; im++)
       fprintf(stderr,"%c:%c,",ppst->sq[aa0pt[im]],ppst->sq[aa1p[im]]);
     fprintf(stderr,"- %3d (%3d:%3d)\n",dmax->score,f_str->aa0ix,lp0);
     */

     f_str->aa0ix += f_str->nmoff;	/* update offset into aa0t */
   }
   /*
      fprintf(stderr," ma_spam returning: %d\n",tot);
   */
   return tot;
}

int m0_spam (unsigned char *aa0, unsigned char *aa1, int n1,
	  struct savestr *dmax, int **pam2,
	  struct f_struct *f_str)
{
   int tpos, lpos, lend, im, ii, nm, ci;
   int     tot, mtot, ctot, pv;
   struct {
     int     start, stop, score;
   } curv, maxv;
   register unsigned char *aa0p, *aa1p;

   lpos = dmax->start;			/* position in library sequence */
   tpos = lpos - dmax->dp + f_str->noff; /* position in query sequence */
   if (tpos > 0) {
     if (lpos-tpos >= 0) {
       lpos = dmax->start -= tpos;    /* force global alignment, reset start*/
       tpos = 0;
     }
     else {
       tpos -= lpos;
       lpos = dmax->start = 0;
     }
   }

   nm = f_str->nm0;
   lend = dmax->stop;
   if (n1 - (lpos + f_str->nmoff-2) < 0 ) {
     lend = dmax->stop = (lpos - tpos) + f_str->nmoff-2;
     if (lend >= n1) lend = n1-1;
   }

   aa1p = &aa1[lpos];
   aa0p = &aa0[tpos];

   curv.start = lpos;

   tot = curv.score = maxv.score = 0;
   for (; lpos <= lend; lpos++) {
     ctot = -10000;
     for (im = 0, ii=0; im < nm; im++,ii+=f_str->nmoff) {
       if ((pv = pam2[aa0p[ii]][*aa1p]) > ctot) {
	 ctot = pv;
       }
     }
     tot += ctot;
     aa0p++; aa1p++;
   }

   /* reset dmax if necessary */

   return tot;
}

/* sconn links up non-overlapping alignments and calculates the score */

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
   int     lstart, tstart, tstop, ltmp, plstop, ptstart, ptstop;

/* sarr[] saves each alignment score/position, and provides a link
   back to the previous alignment that maximizes the score */

   pgap = 0;

/*	sort the score left to right in lib pos */
   kpsort (v, n);

   start = NULL;

/* for the remaining runs, see if they fit */
   for (i = 0, si = 0; i < n; i++) {

/* if the score is less than the gap penalty, it never helps */
      if (v[i]->score < cgap) continue;
      lstart = v[i]->start;

/* put the run in the group */
      sarr[si].vp = v[i];
      sarr[si].score = v[i]->score;
      sarr[si].next = NULL;

/* if it fits, then increase the score */
      for (sl = start; sl != NULL; sl = sl->next) {
	plstop = sl->vp->stop;
	/* if end < start or start > end, add score */
	if (plstop < lstart ) {
	  sarr[si].score = sl->score + v[i]->score + pgap;
	  break;
	}
      }

/*	now recalculate where the score fits - resort the scores */
      if (start == NULL)
	start = &sarr[si];
      else
	for (sj = start, so = NULL; sj != NULL; sj = sj->next) {
	  if (sarr[si].score > sj->score) { /* if new score > best score */
	       sarr[si].next = sj;	/* previous best linked to best */
	       if (so != NULL)		
		 so->next = &sarr[si];	/* old best points to new best */
	       else
		  start = &sarr[si];
	       break;
	    }
	    so = sj;			/* old-best saved in so */
	 }
      si++;
   }

   if (start != NULL)
      return (start->score);
   else
      return (0);
}

void
kssort (struct savestr **v, int n)
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

/* sorts alignments from right to left (back to front) based on stop */

void
krsort (v, n)
struct savestr *v[];
int     n;
{
   int     gap, i, j;
   struct savestr *tmp;

   for (gap = n / 2; gap > 0; gap /= 2)
      for (i = gap; i < n; i++)
	 for (j = i - gap; j >= 0; j -= gap)
	 {
	    if (v[j]->stop > v[j + gap]->stop)
	       break;
	    tmp = v[j];
	    v[j] = v[j + gap];
	    v[j + gap] = tmp;
	 }
}

int  do_walign (unsigned char *aa0, int n0,
		unsigned char *aa1, int n1,
		int frame,
		struct pstruct *ppst, 
		struct f_struct *f_str, 
		int **ares, int *nres, struct a_struct *aln)
{
  int hoff, n10;
  struct rstruct rst;
  int ib;
  unsigned char *aa0t, *aa1p;
  struct savestr *vmptr;

#ifdef TFASTF
  f_str->n10 = n10 = aatran(aa1,f_str->aa1x,n1,frame);
  aa1p = f_str->aa1x;
  aln->qlrev = 0;
  aln->qlfact = 1;
  aln->llfact = aln->llmult = 3;
  aln->frame = 0;
  if (frame > 3) aln->llrev = 1;
#else
  aln->llfact = aln->qlfact = aln->llmult = 1;
  aln->llrev = aln->qlrev = 0;
  aln->frame = 0;
  n10 = n1;
  aa1p = aa1;
#endif

  do_fasta(f_str->aa0, n0, aa1p, n10, ppst, f_str, &rst, &hoff);

  /* the alignment portion takes advantage of the information left
     over in f_str after do_fasta is done.  in particular, it is
     easy to run a modified sconn() to produce the alignments.

     unfortunately, the alignment display routine wants to have
     things encoded as with bd_align and sw_align, so we need to do that.
     */

  if ((aa0t = (unsigned char *)calloc(n0+1,sizeof(unsigned char)))==NULL) {
    fprintf(stderr," cannot allocate aa0t %d\n",n0+1);
    exit(1);
  }

   kssort (f_str->vptr, f_str->nsave);
   f_str->aa0ix = 0;
   if (f_str->nsave > f_str->nm0) f_str->nsave = f_str->nm0;
   for (ib=0; ib < f_str->nm0; ib++) {
     if (f_str->vptr[ib]->score > 0) {
       f_str->vptr[ib]->score = 
	 ma_spam (f_str->aa0, n0, aa1p, f_str->vptr[ib], ppst, f_str);
     }
   }

   /* after ma_spam is over, we need to reset aa0 */
   for (ib = 0; ib < n0; ib++) {
     if (f_str->aa0[ib] >= 32) f_str->aa0[ib] -= 32;
   }

   kssort(f_str->vptr,f_str->nsave);

   for ( ; f_str->nsave > 0; f_str->nsave--) 
     if (f_str->vptr[f_str->nsave-1]->score >0) break;

  *nres = sconn_a (aa0t,n0,f_str,ppst->param_u.fa.cgap,
		   ppst->param_u.fa.pgap, aln);
  free(aa0t);

  *ares = f_str->res;
  return rst.score[0];
}

/* this version of sconn is modified to provide alignment information */

int sconn_a (unsigned char *aa0, int n0, 
	     struct f_struct *f_str,int cgap, int pgap,
	     struct a_struct *aln)
{
   int     i, si, cmpp (), n;
   unsigned char *aa0p;
   int sx, dx, doff;

   struct savestr **v;
   struct slink {
     int     score;
     struct savestr *vp;
     struct slink *snext;
     struct slink *aprev;
   } *start, *sl, *sj, *so, sarr[MAXSAV];
   int     lstop, plstart;
   int *res, nres, tres;

/*	sort the score left to right in lib pos */

   v = f_str->vptr;
   n = f_str->nsave;

   pgap = 0;

   krsort (v, n);	/* sort from left to right in library */

   start = NULL;

/*	for each alignment, see if it fits */

   for (i = 0, si = 0; i < n; i++) {
/*	if the score is less than the join threshold, skip it */

     if (v[i]->score < cgap) continue;

     lstop = v[i]->stop;		/* have right-most lstart */

/*	put the alignment in the group */

     sarr[si].vp = v[i];
     sarr[si].score = v[i]->score;
     sarr[si].snext = NULL;
     sarr[si].aprev = NULL;

/* 	if it fits, then increase the score */
/* start points to a sorted (by total score) list of candidate
   overlaps */

     for (sl = start; sl != NULL; sl = sl->snext) { 
       plstart = sl->vp->start;
       if (plstart > lstop ) {
	 sarr[si].score = sl->score + v[i]->score + pgap;
	 sarr[si].aprev = sl;
	 break;		/* quit as soon as the alignment has been added */
       }
     }

/* now recalculate the list of best scores */
     if (start == NULL)
       start = &sarr[si];	/* put the first one in the list */
     else
       for (sj = start, so = NULL; sj != NULL; sj = sj->snext) {
	 if (sarr[si].score > sj->score) { /* new score better than old */
	   sarr[si].snext = sj;		/* snext best after new score */
	   if (so != NULL)
	     so->snext = &sarr[si];	/* prev_best->snext points to best */
	   else  start = &sarr[si];	/* start points to best */
	   break;			/* stop looking */
	 }
	 so = sj;		/* previous candidate best */
       }
     si++;				/* increment to snext alignment */
   }

   /* we have the best set of alignments, write them to *res */
   if (start != NULL) {
     res = f_str->res;	/* set a destination for the alignment ops */
     tres = nres = 0;	/* alignment op length = 0 */
     aa0p = aa0;	/* point into query (needed for calcons later) */
     aln->min1 = start->vp->start+1;	/* start in library */
     aln->min0 = 1;			/* start in query */
     for (sj = start; sj != NULL; sj = sj->aprev ) {
       doff = (int)(aa0p-aa0) - (sj->vp->start-sj->vp->dp+f_str->noff);
       /*
       fprintf(stderr,"doff: %3d\n",doff);
       */
       for (dx=sj->vp->start,sx=sj->vp->start-sj->vp->dp+f_str->noff;
	    dx <= sj->vp->stop; dx++) {
	 *aa0p++ = f_str->aa0t[sx++];	/* copy residue into aa0 */
	 tres++;			/* bump alignment counter */
	 res[nres++] = 0;		/* put 0-op in res */
       }
       sj->vp->dp -= doff;
       if (sj->aprev != NULL) {
	 if (sj->aprev->vp->start - sj->vp->stop - 1 > 0 )
	 /* put an insert op into res to get to next aligned block */
	   tres += res[nres++] = (sj->aprev->vp->start - sj->vp->stop - 1);
       }
       /*
       fprintf(stderr,"t0: %3d, tx: %3d, l0: %3d, lx: %3d, dp: %3d noff: %3d, score: %3d\n",
	       sj->vp->start - sj->vp->dp + f_str->noff,
	       sj->vp->stop - sj->vp->dp + f_str->noff,
	       sj->vp->start,sj->vp->stop,sj->vp->dp,
	       f_str->noff,sj->vp->score);
       fprintf(stderr,"%3d - %3d: %3d\n",
	       sj->vp->start,sj->vp->stop,sj->vp->score);
       */
       aln->max1 = sj->vp->stop+1;
       aln->max0 = aln->max1 - sj->vp->dp + f_str->noff;
     }

     /*
     fprintf(stderr,"(%3d - %3d):(%3d - %3d)\n",
     aln->min0,aln->max0,aln->min1,aln->max1);
     */

     /* now replace f_str->aa0t with aa0 */
     for (i=0; i<n0; i++) f_str->aa0t[i] = aa0[i];

     return tres;
   }
   else return (0);
}


/* calculate the 100% identical score */
shscore(unsigned char *aa0, int n0, int **pam2, int nsq)
{
  int i, sum;
  for (i=0,sum=0; i<n0; i++)
    if (aa0[i]!=0 && aa0[i]<=nsq) sum += pam2[aa0[i]][aa0[i]];
  return sum;
}

int calcons(unsigned char *aa0, int n0,
	    unsigned char *aa1, int n1,
	    int *res, int nres, int *nc,
	    struct a_struct *aln, struct pstruct pst,
	    char *seqc0, char *seqc1,
	    struct f_struct *f_str)
{
  int i0, i1, nn1, n0t;
  int op, lenc, len_gap, nd, ns, itmp;
  unsigned char *aa1p;
  char *sp0, *sp1;
  int *rp;
  
#ifndef TFASTF
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* first fill in the ends */
  aln->min0--; aln->min1--;
  n0 -= (f_str->nm0-1);

  if (min(aln->min0,aln->min1)<aln->llen || aln->showall==1) {
    /* will we show all the start ?*/ 
    aln->smins=0;
    aln->mins = min(aln->min1,aln->llen/2);
    aancpy(seqc1,(char *)(aa1p+aln->min1-aln->mins),aln->mins,pst);
    aln->smin1 = aln->min1-aln->mins;
    if ((aln->mins-aln->min0)>0) {
      memset(seqc0,' ',aln->mins-aln->min0);
      aancpy(seqc0+aln->mins-aln->min0,(char *)f_str->aa0t,aln->min0,pst);
      aln->smin0 = 0;
    }
    else {
      aancpy(seqc0,(char *)f_str->aa0t+aln->min0-aln->mins,aln->mins,pst);
      aln->smin0 = aln->min0-aln->mins;
    }
  }
  else {
    aln->mins= min(aln->llen/2,min(aln->min0,aln->min1));
    aln->smins=aln->mins;
    aln->smin0=aln->min0;
    aln->smin1=aln->min1;
    aancpy(seqc0,(char *)f_str->aa0t+aln->min0-aln->mins,aln->mins,pst);
    aancpy(seqc1,(char *)aa1p+aln->min1-aln->mins,aln->mins,pst);
  }

/* now get the middle */

  sp0 = seqc0+aln->mins;
  sp1 = seqc1+aln->mins;
  rp = res;
  n0t = lenc = len_gap = aln->nident = aln->ngap_q = aln->ngap_l
    = aln->nfs = op = 0;
  i0 = aln->min0;
  i1 = aln->min1;
  
  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      *sp0 = pst.sq[f_str->aa0t[i0++]];
      *sp1 = pst.sq[aa1p[i1++]];
      n0t++;
      lenc++;
      if (toupper(*sp0) == toupper(*sp1)) aln->nident++;
      sp0++; sp1++;
    }
    else {
      if (op==0) { op = *rp++;}
      if (op>0) {
	*sp0++ = '-';
	*sp1++ = pst.sq[aa1p[i1++]];
	op--;
	len_gap++;
	lenc++;
      }
      else {
	*sp0++ = pst.sq[f_str->aa0t[i0++]];
	*sp1++ = '-';
	op++;
	n0t++;
	len_gap++;
	lenc++;
      }
    }
  }

  *nc = lenc-len_gap;
  /* now we have the middle, get the right end */
  /* ns is amount to be shown */
  /* nd is amount remaining to be shown */
  ns = aln->mins + lenc + aln->llen;
  ns -= (itmp = ns %aln->llen);
  if (itmp>aln->llen/2) ns += aln->llen;
  nd = ns - (aln->mins+lenc);
  if (nd > max(n0t-aln->max0,nn1-aln->max1)) nd = max(n0t-aln->max0,nn1-aln->max1);
  
  if (aln->showall==1) {
    nd = max(n0t-aln->max0,nn1-aln->max1);	/* reset for showall=1 */
    /* get right end */
    /* there isn't any aa0 to get */
    memset(seqc0+aln->mins+lenc,' ',n0t-aln->max0);
    aancpy(seqc1+aln->mins+lenc,(char *)aa1p+aln->max1,nn1-aln->max1,pst);
    /* fill with blanks - this is required to use one 'nc' */
    memset(seqc0+aln->mins+lenc+n0t-aln->max0,' ',nd-(n0t-aln->max0));
    memset(seqc1+aln->mins+lenc+nn1-aln->max1,' ',nd-(nn1-aln->max1));
  }
  else {
    memset(seqc0+aln->mins+lenc,' ',nd);
    if ((nd-(nn1-aln->max1))>0) {
      aancpy(seqc1+aln->mins+lenc,(char *)aa1p+aln->max1,nn1-aln->max1,pst);
      memset(seqc1+aln->mins+lenc+nn1-aln->max1,' ',nd-(nn1-aln->max1));
    }
    else aancpy(seqc1+aln->mins+lenc,(char *)aa1p+aln->max1,nd,pst);
  }
  
  return aln->mins+lenc+nd;
}

/* build an array of match/ins/del - length strings */
int calc_code(const unsigned char *aa0, const int n0,
	      const unsigned char *aa1, const int n1,
	      int *res, int nres,
	      struct a_struct *aln, struct pstruct pst,
	      char *al_str, int al_str_n, struct f_struct *f_str)
{
  int i0, i1, nn1;
  int op, lenc, nd, ns, itmp;
  int p_op, op_cnt;
  const unsigned char *aa1p;
  char tmp_cnt[20];
  char sp0, sp1, *sq;
  int *rp;

  if (pst.ext_sq_set) {
    sq = pst.sqx;
  }
  else {
    sq = pst.sq;
  }

#ifndef TFASTA
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* first fill in the ends */
  aln->min0--; aln->min1--;

  rp = res;
  lenc = aln->nident = aln->ngap_q = aln->ngap_l = aln->nfs = op = p_op = 0;
  op_cnt = 0;

  i0 = aln->min0;
  i1 = aln->min1;
  tmp_cnt[0]='\0';
  
  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {
      if (p_op == 0) { 	op_cnt++;}
      else {
	update_code(al_str,al_str_n-strlen(al_str),p_op,op_cnt);
	op_cnt = 1; p_op = 0;
      }

      op = *rp++;
      lenc++;


      sp0 = sq[aa0[i0++]];
      sp1 = sq[aa1p[i1++]];
      if (toupper(sp0) == toupper(sp1)) aln->nident++;
      else if (pst.dnaseq==1) {
	if ((toupper(sp0) == 'T' && toupper(sp1) == 'U') ||
	    (toupper(sp0)=='U' && toupper(sp1)=='T')) aln->nident++;
	else if (toupper(sp0) == 'N') aln->ngap_q++;
	else if (toupper(sp1) == 'N') aln->ngap_l++;
      }
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {
	if (p_op == 1) { op_cnt++;}
	else {
	  update_code(al_str,al_str_n - strlen(al_str),p_op,op_cnt);
	  op_cnt = 1; p_op = 1;
	}
	op--; lenc++; i1++; aln->ngap_q++;
      }
      else {
	if (p_op == 2) { op_cnt++;}
	else {
	  update_code(al_str,al_str_n - strlen(al_str),p_op,op_cnt);
	  op_cnt = 1; p_op = 2;
	}
	op++; lenc++; i0++; aln->ngap_l++;
      }
    }
  }
  update_code(al_str,al_str_n - strlen(al_str),p_op,op_cnt);

  return lenc;
}

update_code(char *al_str, int al_str_max, int op, int op_cnt) {

  char op_char[4]={"=-+"};
  char tmp_cnt[20];

  sprintf(tmp_cnt,"%c%d",op_char[op],op_cnt);
  strncat(al_str,tmp_cnt,al_str_max);
}

int calc_id(unsigned char *aa0, int n0,
	    unsigned char *aa1, int n1,
	    int *res, int nres,
	    struct a_struct *aln, struct pstruct pst,
	    struct f_struct *f_str)
{
  int i0, i1, nn1, n0t;
  int op, lenc, len_gap, nd, ns, itmp;
  unsigned char *aa1p;
  int sp0, sp1;
  int *rp;
  
#ifndef TFASTF
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* first fill in the ends */
  aln->min0--; aln->min1--;
  n0 -= (f_str->nm0-1);

  /* now get the middle */
  rp = res;
  n0t = lenc = len_gap = aln->nident = aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  i0 = aln->min0;
  i1 = aln->min1;
  
  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      sp0 = pst.sq[f_str->aa0t[i0++]];
      sp1 = pst.sq[aa1p[i1++]];
      n0t++;
      lenc++;
      if (toupper(sp0) == toupper(sp1)) aln->nident++;
    }
    else {
      if (op==0) { op = *rp++;}
      if (op>0) {
	i1++;
	op--;
	len_gap++;
	lenc++;
      }
      else {
	i0++;
	op++;
	n0t++;
	len_gap++;
	lenc++;
      }
    }
  }
  
  return lenc-len_gap;
}


#ifdef PCOMPLIB
#include "p_mw.h"
void
update_params(struct qmng_str *qm_msg, struct pstruct *ppst)
{
  ppst->n0 = qm_msg->n0;
}
#endif
