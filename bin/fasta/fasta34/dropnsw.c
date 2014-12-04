/* copyright (c) 1994, 1995, 1996 William R. Pearson */

/* $Name: fa34t20b3 $ - $Id: dropnsw.c,v 1.16 2002/08/28 21:09:17 wrp Exp $ */

/*
  this is a slower version of dropgsw.c that implements the Smith-Waterman
  algorithm.  It lacks the shortcuts in dropgsw.c that prevent scores less
  than the penalty for the first residue in a gap from being generated.

  Thus, dropnsw.c should be used for tests with very large gap penalties,
  and is more appropriate for programs like prss3, which are interested
  in accurate low scores.
*/

/* the do_walign() code in this file is not thread_safe */
/* init_work(), do_work(), are thread safe */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"

static char *verstr="3.39 May 2001";

struct swstr { int H, E;};

struct f_struct {
  struct swstr *ss;
  struct swstr *f_ss;
  struct swstr *r_ss;
  int *waa0;
  int *waa1;
  int *res;
  double aa0_f[MAXSQ];
  double *kar_p;
};

/* initialize for Smith-Waterman optimal score */

void init_work (unsigned char *aa0, int n0,
		struct pstruct *ppst,
		struct f_struct **f_arg)
{
   int maxn0;
   int *pwaa;
   int i, j, q;
   struct swstr *ss, *f_ss, *r_ss;
   struct f_struct *f_str;
   int *waa;
   int *res;
   int nsq, ip;

   if (ppst->ext_sq_set) {
     nsq = ppst->nsqx; ip = 1;
   }
   else {
     nsq = ppst->nsq; ip = 0;
   }

   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   /* allocate space for the scoring arrays */
   maxn0 = n0 + 2;
   if ((ss = (struct swstr *) calloc (maxn0, sizeof (struct swstr)))
	 == NULL) {
     fprintf (stderr, "cannot allocate ss array %3d\n", n0);
     exit (1);
   }
   ss++;
   f_str->ss = ss;

   if ((f_ss = (struct swstr *) calloc (maxn0, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, "cannot allocate f_ss array %3d\n", n0);
     exit (1);
   }
   f_ss++;
   f_str->f_ss = f_ss;

   if ((r_ss = (struct swstr *) calloc (n0+2, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, "cannot allocate r_ss array %3d\n", n0);
     exit (1);
   }
   r_ss++;
   f_str->r_ss = r_ss;

   if ((waa= (int *)malloc (sizeof(int)*(nsq+1)*n0)) == NULL) {
     fprintf(stderr,"cannot allocate waa0 struct %3d\n",nsq*n0);
     exit(1);
   }
   f_str->waa0 = waa;		/* waa0 => -S dependent pam2[] */

   pwaa = waa;
   for (i=0; i<=nsq; i++) {
     for (j=0;j<n0; j++) {
       *pwaa = ppst->pam2[ip][i][aa0[j]];
       pwaa++;
     }
   }

   if ((waa= (int *)malloc (sizeof(int)*(nsq+1)*n0)) == NULL) {
     fprintf(stderr,"cannot allocate waa1 struct %3d\n",nsq*n0);
     exit(1);
   }
   f_str->waa1 = waa;		/* waa1 has universal pam2[] */

   pwaa = waa;
   for (i=0; i<=nsq; i++) {
     for (j=0;j<n0; j++) {
       *pwaa = ppst->pam2[0][i][aa0[j]];
       pwaa++;
     }
   }

   maxn0 = max(3*n0/2,MIN_RES);
   if ((res = (int *)calloc((size_t)maxn0,sizeof(int)))==NULL) {
     fprintf(stderr,"cannot allocate alignment results array %d\n",maxn0);
     exit(1);
   }
   f_str->res = res;

   *f_arg = f_str;
}

void close_work (unsigned char *aa0, int n0,
		 struct pstruct *ppst, struct f_struct **f_arg)
{
  struct f_struct *f_str;

  f_str = *f_arg;

  if (f_str != NULL) {
    f_str->ss--;
    f_str->f_ss--;
    f_str->r_ss--;
    free(f_str->ss);
    free(f_str->f_ss);
    free(f_str->r_ss);
    free(f_str->waa0);
    free(f_str->waa1);
    free(f_str->res);
    free(f_str);
    *f_arg = NULL;
  }
}


/* pstring1 is a message to the manager, currently 512 */
/*void get_param(struct pstruct *pstr,char *pstring1)*/
void    get_param (struct pstruct *pstr, char *pstring1, char *pstring2)
{
  char *pg_str="Smith-Waterman";

#ifndef GAP_OPEN
   sprintf (pstring1, " Smith-Waterman (%s) function [%s matrix (%d:%d)%s], gap-penalty: %d/%d",
#else
   sprintf (pstring1, " Smith-Waterman (%s) function [%s matrix (%d:%d)%s], open/ext: %d/%d",
#endif
	    verstr,pstr->pamfile, pstr->pam_h,pstr->pam_l, 
	    (pstr->ext_sq_set) ? "xS":"\0",
	    pstr->gdelval, pstr->ggapval);
   if (pstring2 != NULL) {
#ifndef GAP_OPEN
     sprintf(pstring2,"; pg_name: %s\n; pg_ver: %s\n; pg_matrix: %s (%d:%d)%s\n; pg_gap-pen: %d %d\n",
#else
     sprintf(pstring2,"; pg_name: %s\n; pg_ver: %s\n; pg_matrix: %s (%d:%d)%s\n; pg_open-ext: %d %d\n",
#endif
	     pg_str,verstr,pstr->pamfile,pstr->pam_h,pstr->pam_l,
	     (pstr->ext_sq_set) ? "xS":"\0",
	     pstr->gdelval,pstr->ggapval);
   }
}


void do_work (unsigned char *aa0, int n0,
	      unsigned char *aa1, int n1,
	      int frame,
	      struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg, struct rstruct *rst)
{
   register unsigned char *aa0p, *aa1p;
   register struct swstr *ssj;
   struct swstr *ss, *f_ss, *r_ss;
   register int *pwaa;
   int *waa;
   register int i, j;
   int     e, f, h, p;
   int     q, r, m;
   int     score;

   double lambda, H, K;

   rst->escore = 1.0;
   rst->segnum = rst->seglen = 1;

   waa = f_str->waa0;
   ss = f_str->ss;
   f_ss = f_str->f_ss;
   r_ss = f_str->r_ss;

#ifndef GAP_OPEN
   q = -(ppst->gdelval - ppst->ggapval);
#else
   q = -ppst->gdelval;
#endif
   r = -ppst->ggapval;
   m = q + r;

   /* initialize 0th row */
   for (ssj=ss; ssj<&ss[n0]; ssj++) {
     ssj->H = 0;
     ssj->E = -q;
   }

   score = 0;
   aa1p = aa1;
   while (*aa1p) {
     h = p = 0;
     f = -q;
     pwaa = waa + (*aa1p++ * n0);
     for (ssj = ss, aa0p = aa0; ssj < ss+n0; ssj++) {
       if ((h =   h     - m) > (f =   f     - r)) f = h;
       if ((h = ssj->H - m) > (e = ssj->E - r)) e = h;
       h = p + *pwaa++;
       if (h < 0 ) h = 0;
       if (h < f ) h = f;
       if (h < e ) h = e;
       p = ssj->H;
       ssj->H = h;
       ssj->E = e;
       if (h > score) score = h;
     }
   }				/* done with forward pass */

   rst->score[0] = score;

   if(ppst->zsflag == 6 || ppst->zsflag == 16 && 
     (do_karlin(aa1, n1, ppst->pam2[0], ppst,f_str->aa0_f, 
		f_str->kar_p, &lambda, &H)>0)) {
     rst->comp = 1.0/lambda;
     rst->H = H;
   }
  else {rst->comp = rst->H = -1.0;}
}				/* here we should be all done */

void    do_opt (unsigned char *aa0, int n0,
		unsigned char *aa1, int n1,
		int frame,
		struct pstruct *pst, struct f_struct *f_str,
		struct rstruct *rstr)
{
}

int do_walign (unsigned char *aa0, int n0,
		 unsigned char *aa1, int n1,
		 int frame,
		 struct pstruct *ppst, 
		 struct f_struct *f_str, 
		 int **ares, int *nres, struct a_struct *aln)
{
   register unsigned char *aa0p, *aa1p;
   register int *pwaa;
   register int i, j;
   register struct swstr *ssj;
   struct swstr *f_ss, *r_ss, *ss;
   int *res, *waa;
   int e, f, h, p;
   int     q, r, m;
   int     score;
   int cost, I, J, K, L;

   ss = f_str->ss;
   f_ss = f_str->f_ss;
   r_ss = f_str->r_ss;
   res = f_str->res;
   waa = f_str->waa1;

#ifndef GAP_OPEN
   q = -(ppst->gdelval - ppst->ggapval);
#else
   q = -ppst->gdelval;
#endif

   r = -ppst->ggapval;
   m = q + r;

   /* initialize 0th row */
   for (ssj=ss; ssj<ss+n0; ssj++) {
     ssj->H = 0;
     ssj->E = -q;
   }

   score = 0;
   aa1p = aa1;
   i = 0;
   while (*aa1p) {
     h = p = 0;
     f = -q;
     pwaa = waa + (*aa1p++ * n0);
     for (ssj = ss, aa0p = aa0; ssj < ss+n0; ssj++) {
       if ((h =   h     - m) > (f =   f     - r)) f = h;
       if ((h = ssj->H - m) > (e = ssj->E - r)) e = h;
       h = p + *pwaa++;
       if (h < 0 ) h = 0;
       if (h < f ) h = f;
       if (h < e ) h = e;
       p = ssj->H;
       ssj->H = h;
       ssj->E = e;
       if (h > score) {
	 score = h;
	 I = i;
	 J = (int)(ssj-ss);
       }
     }
     i++;
   }				/* done with forward pass */
   if (score <= 0) return 0;

  /* to get the start point, go backwards */
  
  cost = K = L = 0;
  for (ssj=ss+J; ssj>=ss; ssj--) ssj->H= ssj->E= -1;
  
  for (i=I; i>=0; i--) {
    h = f = -1;
    p = (i == I) ? 0 : -1;
    for (ssj=ss+J,aa0p = &aa0[J]; ssj>=ss; ssj--,aa0p--) {
      f = max (f,h-q)-r;
      ssj->E=max(ssj->E,ssj->H-q)-r;
      h = max(max(ssj->E,f),p+ppst->pam2[0][*aa0p][aa1[i]]);
      p = ssj->H;
      ssj->H=h;
      if (h > cost) {
	cost = h;
	K = i;
	L = (int)(ssj-ss);
	if (cost >= score) goto found;
      }
    }
  }
  
found:	

  /*  printf(" %d: L: %3d-%3d/%3d; K: %3d-%3d/%3d\n",score,L,J,n0,K,I,n1); */

/* in the f_str version, the *res array is already allocated at 4*n0/3 */

  aln->max0 = J+1; aln->min0 = L+1; aln->max1 = I+1; aln->min1 = K+1;
  
  ALIGN(&aa1[K-1],&aa0[L-1],I-K+1,J-L+1,ppst->pam2[0],q,r,res,nres,f_str);

/*   DISPLAY(&aa0[L-1],&aa1[K-1],J-L+1,I-K+1,res,L,K,ppst->sq); */

   /* return *res and nres */

  aln->llfact = aln->llmult = aln->qlfact = 1;
  aln->qlrev = aln->llrev = 0;
  aln->frame = 0;

  *ares = f_str->res;
  return score;
}

static int CHECK_SCORE(unsigned char *A, unsigned char *B, int M, int N,
		       int *S, int **W, int G, int H, int *nres);

#define gap(k)  ((k) <= 0 ? 0 : g+h*(k))	/* k-symbol indel cost */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

						/* Append "Delete k" op */
#define DEL(k)				\
{ if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ if (last > 0)				\
    last = sapp[-1] += (k);	 	\
  else					\
    last = *sapp++ = (k);		\
}

#define REP { last = *sapp++ = 0; }		/* Append "Replace" op */

/* align(A,B,M,N,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

static int align(unsigned char *A, unsigned char *B, int M, int N,
		 int tb, int te, int **w, int g, int h, 
		 struct f_struct *f_str)
{
  int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;
  int c1, c2;

{ register int   i, j;
  register int c, e, d, s;
  int m, t, *wa;
  struct swstr *f_ss, *r_ss;

  m = g + h;

  f_ss = f_str->f_ss;
  r_ss = f_str->r_ss;

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0)
    { if (M > 0) INS(M)
      return -gap(M);
    }
  if (M <= 1)
    { if (M <= 0)
        { DEL(N);
          return -gap(N);
        }
      if (tb < te) tb = te;
      midc = (tb-h) - gap(N);
      midj = 0;
      wa = w[A[1]];
      for (j = 1; j <= N; j++)
        { c = -gap(j-1) + wa[B[j]] - gap(N-j);
          if (c > midc)
            { midc = c;
              midj = j;
            }
        }
      if (midj == 0)
        { DEL(N) INS(1) }
      else
        { if (midj > 1) DEL(midj-1)
          REP
          if (midj < N) DEL(N-midj)
        }
      return midc;
    }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  f_ss[0].H = 0;		/*   Compute H(M/2,k) & E(M/2,k) for all k */
  t = -g;
  for (j = 1; j <= N; j++)
    { f_ss[j].H = t = t-h;
      f_ss[j].E = t-g;
    }
  t = tb;
  for (i = 1; i <= midi; i++)
    { s = f_ss[0].H;
      f_ss[0].H = c = t = t-h;
      e = t-g;
      wa = w[A[i]];
      for (j = 1; j <= N; j++)
        { if ((c =   c   - m) > (e =   e   - h)) e = c;
          if ((c = f_ss[j].H - m) > (d = f_ss[j].E - h)) d = c;
          c = s + wa[B[j]];
          if (e > c) c = e;
          if (d > c) c = d;
          s = f_ss[j].H;
          f_ss[j].H = c;
          f_ss[j].E = d;
        }
    }
  f_ss[0].E = f_ss[0].H;

  r_ss[N].H = 0;			/* Reverse phase:                          */
  t = -g;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
  for (j = N-1; j >= 0; j--)
    { r_ss[j].H = t = t-h;
      r_ss[j].E = t-g;
    }
  t = te;
  for (i = M-1; i >= midi; i--)
    { s = r_ss[N].H;
      r_ss[N].H = c = t = t-h;
      e = t-g;
      wa = w[A[i+1]];
      for (j = N-1; j >= 0; j--)
        { if ((c =   c   - m) > (e =   e   - h)) e = c;
          if ((c = r_ss[j].H - m) > (d = r_ss[j].E - h)) d = c;
          c = s + wa[B[j+1]];
          if (e > c) c = e;
          if (d > c) c = d;
          s = r_ss[j].H;
          r_ss[j].H = c;
          r_ss[j].E = d;
        }
    }
  r_ss[N].E = r_ss[N].H;

  midc = f_ss[0].H+r_ss[0].H;		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  for (j = 0; j <= N; j++)
    if ((c = f_ss[j].H + r_ss[j].H) >= midc)
      if (c > midc || f_ss[j].H != f_ss[j].E && r_ss[j].H == r_ss[j].E)
        { midc = c;
          midj = j;
        }
  for (j = N; j >= 0; j--)
    if ((c = f_ss[j].E + r_ss[j].E + g) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
}

/* Conquer: recursively around midpoint */

  if (type == 1)
    { c1 = align(A,B,midi,midj,tb,-g,w,g,h,f_str);
      c2 = align(A+midi,B+midj,M-midi,N-midj,-g,te,w,g,h,f_str);
    }
  else
    { align(A,B,midi-1,midj,tb,0,w,g,h,f_str);
      INS(2);
      align(A+midi+1,B+midj,M-midi-1,N-midj,0,te,w,g,h,f_str);
    }
  return midc;
}

/* Interface and top level of comparator */

int ALIGN(unsigned char *A, unsigned char *B, int M, int N,
	  int **W, int G, int H, int *S, int *NC, struct f_struct *f_str)
{ 
  int c, ck;

  sapp = S;
  last = 0;

  c = align(A,B,M,N,-G,-G,W,G,H,f_str);	/* OK, do it */
  ck = CHECK_SCORE(A,B,M,N,S,W,G,H,NC);
  if (c != ck) printf("Check_score error. %d != %d\n",c,ck);
  return c;
}

/* Alignment display routine */

static char ALINE[51], BLINE[51], CLINE[51];

void DISPLAY(unsigned char *A, unsigned char *B, int M, int N,
	    int *S, int AP, int BP, char *sq)
{ register char *a, *b, *c;
  register int   i,  j, op;
           int   lines, ap, bp;

  i = j = op = lines = 0;
  ap = AP;
  bp = BP;
  a = ALINE;
  b = BLINE;
  c = CLINE;
  while (i < M || j < N)
    { if (op == 0 && *S == 0)
        { op = *S++;
          *a = sq[A[++i]];
          *b = sq[B[++j]];
          *c++ = (*a++ == *b++) ? '|' : ' ';
        }
      else
        { if (op == 0)
            op = *S++;
          if (op > 0)
            { *a++ = ' ';
              *b++ = sq[B[++j]];
              op--;
            }
          else
            { *a++ = sq[A[++i]];
              *b++ = ' ';
              op++;
            }
          *c++ = '-';
        }
      if (a >= ALINE+50 || i >= M && j >= N)
        { *a = *b = *c = '\0';
          printf("\n%5d ",50*lines++);
          for (b = ALINE+10; b <= a; b += 10)
            printf("    .    :");
          if (b <= a+5)
            printf("    .");
          printf("\n%5d %s\n      %s\n%5d %s\n",ap,ALINE,CLINE,bp,BLINE);
	  ap = AP + i;
	  bp = BP + j;
          a = ALINE;
          b = BLINE;
          c = CLINE;
        }
    }
}

/* CHECK_SCORE - return the score of the alignment stored in S */

static int CHECK_SCORE(unsigned char *A, unsigned char *B,
		       int M, int N, int *S, int **w, int g, int h, int *NC)
{ 
  register int   i,  j, op, nc;
  int score;

  score = i = j = op = nc = 0;
  while (i < M || j < N) {
    op = *S++;
    if (op == 0) {
      score = w[A[++i]][B[++j]] + score;
      nc++;
      /*      fprintf(stderr,"=%4d %4d %4d %4d\n",i,j,w[A[i]][B[i]],score); */
    }
    else if (op > 0) {
      score = score - (g+op*h);
      /*      fprintf(stderr,">%4d %4d %4d %4d\n",i,j,-(g+op*h),score); */
      i += op;
      nc += op;
    } else {
      score = score - (g-op*h);
      /*      fprintf(stderr,"<%4d %4d %4d %4d\n",i,j,-(g-op*h),score); */
      j -= op;
      nc -= op;
    }
  }
  *NC = nc;
  return(score);
}


int calcons(unsigned char *aa0, int n0,
	    unsigned char *aa1, int n1,
	    int *res, int nres, int *nc,
	    struct a_struct *aln, struct pstruct pst,
	    char *seqc0, char *seqc1,
	    struct f_struct *f_str)
{
  int i0, i1;
  int op, lenc, nd, ns, itmp;
  char *sp0, *sp1;
  int *rp;
  
  /* first fill in the ends */
  aln->min0--; aln->min1--;

  /* #define LFASTA */
#ifndef LFASTA
  if (min(aln->min0,aln->min1)<aln->llen || aln->showall==1)     /* will we show all the start ?*/
    if (aln->min0>=aln->min1) {                        /* aa0 extends more to left */
      aln->smins=0;
      if (aln->showall==1) aln->mins=aln->min0;
      else aln->mins = min(aln->min0,aln->llcntx);
      aancpy(seqc0,(char *)aa0+aln->min0-aln->mins,aln->mins,pst);
      aln->smin0 = aln->min0-aln->mins;
      if ((aln->mins-aln->min1)>0) {
	memset(seqc1,' ',aln->mins-aln->min1);
	aancpy(seqc1+aln->mins-aln->min1,(char *)aa1,aln->min1,pst);
	aln->smin1 = 0;
      }
      else {
	aancpy(seqc1,(char *)aa1+aln->min1-aln->mins,aln->mins,pst);
	aln->smin1 = aln->min1-aln->mins;
      }
    }
    else {
      aln->smins=0;
      if (aln->showall == 1) aln->mins=aln->min1;
      else aln->mins = min(aln->min1,aln->llcntx);
      aancpy(seqc1,(char *)(aa1+aln->min1-aln->mins),aln->mins,pst);
      aln->smin1 = aln->min1-aln->mins;
      if ((aln->mins-aln->min0)>0) {
	memset(seqc0,' ',aln->mins-aln->min0);
	aancpy(seqc0+aln->mins-aln->min0,(char *)aa0,aln->min0,pst);
	aln->smin0 = 0;
      }
      else {
	aancpy(seqc0,(char *)aa0+aln->min0-aln->mins,aln->mins,pst);
	aln->smin0 = aln->min0-aln->mins;
      }
    }
  else {
    aln->mins= min(aln->llcntx,min(aln->min0,aln->min1));
    aln->smins=aln->mins;
    aln->smin0=aln->min0;
    aln->smin1=aln->min1;
    aancpy(seqc0,(char *)aa0+aln->min0-aln->mins,aln->mins,pst);
    aancpy(seqc1,(char *)aa1+aln->min1-aln->mins,aln->mins,pst);
  }
#else
  aln->smin0 = aln->min0;
  aln->smin1 = aln->min1;
  aln->smins = aln->mins = 0;
#endif

/* now get the middle */

  sp0 = seqc0+aln->mins;
  sp1 = seqc1+aln->mins;
  rp = res;
  lenc = aln->nident = aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  i0 = aln->min0;
  i1 = aln->min1;
  
  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      *sp0 = pst.sq[aa0[i0++]];
      *sp1 = pst.sq[aa1[i1++]];
      lenc++;
      if (toupper(*sp0) == toupper(*sp1)) aln->nident++;
      else if (pst.dnaseq==1 && ((*sp0 == 'T' && *sp1 == 'U') ||
				(*sp0=='U' && *sp1=='T')))
	aln->nident++;
      sp0++; sp1++;
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {
	*sp0++ = '-';
	*sp1++ = pst.sq[aa1[i1++]];
	op--;
	lenc++;
	aln->ngap_q++;
      }
      else {
	*sp0++ = pst.sq[aa0[i0++]];
	*sp1++ = '-';
	op++;
	lenc++;
	aln->ngap_l++;
      }
    }
  }

  *nc = lenc;
/*	now we have the middle, get the right end */

#ifndef LFASTA
  /* how much extra to show at end ? */
  if (!aln->llcntx_flg) {
    ns = aln->mins + lenc + aln->llen;	/* show an extra line? */
    ns -= (itmp = ns %aln->llen);	/* itmp = left over on last line */
    if (itmp>aln->llen/2) ns += aln->llen;  /* more than 1/2 , use another*/
    nd = ns - (aln->mins+lenc);		/* this much extra */
  }
  else nd = aln->llcntx;

  if (nd > max(n0-aln->max0,n1-aln->max1))
    nd = max(n0-aln->max0,n1-aln->max1);
  
  if (aln->showall==1) {
    nd = max(n0-aln->max0,n1-aln->max1);	/* reset for showall=1 */
    /* get right end */
    aancpy(seqc0+aln->mins+lenc,(char *)aa0+aln->max0,n0-aln->max0,pst);
    aancpy(seqc1+aln->mins+lenc,(char *)aa1+aln->max1,n1-aln->max1,pst);
    /* fill with blanks - this is required to use one 'nc' */
    memset(seqc0+aln->mins+lenc+n0-aln->max0,' ',nd-(n0-aln->max0));
    memset(seqc1+aln->mins+lenc+n1-aln->max1,' ',nd-(n1-aln->max1));
  }
  else {
    if ((nd-(n0-aln->max0))>0) {
      aancpy(seqc0+aln->mins+lenc,(char *)aa0+aln->max0,n0-aln->max0,pst);
      memset(seqc0+aln->mins+lenc+n0-aln->max0,' ',nd-(n0-aln->max0));
    }
    else aancpy(seqc0+aln->mins+lenc,(char *)aa0+aln->max0,nd,pst);

    if ((nd-(n1-aln->max1))>0) {
      aancpy(seqc1+aln->mins+lenc,(char *)aa1+aln->max1,n1-aln->max1,pst);
      memset(seqc1+aln->mins+lenc+n1-aln->max1,' ',nd-(n1-aln->max1));
    }
    else aancpy(seqc1+aln->mins+lenc,(char *)aa1+aln->max1,nd,pst);
  }
  
#else	/* LFASTA */
  nd = 0;
#endif
  /* #undef LFASTA */
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

int calc_id(const unsigned char *aa0, const int n0,
	    const unsigned char *aa1, const int n1,
	    int *res, int nres,
	    struct a_struct *aln, struct pstruct pst,
	    struct f_struct *f_str)
{
  int i0, i1, nn1;
  int op, lenc, nd, ns, itmp;
  int sp0, sp1;
  const unsigned char *aa1p;
  int *rp;
  
#ifndef TFASTA
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* back up one step */
  aln->min0--; aln->min1--;

  rp = res;
  lenc = aln->nident = aln->ngap_q = aln->ngap_l = aln->nfs =op = 0;
  i0 = aln->min0;
  i1 = aln->min1;

  while (i0 < aln->max0 || i1 < aln->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      lenc++;
      sp0 = toupper(pst.sq[aa0[i0++]]);
      sp1 = toupper(pst.sq[aa1p[i1++]]);
      if (sp0 == sp1) aln->nident++;
      else if (pst.dnaseq==1 &&
	       ((sp0=='T' && sp1== 'U')||(sp0=='U' && sp1=='T')))
	aln->nident++;
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {op--; lenc++; i1++; aln->ngap_q++;}
      else {op++; lenc++; i0++;	aln->ngap_l++;}
    }
  }
  return lenc;
}

#ifdef PCOMPLIB
#include "p_mw.h"

void
update_params(struct qmng_str *qm_msg, struct pstruct *ppst)
{
  ppst->n0 = qm_msg->n0;
}
#endif
