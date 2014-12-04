/*	initsw.c	general and function-specific initializations */

/* $Name: fa34t20b3 $ - $Id: initsw.c,v 1.44 2002/10/02 19:20:18 wrp Exp $ */

/* copyright (c) 1996, 1997, 1998  William R. Pearson and the U. of Virginia */

/* init??.c files provide function specific initializations */

/* h_init()	- called from comp_lib.c, comp_thr.c to initialize pstruct ppst
   		  which includes the alphabet, and pam matrix

   alloc_pam()	- allocate pam matrix space
   initpam2()	- convert from 1D to 2D pam

   f_initenv()	- set up mngmsg and pstruct defaults
   f_getopt()	- read fasta specific command line options
   f_getarg()	- does nothing for SW - gets # of shuffles for prss3

   resetp()	- reset the parameters, scoring matrix for DNA-DNA

   query_parm()	- does nothing for Smith-Waterman, gets number of shuffles,
   		  window shuffle for prss3.

   last_init()	- some things must be done last

   f_initpam()	- set some parameters based on the pam matrix

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef UNIX
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "defs.h"
#include "param.h"
#ifndef PCOMPLIB
#include "mw.h"
#else
#include "p_mw.h"
#endif
#include "structs.h"

#define XTERNAL
#include "upam.h"
#include "uascii.h"
#undef XTERNAL

char *iprompt1=" query sequence file name: ";
#ifndef PRSS
char *refstr="\nPlease cite:\n T. F. Smith and M. S. Waterman, (1981) J. Mol. Biol. 147:195-197; \n W.R. Pearson (1991) Genomics 11:635-650\n";
char *iprompt0=" SSEARCH searches a sequence database\n\
 using the Smith-Waterman algorithm\n";
char *iprompt2=" database file name: ";
int PgmDID=404;
#else
char *refstr="\nPlease cite:\n W.R. Pearson (1996) Meth. Enzymol. 266:227-258\n";
char *iprompt0=" PRSS compares a query sequence to shuffled sequences\n\
 using the Smith-Waterman algorithm\n";
char *iprompt2=" shuffle file name: ";
int PgmDID=400;
#ifndef SS
#define SS
#endif
#endif


char *verstr="version 3.4t20 Sep 26, 2002";

char *prog_func = "Smith-Waterman";
char *s_optstr = "3E:F:f:g:";

extern void init_ascii(int ext_sq, int *sascii, int dnaseq);

/* Sets defaults assuming a protein sequence */
void h_init (struct pstruct *ppst, char *pgm_name)
{
  int i;

#ifndef SS  
  strncpy(pgm_name, "gsw", MAX_FN);
#else
  strncpy(pgm_name, "ssw", MAX_FN);
#endif

  standard_pam("BL50",ppst);

  ppst->nsq = naa;
  ppst->nsqx = naax;
  for (i=0; i<=ppst->nsqx; i++) {
    ppst->sq[i]=aa[i];		/* sq = aa */
    ppst->hsq[i]=haa[i];	/* hsq = haa */
    ppst->sqx[i]=aax[i];	/* sq = aax */
    ppst->hsqx[i]=haax[i];	/* hsq = haax */
  }
  memcpy(qascii,aascii,sizeof(qascii));

  /* set up the c_nt[] mapping */
  ppst->c_nt[0]=0;
  for (i=1; i<=nnt; i++) {
    ppst->c_nt[i]=gc_nt[i];
    ppst->c_nt[i+nnt]=gc_nt[i]+nnt;
  }
}

/*
 * alloc_pam(): allocates memory for the 2D pam matrix as well
 * as for the integer array used to transmit the pam matrix
 */
void    alloc_pam (int d1, int d2, struct pstruct *ppst)
{
  int     i, *d2p;
  char err_str[128];

  if ((ppst->pam2[0] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
  }

  if ((ppst->pam2[1] = (int **) malloc (d1 * sizeof (int *))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
  }

  if ((pam12 = (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
   }

   if ((pam12x = (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2d pam matrix: %d",d2);
     s_abort (err_str,"");
   }

   for (i = 0, d2p = pam12; i < d1; i++, d2p += d2)
      ppst->pam2[0][i] = d2p;

   for (i = 0, d2p = pam12x; i < d1; i++, d2p += d2)
      ppst->pam2[1][i] = d2p;
}

/*
 *  initpam2(struct pstruct pst): Converts 1-D pam matrix to 2-D
 */
void    initpam2 (struct pstruct *ppst)
{
   int     i, j, k, nsq, pam_x, sa_x;


   nsq = ppst->nsq;
   sa_x = pascii['X'];
   ppst->pam2[0][0][0] = -BIGNUM;
   ppst->pam_h = -1; ppst->pam_l = 1;

   k = 0;
   for (i = 1; i <= nsq; i++) {
     ppst->pam2[0][0][i] = ppst->pam2[0][i][0] = -BIGNUM;
     for (j = 1; j <= i; j++) {
       ppst->pam2[0][j][i] = ppst->pam2[0][i][j] = pam[k++] - ppst->pamoff;
       if (ppst->pam_l > ppst->pam2[0][i][j]) ppst->pam_l =ppst->pam2[0][i][j];
       if (ppst->pam_h < ppst->pam2[0][i][j]) ppst->pam_h =ppst->pam2[0][i][j];
     }
   }

   if (ppst->pam_x <= 0) {
     ppst->pam2[0][sa_x][sa_x]=ppst->pam_x;
     for (i=1; i<=nsq; i++)
       ppst->pam2[0][sa_x][i] = ppst->pam2[0][i][sa_x]=ppst->pam_x;
   }
   else {ppst->pam_x = ppst->pam2[0][sa_x][sa_x];}

   pam_x = ppst->pam_x;

   if (ppst->ext_sq_set) {	/* using extended alphabet */
     /* fill in pam2[1] matrix */
     ppst->pam2[1][0][0] = -BIGNUM;
     /* fill in additional parts of the matrix */
     for (i = 1; i <= nsq; i++) {

       /* -BIGNUM to all matches vs 0 */
       ppst->pam2[0][0][i+nsq] = ppst->pam2[0][i+nsq][0] = 
       ppst->pam2[1][0][i+nsq] = ppst->pam2[1][i+nsq][0] = 
	 ppst->pam2[1][0][i] = ppst->pam2[1][i][0] = -BIGNUM;

       for (j = 1; j <= nsq; j++) {

	 /* replicate pam2[0] to i+nsq, j+nsq */
	 ppst->pam2[0][i+nsq][j] = ppst->pam2[0][i][j+nsq] =
	   ppst->pam2[0][i+nsq][j+nsq] = ppst->pam2[1][i][j] =
	   ppst->pam2[0][i][j];

	 /* set the high portion of pam2[1] to the corresponding value
            of pam2[1][sa_x][j] */

	 ppst->pam2[1][i+nsq][j] = ppst->pam2[1][i][j+nsq]=
	   ppst->pam2[1][i+nsq][j+nsq]=ppst->pam2[0][sa_x][j];
       }
     }
   }
}

/*
 * initenv () : initializes the environment
 */
void
f_initenv (struct mngmsg *m_msg, struct pstruct *ppst, unsigned char  **aa0) {

#ifdef PRSS
  m_msg->aln.llen = 0;
#endif
  m_msg->stages = 1;
  m_msg->last_calc_flg = 0;

  /* seq_type not defined for general FASTA, SSEARCH */
  m_msg->qdnaseq = SEQTYPE_UNK;

  strncpy (m_msg->label, "s-w", sizeof(m_msg->label));
   strncpy (m_msg->f_id0,"sw",3);
   strncpy (m_msg->f_id1,"sw",3);
   ppst->histint = 2;

   /* gap values are now set by standard_pam() */

   ppst->ggapval = -2;
   m_msg->nframe = -1;
   m_msg->qframe = 1;
   ppst->score_ix = 0;
   ppst->sw_flag = 1;
   m_msg->nrelv = 1;		/* number of relevant scores */
   m_msg->srelv = 1;		/* relevant scores in showbest */
   m_msg->arelv = 1;		/* number of alignment scores */
   strncpy (m_msg->alab[0],"s-w opt",20);

   alloc_pam (MAXSQ, MAXSQ, ppst);
}

static int gap_set=0;
static int del_set=0;
static int nframe_set = 0;

void    f_getopt (copt, optarg, m_msg, ppst)
char    copt;
char   *optarg;
struct mngmsg *m_msg;
struct pstruct *ppst;
{
  switch (copt) {
  case '3':
    m_msg->qframe = m_msg->nframe = 1; 
    nframe_set = 1;
    break;
  case 'E': sscanf(optarg,"%lf",&m_msg->e_cut);
    m_msg->e_cut_set = 1;
    break;
  case 'F': sscanf(optarg,"%lf",&m_msg->e_low);
    m_msg->e_cut_set = 1;
    break;
  case 'f':
    sscanf (optarg, "%d", &ppst->gdelval);
    if (ppst->gdelval > 0) ppst->gdelval = -ppst->gdelval;
    del_set = 1;
    break;
  case 'g':
    sscanf (optarg, "%d", &ppst->ggapval);
    if (ppst->ggapval > 0) ppst->ggapval = -ppst->ggapval;
    gap_set = 1;
    break;
  }
}

void
f_lastenv(struct mngmsg *m_msg, struct pstruct *ppst)
{
  if (m_msg->qdnaseq == SEQTYPE_UNK)
    build_xascii(qascii);
  else
   init_ascii(ppst->ext_sq_set,qascii,m_msg->qdnaseq);
}

void
f_getarg (int argc, char **argv, int optind, 
	  struct mngmsg *m_msg, struct pstruct *ppst)
{
  /* this sets up the number of shuffles for PRSS */
#ifdef PRSS
  if (argc - optind == 4) {
      sscanf (argv[optind + 3], "%d", &m_msg->mshow);
      m_msg->mshow_flg = 1;
  }
#endif
}

/* recode has become function specific to accommodate FASTS/M */
int
recode(unsigned char *seq, int n, int *qascii) {
  int i,j;

  for (i=0,j= 0; j < n; j++) {
    if ((seq[i] = qascii[seq[j]]) < NA) i++;
  }
  return i;
}

void
resetp (struct mngmsg *m_msg, struct pstruct *ppst)
{
  int i;

  if (m_msg->qdnaseq == SEQTYPE_DNA) {
    ppst->dnaseq = SEQTYPE_DNA;
    pascii=nascii;

    m_msg->ldnaseq = SEQTYPE_DNA;
    memcpy(lascii,nascii,sizeof(lascii));	/* initialize lib mapping */

    if (!nframe_set) {
#ifndef PRSS
      m_msg->qframe = 2;		/* number of query seq directions */
#else
      m_msg->qframe = 1;		/* number of query seq directions */
#endif
      m_msg->nframe = 1;		/* use frames 0, 1 */
    }
    else {
      m_msg->qframe = m_msg->revcomp+1;
    }
  }
  else {
    ppst->dnaseq = SEQTYPE_PROT;
    pascii=aascii;

    m_msg->ldnaseq = SEQTYPE_PROT;
    memcpy(lascii,aascii,sizeof(lascii));	/* initialize lib mapping */
  }

  /* set extended alphabet */
  init_ascii(ppst->ext_sq_set,qascii,m_msg->qdnaseq);
  init_ascii(ppst->ext_sq_set,lascii,m_msg->ldnaseq);

  if (m_msg->qdnaseq == SEQTYPE_DNA) {

    ppst->histint = 4;
#ifndef GAP_OPEN
    if (!del_set) ppst->gdelval = -16;	/* def. del penalty */
#else
    if (!del_set) ppst->gdelval = -12;	/* def. del penalty */
#endif
    if (!gap_set) ppst->ggapval = -4;

    if (!m_msg->e_cut_set) m_msg->e_cut=2.0;

    ppst->nsq = nnt;
    ppst->nsqx = nntx;

    for (i=0; i<=ppst->nsqx; i++) {
      ppst->hsq[i] = hnt[i];
      ppst->sq[i] = nt[i];
      ppst->hsqx[i] = hntx[i];
      ppst->sqx[i] = ntx[i];
    }

    if (!ppst->pam_set) {
      if (ppst->p_d_set)
	mk_n_pam(npam,nnt,ppst->p_d_mat,ppst->p_d_mis);
      else if (strncmp(ppst->pamfile,"BL50",4)==0) {
	strncpy (ppst->pamfile, "+5/-4", sizeof(ppst->pamfile));
      }
      pam = npam;
    }

    strncpy (m_msg->sqnam, "nt",sizeof(m_msg->sqnam));
    strncpy (m_msg->sqtype, "DNA",sizeof(m_msg->sqtype));
  }
}
/* query_parm()	this function asks for any additional parameters
	that have not been provided.  Could be null. */

void
query_parm (m_msp, ppst)
struct mngmsg *m_msp;
struct pstruct *ppst;
{
   char    qline[40];

#ifdef PRSS
   if (!m_msp->mshow_flg) m_msp->mshow = 200;

   printf(" number of shuffles [%d]? ",m_msp->mshow);
   fflush(stdout);
   if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
   else sscanf(qline,"%d",&m_msp->mshow);

   if (m_msp->aln.llen < 1) {
     printf (" local (window) (w) or uniform (u) shuffle [u]? ");
     if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
     else if (qline[0]=='w' || qline[0]=='W') {
       m_msp->aln.llen = 20;
       printf(" local shuffle window size [%d]? ",m_msp->aln.llen);
       if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
       else sscanf(qline,"%d",&m_msp->aln.llen);
     }
   }

#endif
}

/* last_init() cannot look at aa0, n0, because it is only run once,
   it is not run before each new aa0 search */
void
last_init (struct mngmsg *m_msg, struct pstruct *ppst
#ifdef PCOMPLIB
	   , int nnodes
#endif
	   )
{
  double *kar_p;
  double aa0_f[MAXSQ];

  m_msg->qshuffle = 0;
  m_msg->escore_flg = 0;
  m_msg->shuff_max = 1;
  m_msg->last_calc_flg = 0;

  m_msg->thr_fact = 1;
#ifdef PRSS
   if (!m_msg->mshow_flg) m_msg->mshow = 200;
   if (ppst->zsflag < 10) ppst->zsflag += 10;
#endif

  m_msg->nitt1 = m_msg->qframe-1;

  initpam2(ppst);

   /* once we have a complete pam matrix, we can calculate Lambda and K 
      for "average" sequences */

   kar_p = NULL;
   init_karlin_a(ppst, aa0_f, &kar_p);
   do_karlin_a(ppst->pam2[0], ppst,aa0_f, 
	       kar_p, &m_msg->Lambda, &m_msg->K, &m_msg->H);
   free(kar_p);
}

void
f_initpam (line, ppst)
char   *line;
struct pstruct *ppst;
{}

#ifndef PCOMPLIB
void
qshuffle() {}
#endif

int
last_calc(unsigned char *aa0, unsigned char *aa1, int maxn,
	  struct beststr **bestp_arr, int nbest,
	  struct mngmsg *m_msg, struct pstruct *pst,
	  void **f_str)
{
  return nbest;
}

void sortbest (bptr, nbest, irelv)
struct beststr **bptr;
int nbest, irelv;
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->score[irelv] >= bptr[j + gap]->score[irelv]) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

void 
show_aux(FILE *fp, struct beststr *bptr) {}

void
last_params(unsigned char *aa0, int n0, 
	   struct mngmsg *m_msg,
	   struct pstruct *ppst
#ifdef PCOMPLIB
	   , struct qmng_str *qm_msg
#endif
	   ) {
  /* needed for scaleswn.c */
  ppst->n0 = m_msg->n0;
  m_msg->last_calc_flg = 0;
  m_msg->qshuffle = 0;
  m_msg->escore_flg = 0;
  m_msg->nm0 = 1;
}


