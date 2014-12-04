/*	initfa.c	*/

/* $Name: fa34t20b3 $ - $Id: initfa.c,v 1.48 2002/10/02 19:20:18 wrp Exp $ */

/* copyright (c) 1996, 1997, 1998  William R. Pearson and the U. of Virginia */

/* init??.c files provide function specific initializations */

/* h_init()	- called from comp_lib.c, comp_thr.c to initialize pstruct ppst
   		  which includes the alphabet, and pam matrix

   alloc_pam()	- allocate pam matrix space
   initpam2()	- convert from 1D to 2D pam

   f_initenv()	- set up mngmsg and pstruct defaults
   f_getopt()	- read fasta specific command line options
   f_getarg()	- read ktup

   resetp()	- reset the parameters, scoring matrix for DNA-DNA/DNA-prot

   query_parm()	- ask for ktup
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
#include "structs.h"
#include "param.h"

#ifndef PCOMPLIB
#include "mw.h"
#else
#include "p_mw.h"
#endif

#define XTERNAL
#include "upam.h"
#include "uascii.h"
#undef XTERNAL

#define MAXWINDOW 32

#ifdef FASTA
char *prog_func = "FASTA";
char *iprompt0=" FASTA searches a protein or DNA sequence data bank\n";
int PgmDID=401;
#endif
#ifdef TFASTX
char *prog_func = "TFASTX";
char *iprompt0=" TFASTX compares a protein to a translated DNA data bank\n";
int PgmDID=406;
#endif
#ifdef TFASTY
char *prog_func = "TFASTY";
char *iprompt0=" TFASTY compares a protein to a translated DNA data bank\n";
int PgmDID=406;
#endif
#ifdef FASTX
char *prog_func = "FASTX";
char *iprompt0=" FASTX compares a DNA sequence to a protein sequence data bank\n";
int PgmDID=405;
#endif
#ifdef FASTY
char *prog_func = "FASTY";
char *iprompt0=" FASTY compares a DNA sequence to a protein sequence data bank\n";
int PgmDID=405;
#endif
#ifdef FASTF
char *prog_func = "FASTF";
char *iprompt0=" FASTF compares mixed peptides to a protein databank\n";
int PgmDID=400;
#endif
#ifdef TFASTF
char *prog_func = "TFASTF";
char *iprompt0=" TFASTF compares mixed peptides to a translated DNA data bank\n";
int PgmDID=400;
#endif
#ifdef TFASTA
char *prog_func = "TFASTA";
char *iprompt0=" TFASTA translates and searches a DNA sequence data bank\n";
int PgmDID=402;
#endif
#ifdef FASTS
char *prog_func = "FASTS";
char *iprompt0=" FASTS compares linked peptides to a protein databank\n";
int PgmDID=400;
#endif
#ifdef TFASTS
char *prog_func = "TFASTS";
char *iprompt0=" TFASTS compares linked peptides to a protein databank\n";
int PgmDID=400;
#endif
#ifdef FASTM
char *prog_func = "FASTM";
char *iprompt0=" FASTM compares ordered peptides to a protein databank\n";
int PgmDID=400;
#endif

char *iprompt1=" test sequence file name: ";
char *iprompt2=" database file name: ";
#if defined(FASTX) || defined(TFASTX) || defined(FASTY) || defined(TFASTY)
char *refstr="\nPlease cite:\n Pearson et al, Genomics (1997) 46:24-36\n";
#define NOT_FASTA
#endif
#if defined(FASTS) || defined(TFASTS) || defined(FASTF) || defined(TFASTF) || defined(FASTM)
char *refstr="\nPlease cite:\n Mackey et al. Mol. Cell. Proteomics  (2002) 1:139-147\n";
#define NOT_FASTA
#endif
#ifndef NOT_FASTA
char *refstr="\nPlease cite:\n W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448\n";
#endif

char *verstr="version 3.4t20 Sep 26, 2002";

char   *s_optstr = "13Ac:E:f:F:g:h:j:oy:";

static int mktup=2;

extern int max_workers;

extern void s_abort(char *, char *);
extern void init_ascii(int ext_sq, int *sascii, int dnaseq);
extern int standard_pam(char *smstr, struct pstruct *ppst);
extern void mk_n_pam(int *arr,int siz, int mat, int mis);
extern void init_karlin_a(struct pstruct *, double *, double **);
extern int do_karlin_a(int **, struct pstruct *, double *,
		       double *, double *, double *, double *);

#if defined(TFAST) || defined(FASTX) || defined(FASTY)
extern void aainit(int tr_type, int debug);
#endif

/* Sets defaults assuming a protein sequence */
void h_init (struct pstruct *ppst, char *pgm_name )
{
  int i;

#ifdef FASTA	/* FASTA */
   strncpy (pgm_name,"fa",MAX_FN);
#endif
#ifdef FASTX		/* FASTX */
   strncpy (pgm_name,"fx",MAX_FN);
#endif
#ifdef TFASTX		/* TFASTX */
   strncpy (pgm_name,"tfx",MAX_FN);
#endif
#ifdef FASTY		/* FASTX */
   strncpy (pgm_name,"fy",MAX_FN);
#endif
#ifdef TFASTY		/* TFASTX */
   strncpy (pgm_name,"tfy",MAX_FN);
#endif
#ifdef TFASTA		/* TFASTA */
   strncpy (pgm_name,"tfa",MAX_FN);
#endif
#ifdef FASTF		/* FASTF */
   strncpy (pgm_name,"ff",MAX_FN);
#endif
#ifdef TFASTF		/* TFASTF */
   strncpy (pgm_name,"tff",MAX_FN);
#endif
#ifdef FASTS		/* FASTS */
   strncpy (pgm_name,"fs",MAX_FN);
#endif
#ifdef TFASTS		/* FASTS */
   strncpy (pgm_name,"tfs",MAX_FN);
#endif
#ifdef FASTM		/* FASTS */
   strncpy (pgm_name,"fm",MAX_FN);
#endif

#if !defined(FASTX) && !defined(FASTY)
  memcpy(qascii,aascii,sizeof(qascii));
#else
  memcpy(qascii,nascii,sizeof(qascii));
#endif

#if defined(FASTF) || defined(TFASTF) || defined(FASTS) || defined(TFASTS) || defined(FASTM)
  qascii[','] = ESS;

#if defined(FASTF) || defined(FASTS) || defined(FASTM)
  standard_pam("M20",ppst);
#else		/* TFASTF */
  standard_pam("M10",ppst);
#endif
#else
  standard_pam("BL50",ppst);
#endif  

  ppst->nsq = naa;
  ppst->nsqx = naax;
  for (i=0; i<=ppst->nsqx; i++) {
    ppst->sq[i] = aa[i];
    ppst->hsq[i] = haa[i];
    ppst->sqx[i]=aax[i];	/* sq = aa */
    ppst->hsqx[i]=haax[i];	/* hsq = haa */
  }

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
void
alloc_pam (int d1, int d2, struct pstruct *ppst)
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

  if ((d2p = pam12 = (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2D pam matrix: %d",d1);
     s_abort (err_str,"");
   }

   for (i = 0; i < d1; i++, d2p += d2)
      ppst->pam2[0][i] = d2p;

   if ((d2p=pam12x= (int *) malloc (d1 * d2 * sizeof (int))) == NULL) {
     sprintf(err_str,"Cannot allocate 2d pam matrix: %d",d2);
     s_abort (err_str,"");
   }

   for (i = 0;  i < d1; i++, d2p += d2)
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
#if defined(FASTS) || defined(TFASTS) || defined(FASTF) || defined(TFASTF) || defined(FASTM)
   ppst->pam_x = 0;
#endif

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

/*  function specific initializations */
void
f_initenv (struct mngmsg *m_msg, struct pstruct *ppst, unsigned char **aa0) {

  m_msg->stages = 2;
  m_msg->last_calc_flg=0;

  /* seq_type not defined for general FASTA, SSEARCH */
  m_msg->qdnaseq = SEQTYPE_UNK;
#if defined(TFAST)
  m_msg->qdnaseq = SEQTYPE_PROT;	/* must be protein */
#endif
#if defined(FASTF)
  m_msg->qdnaseq = SEQTYPE_PROT;	/* must be protein */
#else
#endif
#if defined(FASTX) || defined(FASTY)
  m_msg->qdnaseq = SEQTYPE_DNA;	/* must be DNA */
#endif

#ifdef FASTA	/* FASTA */
   strncpy (m_msg->f_id0,"fa",3);
   strncpy (m_msg->f_id1,"sw",3);
#endif
#ifdef FASTX		/* FASTX */
   strncpy(m_msg->f_id0,"fx",3);
   strncpy(m_msg->f_id1,"sx",3);
#endif
#ifdef FASTY		/* FASTY */
   strncpy(m_msg->f_id0,"fy",3);
   strncpy(m_msg->f_id1,"sy",3);
#endif
#ifdef TFASTX		/* TFASTX */
   strncpy (m_msg->f_id0,"tx",3);
   strncpy (m_msg->f_id1,"sx",3);
#endif
#ifdef TFASTY		/* TFASTX */
   strncpy (m_msg->f_id0,"ty",3);
   strncpy (m_msg->f_id1,"sy",3);
#endif
#ifdef TFASTA		/* TFASTA */
   strncpy (m_msg->f_id0,"tf",3);
   strncpy (m_msg->f_id1,"tf",3);
#endif
#ifdef TFASTF		/* TFASTF */
   strncpy (m_msg->f_id0,"ft",3);
   strncpy (m_msg->f_id1,"ft",3);
#endif
#ifdef FASTF		/* FASTF */
   strncpy (m_msg->f_id0,"ff",3);
   strncpy (m_msg->f_id1,"ff",3);
#endif
#ifdef FASTS		/* FASTF */
   strncpy (m_msg->f_id0,"fs",3);
   strncpy (m_msg->f_id1,"fs",3);
#endif
#ifdef TFASTS		/* FASTF */
   strncpy (m_msg->f_id0,"ts",3);
   strncpy (m_msg->f_id1,"ts",3);
#endif

   strncpy (m_msg->alab[0],"initn",20);
   strncpy (m_msg->alab[1],"init1",20);
   strncpy (m_msg->alab[2],"opt",20);

   ppst->score_ix = 0;
   ppst->histint = 2;
   m_msg->qframe = 1;
#if !defined(TFASTX) && !defined(TFASTY) && !defined(TFASTF) && !defined(TFASTS) && !defined(TFASTA)
   ppst->sw_flag = 1;
   m_msg->nframe = -1;
#else
   if (strncmp(ppst->pamfile,"BL50",4)==0) {
#ifndef GAP_OPEN
     ppst->gdelval = -16;
#else
     ppst->gdelval = -14;
#endif
     ppst->ggapval = -2;
   }
#if defined(TFASTX) || defined(TFASTY)
   m_msg->nframe = 2;
   ppst->sw_flag = 1;
#else
   m_msg->nframe = 6;
   ppst->sw_flag = 0;
#endif
#endif
   m_msg->nrelv = 3;		/* number of relevant scores */
#if !defined(FASTF) && !defined(TFASTF) && !defined(FASTS) && !defined(TFASTS)
   m_msg->srelv = 1;		/* relevant scores in showbest */
   strncpy (m_msg->label, " opt", sizeof(m_msg->label));
#else
   m_msg->srelv = 2;
   strncpy (m_msg->label, "initn init1", sizeof(m_msg->label));
#endif
   m_msg->arelv = 3;		/* number of relevant scores */

/* see param.h for the definition of all these */

   m_msg->qshuffle = 0;
   mktup = 2;
   ppst->param_u.fa.bestscale = 300;
   ppst->param_u.fa.bestoff = 36;
   ppst->param_u.fa.bkfact = 6;
   ppst->param_u.fa.scfact = 3;
   ppst->param_u.fa.bktup = 2;
   ppst->param_u.fa.ktup = 0;
   ppst->param_u.fa.bestmax = 50;
   ppst->param_u.fa.pamfact = 1;
   ppst->param_u.fa.altflag = 0;
   ppst->param_u.fa.optflag = 1;
   ppst->param_u.fa.iniflag = 0;
   ppst->param_u.fa.optcut = 0;
   ppst->param_u.fa.optcut_set = 0;
   ppst->param_u.fa.cgap = 0;
   ppst->param_u.fa.optwid = MAXWINDOW;
   
   alloc_pam (MAXSQ, MAXSQ, ppst);
}

/*  switches for fasta only */

static int gap_set=0;
static int del_set=0;
static int shift_set=0;
static int subs_set=0;
static int sw_flag_set=0;
static int nframe_set=0;
static int wid_set=0;

void
f_getopt (copt, optarg, m_msg, ppst)
char    copt;
char   *optarg;
struct mngmsg *m_msg;
struct pstruct *ppst;
{
  switch (copt) {
  case '1':
     ppst->param_u.fa.iniflag=1;
     break;
  case '3':
    nframe_set = 1;
#ifdef TFASTA
    m_msg->nframe = 3; break;
#else
    m_msg->nframe = 1;	/* for TFASTXY */
    m_msg->qframe = 1;  /* for FASTA, FASTX */
    break;
#endif
  case 'A':
    ppst->sw_flag= 1;
    sw_flag_set = 1;
    break;
  case 'c':
    sscanf (optarg, "%d", &ppst->param_u.fa.optcut);
    ppst->param_u.fa.optcut_set = 1;
    break;
  case 'E':
    sscanf(optarg,"%lf",&m_msg->e_cut);
    m_msg->e_cut_set = 1;
    break;
  case 'F':
    sscanf(optarg,"%lg",&m_msg->e_low);
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
  case 'h':
    sscanf (optarg, "%d", &ppst->gshift);
    if (ppst->gshift > 0) ppst->gshift = -ppst->gshift;
    shift_set = 1;
    break;
  case 'j':
    sscanf (optarg, "%d", &ppst->gsubs);
    subs_set = 1;
    break;
  case 'o':
    ppst->param_u.fa.optflag = 0;
    m_msg->nrelv = 2;
    break;
  case 'y':
    sscanf (optarg, "%d", &ppst->param_u.fa.optwid);
    wid_set = 1;
    break;
  }
}

void
f_lastenv (struct mngmsg *m_msg, struct pstruct *ppst)
{
  if (m_msg->qdnaseq == SEQTYPE_UNK)
    build_xascii(qascii);

/* this check allows lc DNA sequence queries with FASTX */
#if defined(FASTA) || defined(TFASTA) && !defined(FASTS) && !defined(FASTM) && !defined(FASTF)
  else
   init_ascii(ppst->ext_sq_set,qascii,m_msg->qdnaseq);
#endif
}

void
f_getarg (int argc, char **argv, int optind,
	  struct mngmsg *m_msg, struct pstruct *ppst)
{
#if !defined(FASTF) && !defined (TFASTF)
   if (argc - optind == 4)
      sscanf (argv[optind + 3], "%d", &ppst->param_u.fa.ktup);
   else
#endif
      ppst->param_u.fa.ktup = -ppst->param_u.fa.bktup;

}

/* recode has become function specific to accommodate FASTS/M */
int
recode(unsigned char *seq, int n, int *qascii) {
  int i,j;

#if defined(FASTS) || defined(FASTM)
  qascii[',']=ESS;
#endif

  for (i=0,j= 0; j < n; j++) {
    if ((seq[i] = qascii[seq[j]]) < NA) i++;
  }
  
  return i;
}

/* here we have the query sequence, all the command line options,
   but we need to set various parameter options based on the type
   of the query sequence (m_msg->qdnaseq = 0:protein/1:DNA) and
   the function (FASTA/FASTX/TFASTA)
*/

#if !defined(FASTX) && !defined(FASTY)

/* this resetp is for conventional a FASTA/TFASTXYZ search */
void
resetp (struct mngmsg *m_msg, struct pstruct *ppst)
{
  int i;

#if defined(TFAST)
  if (m_msg->qdnaseq==SEQTYPE_DNA) {
    fprintf(stderr," %s compares a protein to a translated\n\
DNA sequence library.  Do not use a DNA query/scoring matrix.\n",prog_func);
    exit(1);
  }

  ppst->dnaseq = 0;	/* force to protein */
  pascii=aascii;

  m_msg->ldnaseq = SEQTYPE_DNA;
  memcpy(lascii,nascii,sizeof(lascii));	/* initialize lib mapping */
  /* no init_ascii() because we translate lower case library sequences */
#else		/* !defined(TFAST) - normal aa:aa or dna:dna */
  if (m_msg->qdnaseq == SEQTYPE_DNA) {
    ppst->dnaseq = 1;
    pascii=&nascii[0];

    m_msg->ldnaseq = SEQTYPE_DNA;
    memcpy(lascii,nascii,sizeof(lascii));	/* initialize lib mapping */

    if (!nframe_set) {
      m_msg->qframe = 2;		/* number of query seq directions */
      m_msg->nframe = 1;		/* use frames 0, 1 */
    }
    else m_msg->qframe = m_msg->revcomp+1;
  }
  else {
    ppst->dnaseq = 0;
    pascii=&aascii[0];

    m_msg->ldnaseq = SEQTYPE_PROT;
    memcpy(lascii,aascii,sizeof(lascii));	/* initialize lib mapping */
  }
  /* set extended alphabet */
  init_ascii(ppst->ext_sq_set,lascii,m_msg->ldnaseq);
#endif

  /* change settings for DNA search */
  if (m_msg->qdnaseq == SEQTYPE_DNA) {

    if (!wid_set) ppst->param_u.fa.optwid = 16;
    ppst->histint = 4;
    if (!del_set) {
#ifndef GAP_OPEN
      ppst->gdelval = -16;	/* def. del penalty */
#else
      ppst->gdelval = -12;	/* def. open penalty */
#endif
    }
    if (!gap_set) ppst->ggapval = -4;	/* def. gap penalty */

    /* these parameters are used to scale optcut, they should be replaced
       by statistically based parameters */
    ppst->param_u.fa.bestscale = 80;
    ppst->param_u.fa.bkfact = 5;
    ppst->param_u.fa.scfact = 1;
    ppst->param_u.fa.bktup = 6;
    ppst->param_u.fa.bestmax = 80;
    ppst->param_u.fa.bestoff = 45;

    /* largest ktup */
    mktup = 6;

    /* use Ecut = 2.0 for DNA */
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
#if !defined(FASTS) && !defined(FASTM)
      else if (strncmp(ppst->pamfile,"BL50",4)==0) {
	strncpy (ppst->pamfile, "+5/-4", sizeof(ppst->pamfile));
#else
      else if (strncmp(ppst->pamfile,"MD_20",5)==0) {
	strncpy (ppst->pamfile, "+2/-2", sizeof(ppst->pamfile));
	ppst->p_d_mat = +2;
	ppst->p_d_mis = -2;
	mk_n_pam(npam,nnt,ppst->p_d_mat,ppst->p_d_mis);
#endif
      }
      pam = npam;
    }

    strncpy (m_msg->sqnam, "nt",sizeof(m_msg->sqnam));
    strncpy (m_msg->sqtype, "DNA",sizeof(m_msg->sqtype));
    if (ppst->param_u.fa.pamfact >= 0)
      ppst->param_u.fa.pamfact = 0;
    if (ppst->param_u.fa.ktup < 0)
      ppst->param_u.fa.ktup = -ppst->param_u.fa.bktup;
    if (!sw_flag_set) ppst->sw_flag = 0;
  }	/* end DNA reset */

#if defined(TFASTX) || defined(TFASTY)
  if (strncmp(ppst->pamfile,"BL50",4)==0) {
    if (!del_set) {
#ifndef GAP_OPEN
      ppst->gdelval = -16;
#else
      ppst->gdelval = -14;
#endif
    }
    if (!gap_set) ppst->ggapval = -2;
  }
  if (!shift_set) ppst->gshift = -20;
  if (!subs_set) {
#ifndef GAP_OPEN
    ppst->gsubs = ppst->gdelval - 4;
#else
    ppst->gsubs = ppst->gdelval + ppst->ggapval - 4;
#endif
  }
  if (!m_msg->e_cut_set) m_msg->e_cut=5.0;
#endif
#if defined(TFASTF)
  if (!m_msg->e_cut_set) m_msg->e_cut=5.0;
#endif
}
#else
/* this resetp is for conventional a FASTXY search */
void
resetp (struct mngmsg *m_msg, struct pstruct *ppst)
{

  if (m_msg->qdnaseq != SEQTYPE_DNA) {
    fprintf(stderr," FASTX/Y compares a DNA sequence to a protein database\n");
    fprintf(stderr," Use a DNA query\n");
    exit(1);
  }

  ppst->dnaseq = 0;	/* protein parameters */
  pascii = &aascii[0];

  m_msg->ldnaseq = SEQTYPE_PROT;	/* protein library */
  memcpy(lascii,aascii,sizeof(lascii));

  init_ascii(ppst->ext_sq_set,lascii,m_msg->ldnaseq);

  if (!nframe_set) {
    m_msg->qframe = 2;	/* two query frames */
    m_msg->nframe = 1;
  }
  if (!wid_set) {
    if (ppst->param_u.fa.ktup==1) ppst->param_u.fa.optwid = 32;
    else ppst->param_u.fa.optwid = 16;
  }
  if (strncmp(ppst->pamfile,"BL50",4)==0) {
    if (!del_set) {
#ifndef GAP_OPEN
      ppst->gdelval = -16;
#else
      ppst->gdelval = -14;
#endif
    }
    if (!gap_set) ppst->ggapval = -2;
  }
  if (!shift_set) ppst->gshift = -20;
  if (!subs_set) {
#ifndef GAP_OPEN
    ppst->gsubs = ppst->gdelval - 4;
#else
    ppst->gsubs = ppst->gdelval + ppst->ggapval - 4;
#endif
  }
  if (!m_msg->e_cut_set) m_msg->e_cut=5.0;

  strncpy (m_msg->sqnam, "aa",sizeof(m_msg->sqnam));
  strncpy (m_msg->sqtype, "protein",sizeof(m_msg->sqtype));
}
#endif

/* query_parm()	this function asks for any additional parameters
	that have not been provided.  Could be null. */

void
query_parm (struct mngmsg *m_msp, struct pstruct *ppst)
{
   char    qline[40];

   if (ppst->param_u.fa.ktup < 0)
      ppst->param_u.fa.ktup = -ppst->param_u.fa.ktup;

#if !defined (FASTF) && !defined(TFASTF)
   if (ppst->param_u.fa.ktup == 0) {
      printf (" ktup? (1 to %d) [%d] ", mktup, ppst->param_u.fa.bktup);
      if (fgets (qline, sizeof(qline), stdin) == NULL) exit (0);
      else sscanf(qline,"%d",&ppst->param_u.fa.ktup);
   }
#endif

   if (ppst->param_u.fa.ktup == 0)
     ppst->param_u.fa.ktup = ppst->param_u.fa.bktup;
}

/* last_init() cannot look at aa0, n0, because it is only run once,
   it is not run before each new aa0 search */
void
last_init (struct mngmsg *m_msg, struct pstruct *ppst
#ifdef PCOMPLIB
	   ,int nnodes
#endif
	   )
{
  int ix_l, ix_i, i;
  double *kar_p;
  double aa0_f[MAXSQ];

#if defined(FASTF) || defined(TFASTF) || defined(FASTS) || defined(TFASTS)
  m_msg->nohist = 1;
  m_msg->shuff_max = 2000;
#ifndef PCOMPLIB
  ppst->shuff_node = m_msg->shuff_max/max_workers;
#else
  ppst->shuff_node = m_msg->shuff_max/nnodes;
#endif
#endif

#ifndef PCOMPLIB
#if defined(FASTX) || defined(FASTY) || defined(TFAST) 
  /* set up translation tables: faatran.c */
  aainit(ppst->tr_type,ppst->debug_lib);
#endif
#endif

/* a sanity check */
#if !defined(TFAST)
   if (m_msg->revcomp && (ppst->dnaseq!=1)) {
     fprintf(stderr," cannot reverse complement protein\n");
     m_msg->revcomp = 0;
   }
#endif

   if (ppst->param_u.fa.ktup < 0)
      ppst->param_u.fa.ktup = -ppst->param_u.fa.ktup;

   if (ppst->param_u.fa.ktup < 1 || ppst->param_u.fa.ktup > mktup) {
      fprintf(stderr," warning ktup = %d out of range [1..%d], reset to %d\n",
	      ppst->param_u.fa.ktup, mktup, ppst->param_u.fa.bktup);
      ppst->param_u.fa.ktup = ppst->param_u.fa.bktup;
   }

#ifndef TFASTA
   if (!ppst->sw_flag) strncpy(m_msg->f_id1,"fa",3);
#else
   m_msg->revcomp *= 3;
   if (m_msg->nframe == 3) m_msg->nframe += m_msg->revcomp;
#endif

#if defined(TFASTX) || defined(TFASTY)
   if (m_msg->nframe == 1) m_msg->nframe += m_msg->revcomp;
#endif

#if !defined(TFAST)
  /* for fasta/fastx searches, itt iterates the the query strand */
  m_msg->nitt1 = m_msg->qframe-1;
#else
  /* for tfasta/tfastxy searches, itt iterates the library frames */
  m_msg->nitt1 = m_msg->nframe-1;
#endif

   if (ppst->param_u.fa.ktup>=2 && !wid_set) {
     ppst->param_u.fa.optwid=16;
#if !defined(TFASTA) && !defined(TFASTX) && !defined(TFASTY) && !defined(TFASTF)
#if !defined(FASTX) && !defined(FASTY)
     m_msg->thr_fact = 16;
#else
     m_msg->thr_fact = 8;
#endif
#else
     m_msg->thr_fact = 4;
#endif
   }
   else m_msg->thr_fact = 4;

   if (ppst->param_u.fa.iniflag) {
     ppst->score_ix = 1;
     strncpy (m_msg->label, "initn init1", sizeof(m_msg->label));
   }
   else if (ppst->param_u.fa.optflag) {
     ppst->score_ix = 2;
     m_msg->stages = 1;
   }

   initpam2(ppst);

#if defined(FASTS) || defined(TFASTS)
   if (m_msg->qdnaseq == SEQTYPE_PROT) {
     /* code to make 'L'/'I' identical scores */
     ix_l = pascii['L'];
     ix_i = pascii['I'];
     ppst->pam2[0][ix_l][ix_i] = ppst->pam2[0][ix_i][ix_l] =
       ppst->pam2[0][ix_l][ix_l] = ppst->pam2[0][ix_i][ix_i] =
       (ppst->pam2[0][ix_l][ix_l]+ppst->pam2[0][ix_i][ix_i]+1)/2;
     for (i=1; i<=ppst->nsq; i++) {
       ppst->pam2[0][i][ix_i] = ppst->pam2[0][i][ix_l] =
	 (ppst->pam2[0][i][ix_l]+ppst->pam2[0][i][ix_i]+1)/2;
       ppst->pam2[0][ix_i][i] = ppst->pam2[0][ix_l][i] =
	 (ppst->pam2[0][ix_i][i]+ppst->pam2[0][ix_l][i]+1)/2;
     }

     /* code to make 'Q'/'K' identical scores */
     if (!shift_set) {
       ix_l = pascii['Q'];
       ix_i = pascii['K'];
       ppst->pam2[0][ix_l][ix_i] = ppst->pam2[0][ix_i][ix_l] =
	 ppst->pam2[0][ix_l][ix_l] = ppst->pam2[0][ix_i][ix_i] =
	 (ppst->pam2[0][ix_l][ix_l]+ppst->pam2[0][ix_i][ix_i]+1)/2;
       for (i=1; i<=ppst->nsq; i++) {
	 ppst->pam2[0][i][ix_i] = ppst->pam2[0][i][ix_l] =
	   (ppst->pam2[0][i][ix_l]+ppst->pam2[0][i][ix_i]+1)/2;
	 ppst->pam2[0][ix_i][i] = ppst->pam2[0][ix_l][i] =
	   (ppst->pam2[0][ix_i][i]+ppst->pam2[0][ix_l][i]+1)/2;
       }
     }
   }
#endif

   /* once we have a complete pam matrix, we can calculate Lambda and K 
      for "average" sequences */
   kar_p = NULL;
   init_karlin_a(ppst, aa0_f, &kar_p);
   do_karlin_a(ppst->pam2[0], ppst, aa0_f,
	       kar_p, &m_msg->Lambda, &m_msg->K, &m_msg->H);
   free(kar_p);

#if defined(FASTF) || defined(TFASTF) || defined(FASTS) || defined(TFASTS) || defined(FASTM)
   if (ppst->ext_sq_set) {
     fprintf(stderr," -S not available on [t]fast[fs]\n");
     ppst->ext_sq_set = 0;

     /* reset sascii to ignore -S, map lc */
     init_ascii(0,lascii,0);
   }
#endif
}

void
f_initpam (line, ppst)
char   *line;
struct pstruct *ppst;
{
   if (sscanf (line, " %d %d %d %d %d %d %d", &ppst->param_u.fa.scfact,
	       &ppst->param_u.fa.bestoff, &ppst->param_u.fa.bestscale,
	       &ppst->param_u.fa.bkfact, &ppst->param_u.fa.bktup,
	       &ppst->param_u.fa.bestmax, &ppst->histint) != 7)
   {
      printf ("  bestcut parameters - bad format\n");
      exit (1);
   }
}

/* sortbest has now become comparison function specific so that we can use
   a different comparison for fasts/f 
*/
#if !defined(FASTS) && !defined (TFASTS) && !defined(FASTM)
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

show_aux() {}

#else
void sortbest (bptr, nbest, irelv)
struct beststr **bptr;
int nbest, irelv;
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->escore < bptr[j + gap]->escore) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

#ifndef PCOMPLIB
/* this shuffle is for FASTS */
/* convert ',' -> '\0', shuffle each of the substrings */
qshuffle(unsigned char *aa0, int n0, int nm0)
{
  unsigned char **aa0start, *aap, tmp;
  int i,j,k, ns;

  if ((aa0start=(unsigned char **)calloc(nm0+1,
					 sizeof(unsigned char *)))==NULL) {
    fprintf(stderr,"cannot calloc for qshuffle %d\n",nm0);
    exit(1);
  }
  aa0start[0]=aa0;
  for (k=1,i=0; i<n0; i++) {
    if (aa0[i]==EOSEQ || aa0[i]==ESS) {
      aa0[i]='\0';
      aa0start[k++] = &aa0[i+1];
    }
  }  

  /* aa0start has the beginning of each substring */
  for (k=0; k<nm0; k++) {
    aap=aa0start[k];
    ns = strlen((char *)aap);
    for (i=ns; i>1; i--) {
      j = nrand(i);
      tmp = aap[j];
      aap[j] = aap[i-1];
      aap[i-1] = tmp;
    }
    aap[ns] = 0;
  }

  for (k=1; k<nm0; k++) {
/*  aap = aa0start[k];
    while (*aap) fputc(pst.sq[*aap++],stderr);
    fputc('\n',stderr);
*/
    aa0start[k][-1]=ESS;
  }

  free(aa0start);
}
#endif

/* show additional best_str values */
show_aux(FILE *fp, struct beststr *bptr) {
  fprintf(fp," %d %d",bptr->segnum,bptr->seglen);
}
#endif

void
last_params(unsigned char *aa0, int n0, 
	    struct mngmsg *m_msg,
	    struct pstruct *ppst
#ifdef PCOMPLIB
	    , struct qmng_str *qm_msg
#endif
	    ) {
  int i;

  ppst->n0 = m_msg->n0;

#if defined(FASTF) || defined(TFASTF) || defined(FASTS) || defined(TFASTS)
  m_msg->nm0 = 1;
  for (i=0; i<n0; i++)
    if (aa0[i]==EOSEQ || aa0[i]==ESS) m_msg->nm0++;

  if (m_msg->nm0 > 10) m_msg->escore_flg = 0;
  else m_msg->escore_flg = 1;

  if (m_msg->escore_flg && (ppst->zsflag&1)) {
    m_msg->last_calc_flg = 0;
    m_msg->qshuffle = 0;
  }
  else {	/* need random query, second set of 2000 scores */
    m_msg->last_calc_flg = 1;
    m_msg->qshuffle = 1;
  }
#else
    m_msg->last_calc_flg = 0;
    m_msg->qshuffle = 0;
    m_msg->escore_flg = 0;
    m_msg->nm0 = 1;
#endif

#ifdef PCOMPLIB
  qm_msg->nm0 = m_msg->nm0;
  qm_msg->escore_flg = m_msg->escore_flg;
  qm_msg->qshuffle = m_msg->qshuffle;
#endif
}
