/* $Name: fa34t20b3 $ - $Id: param.h,v 1.19 2002/08/02 00:49:34 wrp Exp $ */

#define MAXSQ 50

/* Concurrent read version */

struct fastr {
  int ktup;
  int cgap;
  int pgap;
  int pamfact;
  int scfact;
  int bestoff;
  int bestscale;
  int bkfact;
  int bktup;
  int bestmax;
  int altflag;
  int optflag;
  int iniflag;
  int optcut;
  int optcut_set;
  int optwid;
};

struct prostr {
    int gopen;
    int gextend;
    int width;
};

struct pstruct		/* parameters */
{
  int n0;	/* length of query sequence, used for statistics */
  int gdelval;	/* value for first residue in gap (-12) */
  int ggapval;	/* value for additional residues in gap (-2) */
  int gshift;	/* frameshift for fastx, fasty */
  int gsubs;	/* nt substitution in fasty */
  int p_d_mat;	/* dna match penalty */
  int p_d_mis;	/* dna mismatch penalty */
  int p_d_set;	/* using match/mismatch */
  int score_ix;	/* index to sorted score */
  int zsflag;	/* use scalebest() */
  int zs_win;
  int histint;		/* histogram interval */
  char sq[MAXSQ+1];
  int hsq[MAXSQ+1];
  int nsq;		/* length of normal sq */
  int ext_sq_set;	/* flag for using extended alphabet */
  char sqx[MAXSQ];
  int hsqx[MAXSQ+1];
  int c_nt[MAXSQ+1];
  int nsqx;	/* length of extended sq */
  int dnaseq;	/* -1 = not set (protein); 0 = protein; 1 = DNA; 2 = other */
  int debug_lib;
  int tr_type;	/* codon table */
  int sw_flag;
  char pamfile[40];	/* pam file type */
  int pam_set;
  int **pam2[2];
  int pamoff;	/* offset for pam values */
  int pam_l, pam_h, pam_x;	/* lowest, highest pam value */
  int maxlen;
  long zdb_size; /* force database size */
  union {
    struct fastr fa;
    struct prostr pr;
  } param_u;
  int pseudocts;
  int shuff_node;
};

/* Result structure - do not remove */
struct rstruct
{
  int score[3];
  double comp;
  double H;
  double escore;
  int segnum;
  int seglen;
};

#ifndef PCOMPLIB
struct thr_str {
  int worker;
  void *status;
  int max_work_buf;
  int qframe;
  struct pstruct *pst;
  int qshuffle;
  unsigned char *aa0;
  int n0;
  int nm0;
};

/* this structure passes library sequences to the worker threads
   and returns scores */

struct buf_str {
  int n1;
  int *n1tot_p;
  unsigned char *aa1b;
#ifndef USE_FSEEKO
  long lseek;
#else
  off_t lseek;
#endif
  struct lmf_str *m_file_p;
  int cont;
  int qframe;
  int frame;
  int nsfnum;
  int sfnum[10];
  char libstr[20];	/* set to MAX_UID */
  struct rstruct rst;
  int r_score, qr_score;
  double r_escore, qr_escore;
};

struct buf_head {
  int buf_cnt;
  int have_results;
  unsigned char *start;
  struct buf_str *buf;
};
#endif

/* this definition must be the same as in structs.h */
#ifndef A_STRUCT
#define A_STRUCT
struct a_struct {
  int min0, max0, min1, max1;
  int smin0, smin1, smins;
  int mins;
  int llen, llcntx, llcntx_flg, showall;
  int qlrev, qlfact;
  int llrev, llfact, llmult;
  int frame;
  int a_len;			/* consensus alignment length */
  int nident, ngap_q, ngap_l, nfs;	/* number of identities, gaps in q, l */
  long d_start0,d_stop0;
  long d_start1,d_stop1;
};
#endif
