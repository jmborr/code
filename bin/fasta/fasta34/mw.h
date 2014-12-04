/* Concurrent read version */

/* $Name: fa34t20b3 $ - $Id: mw.h,v 1.10 2002/08/02 00:49:33 wrp Exp $ */

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
typedef long FSEEK_T;
#else
typedef off_t FSEEK_T;
#endif
#endif

struct beststr {
  int n1;		/* sequence length */
  int *n1tot_p;		/* pointer (or NULL) to long sequence length */
  int score[3];		/* score */
  double comp;
  double H;
  double zscore;
  double escore;
  int segnum;
  int seglen;
  struct lmf_str *m_file_p;
  FSEEK_T lseek;
  char libstr[MAX_UID];
  int cont;
  int frame;
  int nsfnum;
  int sfnum[10];
  long loffset;
  struct a_struct *aln_d;	/* these values are used by -m9 */
  float percent, gpercent;
};

struct stat_str {
  int score;
  int n1;
  double comp;
  double H;
  double escore;
  int segnum;
  int seglen;
};


