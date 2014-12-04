/* Concurrent read version */

/* $Name: fa34t20b3 $ - $Id: p_mw.h,v 1.12 2002/08/28 21:09:18 wrp Exp $ */

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
typedef long FSEEK_T;
#else
typedef off_t FSEEK_T;
#endif
#endif

struct beststr {
  int n1;		/* sequence number */
  int score[3];		/* score */
  int rscore;	/* score from shuffled sequence */
  int swscore;	/* optimal score from alignment */
  double comp;	/* karlin 1/lambda comp.parameter */
  double H;	/* karlin H information content */
  double zscore;
  double escore;
  double r_escore;
  int segnum;
  int seglen;
  int  lib;
  FSEEK_T lseek;
  int cont;
  int frame;
  int m_seqnm;
  int seqnm;
  int wrkr;
  struct sql *desptr;
  struct a_struct *aln_d;
  char *aln_code;
  char aln_code_n;
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

/* this structure passes library sequences to the worker threads
   and returns scores */

#include "w_mw.h"

/*
struct pbuf_head {
  int buf_cnt;
  unsigned char *start;
  struct sqs2 *buf;
};
*/
