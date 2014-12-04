
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa34t20b3 $ - $Id: showrss.c,v 1.5 2001/07/10 18:03:42 wrp Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#ifndef PCOMPLIB
#include "mw.h"
#else
#include "p_mw.h"
#endif

#include "structs.h"
#include "param.h"

extern double
zs_to_E(float zs, int n1, int isdna, long entries,struct db_str db);
extern double zs_to_bit(float zs, int n0, int n1);
extern double zs_to_p(float zs);

extern float (*find_zp)(int score, int length, double comp, void *);

void showbest (FILE *fp, int *s_save, int n1, int qlib, struct mngmsg *m_msg,
	       struct pstruct pst, struct db_str db,
	       void *pstat, char *gstring2)
{
  float zs;
  int score;
  char *rlabel;

  if ((rlabel=strrchr(m_msg->label,' '))==NULL) rlabel = m_msg->label;

  score = s_save[pst.score_ix];
  zs = (*find_zp)(score,n1,1.0,pstat);
  fprintf(fp,"\n PRSS34 - %d shuffles; ",m_msg->mshow);
  if (m_msg->aln.llen > 0)
    fprintf(fp," window shuffle, window size: %d\n",m_msg->aln.llen);
  else
    fprintf(fp," uniform shuffle\n");

  fprintf(fp," unshuffled %s score: %d;  bits(s=%d|n_l=%d): %4.1f p(%d) < %g\n",
	  rlabel,score,score, n1,zs_to_bit(zs,m_msg->n0,n1),score,zs_to_p(zs));

  fprintf(fp,"For %ld sequences, a score >= %d is expected %4.4g times\n\n", 
	  db.entries,score,zs_to_E(zs,n1,0l,db.entries,db)); 
}

void showalign (FILE *fp, unsigned char *aa0, unsigned char *aa1, int maxn,
		struct beststr **bptr, int nbest,int qlib, struct mngmsg m_msg,
		struct pstruct pst, void *f_str, char *gstring2)
{
}

void
aancpy(char *to, char *from, int count,
       struct pstruct pst)
{
  char *tp;
  tp=to;
  while (count-- && *from) {
    if (*from <= pst.nsq) *tp++ = pst.sq[*(from++)];
    else *tp++ = *from++;
  }
  *tp='\0';
}

