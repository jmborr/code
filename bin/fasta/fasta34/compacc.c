
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Name: fa34t20b3 $ - $Id: compacc.c,v 1.29 2002/10/10 00:28:53 wrp Exp $ */

/* Concurrent read version */

#include <stdio.h>
#include <stdlib.h>
#if defined(UNIX) || defined(MSDOS)
#include <sys/types.h>
#endif

#include <limits.h>

#include <string.h>
#include <time.h>
#include <math.h>

#include "defs.h"
#include "param.h"
#include "structs.h"

#ifndef PCOMPLIB
#include "mw.h"
#else
#include "p_mw.h"
#endif

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

#ifdef PCOMPLIB
#include "msg.h"
extern int nnodes;
#ifdef PVM_SRC
#include "pvm3.h"
extern int pinums[],hosttid;
#endif
#ifdef MPI_SRC
#include "mpi.h"
#endif
#endif

extern time_t tdone, tstart;		/* Timing */
extern void abort ();
extern void ptime ();


int
reset_maxn(struct mngmsg *m_msg) {
  int maxn;

  maxn = MAXTOT - m_msg->n0 -2;	/* maxn = max library sequence space */
  /* reduce maxn if requested */
  if (m_msg->maxn > 0 && m_msg->maxn < maxn) maxn = m_msg->maxn;

  if (m_msg->qdnaseq==m_msg->ldnaseq || m_msg->qdnaseq==SEQTYPE_DNA) {/* !TFAST */
   if (m_msg->n0>(MAXTST+MAXLIB)/3) {
      fprintf(stderr," query sequence is too long %d %s\n",m_msg->n0,m_msg->sqnam);
      exit(1);
    }
    m_msg->loff = m_msg->n0;
    m_msg->maxt3 = maxn-m_msg->loff;
  }
  else {	/* TFAST */
    if (m_msg->n0 > MAXTST) {
      fprintf(stderr," query sequence is too long %d %s\n",m_msg->n0,m_msg->sqnam);
      exit(1);
    }

    if (m_msg->n0*6 > maxn ) {	/* n0 for the complement + n0*3 for
                                   the three frames */
      if (m_msg->n0*7+2 < MAXTST+MAXLIB) { /* m_msg0*6 + m_msg0 */
	fprintf(stderr,
		" query sequence too long for library segment: %d - resetting to %d\n",
	      maxn,m_msg->n0*6);
	maxn = m_msg->maxn = m_msg->n0*6;
      }
      else {
	fprintf(stderr," query sequence too long for translated search: %d %s\n",
	      m_msg->n0,m_msg->sqnam);
	exit(1);
      }
    }

    /* set up some constants for overlaps */
    m_msg->loff = 3*m_msg->n0;
    m_msg->maxt3 = maxn-m_msg->loff-3;
    m_msg->maxt3 -= m_msg->maxt3%3;
    m_msg->maxt3++;

    maxn = maxn - 3; maxn -= maxn%3; maxn++;
  }
  return maxn;
}


int
scanseq(unsigned char *seq, int n, char *str) {
  int tot,i;
  char aaray[128];		/* this must be set > nsq */
	
  for (i=0; i<128; i++)  aaray[i]=0;
  for (i=0; (size_t)i < strlen(str); i++) aaray[qascii[str[i]]]=1;
  for (i=tot=0; i<n; i++) tot += aaray[seq[i]];
  return tot;
}

void selectbest(bptr,k,n)	/* k is rank in array */
     struct beststr **bptr;
     int k,n;
{
  int v, i, j, l, r;
  struct beststr *tmptr;

  l=0; r=n-1;

  while ( r > l ) {
    v = bptr[r]->score[0];
    i = l-1;
    j = r;
    do {
      while (bptr[++i]->score[0] > v) ;
      while (bptr[--j]->score[0] < v) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

void selectbestz(bptr,k,n)	/* k is rank in array */
     struct beststr **bptr;
     int k,n;
{
  int i, j, l, r;
  struct beststr *tmptr;
  double v;

  l=0; r=n-1;

  while ( r > l ) {
    v = bptr[r]->zscore;
    i = l-1;
    j = r;
    do {
      while (bptr[++i]->zscore > v) ;
      while (bptr[--j]->zscore < v) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

void sortbestz(bptr, nbest)
struct beststr **bptr;
int nbest;
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->zscore >= bptr[j + gap]->zscore) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

void sortbeste(bptr, nbest)
struct beststr **bptr;
int nbest;
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->escore <= bptr[j + gap]->escore) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

extern double zs_to_Ec(double zs, long entries);

/*
extern double ks_dev;
extern int ks_df;
*/
extern char hstring1[];

void
prhist(FILE *fd, struct mngmsg m_msg,
       struct pstruct pst, 
       struct hist_str hist, 
       int nstats,
       struct db_str ntt,
       char *gstring2)
{
  int i,j,hl,hll, el, ell, ev;
  char hline[80], pch, *bp;
  int mh1, mht;
  int maxval, maxvalt, dotsiz, ddotsiz,doinset;
  double cur_e, prev_e, f_int;
  double max_dev, x_tmp;
  double db_tt;
  int n_chi_sq, cum_hl, max_i;


  fprintf(fd,"\n");
  
  if (pst.zsflag < 0 || nstats <= 10) {
    fprintf(fd, "%7ld residues in %5ld sequences\n", ntt.length,ntt.entries);
    fprintf(fd,"\n%s\n",gstring2);
    return;
  }

  max_dev = 0.0;
  mh1 = hist.maxh-1;
  mht = (3*hist.maxh-3)/4 - 1;

  if (!m_msg.nohist && mh1 > 0) {
    for (i=0,maxval=0,maxvalt=0; i<hist.maxh; i++) {
      if (hist.hist_a[i] > maxval) maxval = hist.hist_a[i];
      if (i >= mht &&  hist.hist_a[i]>maxvalt) maxvalt = hist.hist_a[i];
    }
    n_chi_sq = 0;
    cum_hl = -hist.hist_a[0];
    dotsiz = (maxval-1)/60+1;
    ddotsiz = (maxvalt-1)/50+1;
    doinset = (ddotsiz < dotsiz && dotsiz > 2);

    if (pst.zsflag>=0)
      fprintf(fd,"       opt      E()\n");
    else 
      fprintf(fd,"     opt\n");

    prev_e =  zs_to_Ec((double)(hist.min_hist-hist.histint/2),hist.entries);
    for (i=0; i<=mh1; i++) {
      pch = (i==mh1) ? '>' : ' ';
      pch = (i==0) ? '<' : pch;
      hll = hl = hist.hist_a[i];
      if (pst.zsflag>=0) {
	cum_hl += hl;
	f_int = (double)(i*hist.histint+hist.min_hist)+(double)hist.histint/2.0;
	cur_e = zs_to_Ec(f_int,hist.entries);
	ev = el = ell = (int)(cur_e - prev_e + 0.5);
	if (hl > 0  && i > 5 && i < (90-hist.min_hist)/hist.histint) {
	  x_tmp  = fabs(cum_hl - cur_e);
	  if ( x_tmp > max_dev) {
	    max_dev = x_tmp;
	    max_i = i;
	  }
	  n_chi_sq++;
	}
	if ((el=(el+dotsiz-1)/dotsiz) > 60) el = 60;
	if ((ell=(ell+ddotsiz-1)/ddotsiz) > 40) ell = 40;
	fprintf(fd,"%c%3d %5d %5d:",
		pch,(i<mh1)?(i)*hist.histint+hist.min_hist :
		mh1*hist.histint+hist.min_hist,hl,ev);
      }
      else fprintf(fd,"%c%3d %5d :",
		   pch,(i<mh1)?(i)*hist.histint+hist.min_hist :
		   mh1*hist.histint+hist.min_hist,hl);

      if ((hl=(hl+dotsiz-1)/dotsiz) > 60) hl = 60;
      if ((hll=(hll+ddotsiz-1)/ddotsiz) > 40) hll = 40;
      for (j=0; j<hl; j++) hline[j]='='; 
      if (pst.zsflag>=0) {
	if (el <= hl ) {
	  if (el > 0) hline[el-1]='*';
	  hline[hl]='\0';
	}
	else {
	  for (j = hl; j < el; j++) hline[j]=' ';
	  hline[el-1]='*';
	  hline[hl=el]='\0';
	}
      }
      else hline[hl] = 0;
      if (i==1) {
	for (j=hl; j<10; j++) hline[j]=' ';
	sprintf(&hline[10]," one = represents %d library sequences",dotsiz);
      }
      if (doinset && i == mht-2) {
	for (j = hl; j < 10; j++) hline[j]=' ';
	sprintf(&hline[10]," inset = represents %d library sequences",ddotsiz);
      }
      if (i >= mht&& doinset ) {
	for (j = hl; j < 10; j++) hline[j]=' ';
	hline[10]=':';
	for (j = 11; j<11+hll; j++) hline[j]='=';
	hline[11+hll]='\0';
	if (pst.zsflag>=0) {
	  if (ell <= hll) hline[10+ell]='*';
	  else {
	    for (j = 11+hll; j < 10+ell; j++) hline[j]=' ';
	    hline[10+ell] = '*';
	    hline[11+ell] = '\0';
	  }
	}
      }

      fprintf(fd,"%s\n",hline);
      prev_e = cur_e;
    }
  }

  if (ntt.carry==0) {
    fprintf(fd, "%7ld residues in %5ld sequences\n", ntt.length, ntt.entries);
  }
  else {
    db_tt = (double)ntt.carry*(double)LONG_MAX + (double)ntt.length;
    fprintf(fd, "%.0f residues in %5ld library sequences\n", db_tt, ntt.entries);
  }

  if (pst.zsflag>=0) {
    if (MAXSTATS < hist.entries)
#ifdef SAMP_STATS
      fprintf(fd," statistics sampled from %d to %ld sequences\n",
	      MAXSTATS,hist.entries);
#else
      fprintf(fd," statistics extrapolated from %d to %ld sequences\n",
	      MAXSTATS,hist.entries);
#endif
    /*    summ_stats(stat_info); */
    fprintf(fd," %s\n",hist.stat_info);
    if (!m_msg.nohist && cum_hl > 0)
      fprintf(fd," Kolmogorov-Smirnov  statistic: %6.4f (N=%d) at %3d\n",
	      max_dev/(float)cum_hl, n_chi_sq,max_i*hist.histint+hist.min_hist);
    if (m_msg.markx==10) {
      while ((bp=strchr(hist.stat_info,'\n'))!=NULL) *bp=' ';
      if (cum_hl <= 0) cum_hl = -1;
      sprintf(hstring1,"; mp_extrap: %d %ld\n; mp_stats: %s\n; mp_KS: %6.4f (N=%d) at %3d\n",
	      MAXSTATS,hist.entries,hist.stat_info,max_dev/(float)cum_hl, n_chi_sq,max_i*hist.histint+hist.min_hist);
    }
  }
  fprintf(fd,"\n%s\n",gstring2);
  fflush(fd);
}

extern char prog_name[], *verstr;

void s_abort (char *p,  char *p1)
{
  int i;

  fprintf (stderr, "\n***[%s] %s%s***\n", prog_name, p, p1);
#ifdef PCOMPLIB
#ifdef PVM_SRC
  for (i=FIRSTNODE; i< nnodes; i++) pvm_kill(pinums[i]);
  pvm_exit();
#endif
#ifdef MPI_SRC
  MPI_Abort(MPI_COMM_WORLD,1);
  MPI_Finalize();
#endif
#endif
  exit (1);
}

#ifndef MPI_SRC
void w_abort (p, p1)
char *p, *p1;
{
  fprintf (stderr, "\n***[%s] %s%s***\n\n", prog_name, p, p1);
  exit (1);
}
#endif

#ifndef PCOMPLIB
/* copies from from to to shuffling */

extern int nrand(int);

void
shuffle(unsigned char *from, unsigned char *to, int n)
{
  int i,j; unsigned char tmp;

  if (from != to) memcpy((void *)to,(void *)from,n);

  for (i=n; i>1; i--) {
    j = nrand(i);
    tmp = to[j];
    to[j] = to[i-1];
    to[i-1] = tmp;
  }
  to[n] = 0;
}

/* copies from from to from shuffling, ieven changed for threads */
void
wshuffle(unsigned char *from, unsigned char *to, int n, int wsiz, int *ieven)
{
  int i,j, k, mm; 
  unsigned char tmp, *top;

  memcpy((void *)to,(void *)from,n);
	
  mm = n%wsiz;

  if (*ieven) {
    for (k=0; k<(n-wsiz); k += wsiz) {
      top = &to[k];
      for (i=wsiz; i>1; i--) {
	j = nrand(i);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[n-mm];
    for (i=mm; i>1; i--) {
      j = nrand(i);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    *ieven = 0;
  }
  else {
    for (k=n; k>=wsiz; k -= wsiz) {
      top = &to[k-wsiz];
      for (i=wsiz; i>1; i--) {
	j = nrand(i);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[0];
    for (i=mm; i>1; i--) {
      j = nrand(i);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    *ieven = 1;
  }
  to[n] = 0;
}

#endif

int
sfn_cmp(int *q, int *s)
{
  if (*q == *s) return *q;
  while (*q && *s) {
    if (*q == *s) return *q;
    else if (*q < *s) q++;
    else if (*q > *s) s++;
  }
  return 0;
}

#ifndef MPI_SRC
void
revcomp(unsigned char *seq, int n, int *c_nt)
{
  unsigned char tmp;
  int i, ni;

  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = c_nt[seq[i]];
    seq[i] = c_nt[seq[ni]];
    seq[ni] = tmp;
  }
  if ((n%2)==1) {
    i = n/2;
    seq[i] = c_nt[seq[i]];
  }
}
#endif

#ifdef PCOMPLIB

/* init_stage2 sets up the data structures necessary to send a subset
   of sequences to the nodes, and then collects the results
*/

/* wstage2[] FIRSTNODE .. nnodes has the next sequence to be do_opt()/do_walign()ed */
/* wstage2p[] is a list of sequence numbers/frames, to be sent to workers */
/* wstage2b[] is a list of bptr's that shares the index with wstage2p[] */

static int  wstage2[MAXWRKR +1];	/* count of second stage scores */
static struct stage2_str  *wstage2p[MAXWRKR+1]; /* list of second stage sequences */
static int  wstage2i[MAXWRKR+1];	/* index into second stage sequences */
static struct beststr *bbptr,
     **wstage2b[MAXWRKR+1];	/* reverse pointers to bestr */

void
do_stage2(struct beststr **bptr, int nbest, struct mngmsg m_msg0,
	  int s_func, struct qmng_str *qm_msp) {

  int i, is, ib, iw, nres;
  int node, snode, node_done;
  int bufid, numt, tid;
  char errstr[120];
  struct comstr2 bestr2[BFR2+1];	/* temporary structure array */
  char *seqc_buff, *seqc;
  int seqc_buff_len;
#ifdef MPI_SRC
  MPI_Status mpi_status;
#endif

  /* initialize the counter for each worker to 0 */
  for (iw = FIRSTNODE; iw < nnodes; iw++) wstage2[iw] = 0;

  /* for each result, bump the counter for the worker that has
     the sequence */
  for (ib = 0; ib < nbest; ib++ ) { wstage2[bptr[ib]->wrkr]++;  }

  /* now allocate enough space to send each worker a 
     list of its sequences stage2_str {seqnm, frame} */
  for (iw = FIRSTNODE; iw < nnodes; iw++) {
    if (wstage2[iw]>0) {
      if ((wstage2p[iw]=
	   (struct stage2_str *)
	   calloc(wstage2[iw],sizeof(struct stage2_str)))==NULL) {
	sprintf(errstr," cannot allocate sequence listp %d %d",
		iw,wstage2[iw]);
	s_abort(errstr,"");
      }

      /* allocate space to remember the bptr's for each result */
      if ((wstage2b[iw]=(struct beststr **)
	   calloc(wstage2[iw],sizeof(struct beststr *)))==NULL) {
	sprintf(errstr," cannot allocate sequence listb %d %d",
		iw,wstage2[iw]);
	s_abort(errstr,"");
      }
      wstage2i[iw]=0;
    }
    else {
      wstage2p[iw] = NULL;
      wstage2b[iw] = NULL;
    }
  }

  /* for each result, set wstage2p[worker][result_index_in_worker] */
  for (is = 0; is < nbest; is++) {
    iw=bptr[is]->wrkr;
    wstage2p[iw][wstage2i[iw]].seqnm = bptr[is]->seqnm;
    wstage2p[iw][wstage2i[iw]].frame = bptr[is]->frame;
    wstage2b[iw][wstage2i[iw]] = bptr[is];
    wstage2i[iw]++;
  }

  /* at this point, wstage2i[iw] should equal wstage2[iw] */
  node_done = 0;
  for (node = FIRSTNODE; node < nnodes; node++) {

    /* if a worker has no results, move on */
    if (wstage2[node]<=0) { node_done++; continue;}

    qm_msp->slist = wstage2[node];	/* set number of results to return */
    qm_msp->s_func = s_func;		/* set s_funct for do_opt/do_walign */
#ifdef PVM_SRC
    pvm_initsend(PvmDataRaw);
    pvm_pkbyte((char *)qm_msp,sizeof(struct qmng_str),1);
    pvm_send(pinums[node],MSEQTYPE);	/* send qm_msp */
    pvm_initsend(PvmDataRaw);	/* send the list of seqnm/frame */
    pvm_pkbyte((char *)wstage2p[node],wstage2[node]*sizeof(struct stage2_str),1);
    pvm_send(pinums[node],LISTTYPE);
#endif
#ifdef MPI_SRC
    MPI_Send(qm_msp,sizeof(struct qmng_str),MPI_BYTE,node,MSEQTYPE,
	     MPI_COMM_WORLD);
    MPI_Send((char *)wstage2p[node],wstage2[node]*
	     sizeof(struct stage2_str),MPI_BYTE,node,LISTTYPE,
	     MPI_COMM_WORLD);
#endif
  }
	
  /* all the workers have their list of sequences */
  /* reset the index of results to obtain */
  for (iw = 0; iw < nnodes; iw++) wstage2i[iw]=0;
	
  while (node_done < nnodes-FIRSTNODE) {
#ifdef PVM_SRC
    bufid = pvm_recv(-1,LISTRTYPE);	/* wait for results */
    pvm_bufinfo(bufid,NULL,NULL,&tid);
    /* get a chunk of comstr2 results */
    pvm_upkbyte((char *)&bestr2[0],sizeof(struct comstr2)*(BFR2+1),1);
    snode = (iw=tidtonode(tid));
    pvm_freebuf(bufid);
#endif
#ifdef MPI_SRC
    MPI_Recv((char *)&bestr2[0],sizeof(struct comstr2)*(BFR2+1),
	     MPI_BYTE,MPI_ANY_SOURCE,LISTRTYPE,MPI_COMM_WORLD,
	     &mpi_status);
    snode = mpi_status.MPI_SOURCE;
    iw = snode;
#endif

    if (m_msg0.show_code) {
#ifdef PVM_SRC
      bufid = pvm_recv(tid,CODERTYPE);
      pvm_upkint(&seqc_buff_len,1,1);	/* get the code string length */
#endif
#ifdef MPI_SRC
      MPI_Recv((char *)&seqc_buff_len,1,MPI_INT, snode,
	       CODERTYPE,MPI_COMM_WORLD, &mpi_status);
#endif

      if (seqc_buff_len > 0) {		/* allocate space for it */
	if ((seqc=seqc_buff=calloc(seqc_buff_len,sizeof(char)))==NULL) {
	  fprintf(stderr,"Cannot allocate seqc_buff: %d\n",seqc_buff_len);
	  seqc_buff_len=0;
	}
	else {
#ifdef PVM_SRC
	  pvm_upkbyte(seqc_buff,seqc_buff_len*sizeof(char),1);
#endif
#ifdef MPI_SRC
	  MPI_Recv((char *)seqc_buff,seqc_buff_len*sizeof(char),
	     MPI_BYTE,snode,CODERTYPE,MPI_COMM_WORLD, &mpi_status);
#endif
	}
      }
#ifdef PVM_SRC
      pvm_freebuf(bufid);
#endif
    }

    /* get number of results in this message */
    nres = bestr2[BFR2].seqnm & ~FINISHED;
    /* check to see if finished */
    if (bestr2[BFR2].seqnm&FINISHED) {node_done++;}

    /* count through results from a specific worker */
    for (i=0,is=wstage2i[iw]; i < nres; i++,is++) {
      /* get the (saved) bptr for this result */
      bbptr=wstage2b[iw][is];
      /* consistency check seqnm's must agree */
      if (wstage2p[iw][is].seqnm ==  bbptr->seqnm) {
	if (m_msg0.last_calc_flg) {
	  bbptr->escore = bestr2[i].escore;
	  bbptr->segnum = bestr2[i].segnum;
	  bbptr->seglen = bestr2[i].seglen;
	}
	if (m_msg0.stages > 1) {
	  bbptr->score[0] = bestr2[i].score[0];
	  bbptr->score[1] = bestr2[i].score[1];
	  bbptr->score[2] = bestr2[i].score[2];
	}
	if (m_msg0.markx == 9) {
	  /* get score, alignment information, percents */
	  bbptr->swscore = bestr2[i].swscore;
	  memcpy(bbptr->aln_d,&bestr2[i].aln_d,sizeof(struct a_struct));
	  bbptr->percent = bestr2[i].percent;
	  bbptr->gpercent = bestr2[i].gpercent;

	  if (m_msg0.show_code) {	/* if show code */
	    bbptr->aln_code_n = bestr2[i].aln_code_n;	/* length of encoding */
	    if (bbptr->aln_code_n > 0) {
	      if ((bbptr->aln_code = 
		   (char *)calloc(bbptr->aln_code_n+1,sizeof(char)))==NULL) {
		fprintf(stderr,"cannot allocate seq_code: %d\n",
			bbptr->aln_code_n);
		bbptr->aln_code_n = 0;
	      }
	      else {
		strncpy(bbptr->aln_code,seqc,bbptr->aln_code_n);
		bbptr->aln_code[bbptr->aln_code_n]='\0';
		seqc += bbptr->aln_code_n+1;
	      }
	    }
	  }
	}
      }
      else fprintf(stderr,"phase error in phase II return %d %d", iw,i);
    }
    if (seqc_buff != NULL) {
      free(seqc_buff);
      seqc_buff = NULL;
    }
    wstage2i[iw] += nres;
  }

  for (iw=FIRSTNODE; iw < nnodes; iw++) {
    if ((void *)wstage2p[iw]!=NULL) free((void *)wstage2p[iw]);
    if ((void *)wstage2b[iw]!=NULL) free((void *)wstage2b[iw]);
  }
}

#endif
