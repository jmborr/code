/*	doinit.c	general and function-specific initializations */

/* copyright (c) 1996, 1997, 1998  William R. Pearson and the U. of Virginia */

/* $Name: fa34t20b3 $ - $Id: doinit.c,v 1.31 2002/08/28 21:09:17 wrp Exp $ */

/* this file performs general initializations of search parameters

   In addition, it calls several functions in init??.c that provide
   program-specific initializations:

   f_initenv()	- called from initenv()
   f_getopt()	- called from initenv() during a getopt() scan
   f_getarg()	- called from initenv() after the getopt() scan

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __MWERKS__
#define getenv mgetenv
#endif

#include "defs.h"
#include "param.h"
#include "upam.h"

#include "structs.h"

#define XTERNAL
#include "uascii.h"
#undef XTERNAL

extern char *s_optstr;
extern int optind;		/* used by getopt() */

#ifdef PCOMPLIB
#define PARALLEL
#include "p_mw.h"
extern char pgmdir[];
extern char managepgm[];
extern char workerpgm[];
extern int max_buf_cnt;
#define MAX_WORKERS MAXWRKR
#endif

char prog_name[MAX_FN];

extern void f_initenv(struct mngmsg *, struct pstruct *, unsigned char **);
extern void f_lastenv(struct mngmsg *, struct pstruct *);
extern int initpam(char *, struct pstruct *);
extern void f_getopt(char, char *, struct mngmsg *, struct pstruct *);
extern void f_getarg(int, char **, int, struct mngmsg *, struct pstruct *);
extern int standard_pam(char *smstr, struct pstruct *ppst);

int optcnt;
int max_workers=MAX_WORKERS;
#ifdef PCOMPLIB
int worker_1=0;
int worker_n=0;
#endif
extern char *optarg;

/* initenv ()  initializes the environment */
void initenv (int argc, char **argv, struct mngmsg *m_msg, 
		 struct pstruct *ppst, unsigned char **aa0)
{
   char   *cptr, *getenv (), *bp, ctmp;
   int     copt, getopt ();
   char smstr[MAX_FN];		/* smatrix file name */

   /* options for all search functions */
   char   *g_optstr = "ab:BC:d:DHiJ:K:Il:Lm:M:nN:O:pQqr:R:Ss:t:T:v:w:W:x:X:z:Z:";
   char    optstring[MAX_STR];

/*  these initializations will be used by all functions */

   /* prog_name[] is only used for error messages */
   strncpy(prog_name,argv[0],sizeof(prog_name));

#ifdef PARALLEL
   if ((cptr = getenv ("MANAGEPGM")) != NULL) strncpy (managepgm, cptr, 120);
   if ((cptr = getenv ("WORKERPGM")) != NULL) strncpy (workerpgm, cptr, 120);
   if ((cptr = getenv ("PGMDIR")) != NULL) strncpy (pgmdir, cptr, 120);
#endif

   m_msg->ltitle[0] = '\0';

   if ((cptr=getenv("FASTLIBS"))!=NULL)
     strncpy(m_msg->flstr,cptr,MAX_FN);
   else m_msg->flstr[0]='\0';

   m_msg->hist.hist_a = NULL;
   m_msg->outfile[0] = '\0';
   m_msg->ldnaseq = 0;	/* library is protein */
   m_msg->n1_low = 0;
   m_msg->n1_high = BIGNUM;
   m_msg->ql_off = 1;	/* start with first query sequence */

   m_msg->pamd1 = MAXSQ;
   m_msg->pamd2 = MAXSQ;

   ppst->tr_type = 0;
   ppst->debug_lib = 0;
   m_msg->nshow = 20;
#ifdef PCOMPLIB
   m_msg->nohist = 1;
   m_msg->mshow = 20;
#else
   m_msg->nohist = 0;
   m_msg->mshow = 50;
#endif
   m_msg->ashow = -1;
   m_msg->nmlen = DEF_NMLEN;
   m_msg->z_bits = 1;
   m_msg->mshow_flg = 0;
   m_msg->aln.llen = 60;
   m_msg->aln.llcntx = 30;
   m_msg->aln.llcntx_flg = 0;
   m_msg->e_cut = 10.0;
   m_msg->e_low = 0.0;
   m_msg->e_cut_set = 0;
   m_msg->revcomp = 0;
   m_msg->self = 0;
   m_msg->long_info = 0;
   m_msg->maxn = 0;
   m_msg->dupn = SEQDUP;
   m_msg->dfile[0] = '\0';
   m_msg->tname[0] = '\0';
   m_msg->lname[0] = '\0';
   m_msg->show_code = 0;
   m_msg->aln.showall = 0;
   m_msg->markx = 0;
   m_msg->sq0off = m_msg->sq1off = 1;
   strncpy(m_msg->sqnam,"aa",4);
   strncpy(m_msg->sqtype,"protein",10);
   
   ppst->zsflag = 1;
   ppst->zs_win = 0;
   ppst->pam_x = 1;  /* set >0 to use pam['X']['X'] value */
   ppst->pam_set = 0;
   ppst->p_d_set = 0;
   ppst->zdb_size = -1;
   ppst->dnaseq = 0;	/* default is protein */
   ppst->pamoff = 0;
   ppst->ext_sq_set = 0;

   if ((cptr = getenv ("SMATRIX")) != NULL)
   {
      strncpy (smstr, cptr, MAX_FN);
      strncpy (ppst->pamfile, smstr, MAX_FN);
   }

   f_initenv (m_msg, ppst, aa0);

   strncpy (optstring, g_optstr, sizeof (optstring));
   strncat (optstring, s_optstr, sizeof (optstring));

   while ((copt = getopt (argc, argv, optstring)) != EOF)
   {
      if (strchr (g_optstr, copt) != NULL)
      {
	switch (copt) {  /* switches for all options */
	case 'a': m_msg->aln.showall = 1; break;
	case 'B': m_msg->z_bits = 0; break;
	case 'b':
	  if (optarg[0] == '$') {
	    m_msg->mshow = -1;
	    m_msg->e_cut = 10000000.0;
	    break;
	  }
	  else sscanf (optarg, "%d", &m_msg->mshow);
	  m_msg->e_cut = 10000000.0;
	  m_msg->e_cut_set = 1;
	  m_msg->mshow_flg = 1;
	  break;
	case 'C': sscanf(optarg,"%d",&m_msg->nmlen);
	  if (m_msg->nmlen > MAX_UID-1) m_msg->nmlen = MAX_UID-1;
	  break;
	case 'd': sscanf(optarg,"%d",&m_msg->ashow);
	  if (m_msg->ashow > m_msg->mshow) m_msg->mshow=m_msg->ashow;
	  /* m_msg->ashow_flg = 1; (ashow_flg not in structs.h, not used)*/
	  break;
	case 'D': ppst->debug_lib = 1;
	  break;
	case 'H':
#ifndef PCOMPLIB
	  m_msg->nohist = 1; break;
#else
	  m_msg->nohist = 0; break;
#endif
	case 'i':
	  m_msg->revcomp = 1; break;
#ifdef PARALLEL
	case 'I':
	  m_msg->self = 1; break;
	case 'J':
	  sscanf(optarg,"%d",&m_msg->ql_off);
	  break;
	case 'K':
	  sscanf(optarg,"%d",&max_buf_cnt);
	  break;
#endif
	case 'l':
	  strncpy(m_msg->flstr,optarg,MAX_FN);
	  break;
	case 'L':
	  m_msg->long_info = 1; break;
	case 'M':
	  sscanf(optarg,"%d-%d",&m_msg->n1_low,&m_msg->n1_high);
	  if (m_msg->n1_low < 0) {
	    m_msg->n1_high = -m_msg->n1_low;
	    m_msg->n1_low = 0;
	  }
	  if (m_msg->n1_high == 0) m_msg->n1_high = BIGNUM;
	  if (m_msg->n1_low > m_msg->n1_high) {
	    fprintf(stderr," low cutoff %d greater than high %d\n",
		    m_msg->n1_low, m_msg->n1_high);
	    m_msg->n1_low = 0;
	    m_msg->n1_high = BIGNUM;
	  }
	  break;
	case 'm':
	  sscanf(optarg,"%d",&m_msg->markx);
	  sscanf(optarg,"%d%c",&m_msg->markx,&ctmp);
	  if (m_msg->markx==9 && ctmp=='c') {
	    m_msg->show_code = 1;
	  }
	  if (m_msg->markx > 6 && m_msg->markx != 10 && m_msg->markx != 9)
	    m_msg->markx = 0;
	  break;
	case 'n':
	  m_msg->qdnaseq = 1;
	  memcpy(qascii,nascii,sizeof(qascii));
	  strncpy(m_msg->sqnam,"nt",4);
	  break;
	case 'N':
	  sscanf(optarg,"%d",&m_msg->maxn);
	  break;
	case 'p':
	  m_msg->qdnaseq = 0;
	  ppst->dnaseq = 0;
	  strncpy(m_msg->sqnam,"aa",4);
	  break;
	case 'O':
	  strncpy(m_msg->outfile,optarg,MAX_FN);
	  break;
	case 'q':
	case 'Q':
	  m_msg->quiet = 1;
	  break;
	case 'r':
	  sscanf(optarg,"%d/%d",&ppst->p_d_mat,&ppst->p_d_mis);
	  if (ppst->p_d_mat > 0 && ppst->p_d_mis < 0) {
	    ppst->p_d_set = 1;
	    strncpy(ppst->pamfile,optarg,40);
	  }
	  break;
	case 'R':
	  strncpy (m_msg->dfile, optarg, MAX_FN);
	  break;
	case 's':
	  strncpy (smstr, optarg, MAX_FN);
	  smstr[MAX_FN-1]='\0';
	  if ((bp=strchr(smstr,'-'))!=NULL) {
	    ppst->pamoff=atoi(bp+1);
	    *bp = '\0';
	  }
	  else if ((bp=strchr(smstr,'+'))!=NULL) {
	    ppst->pamoff= -atoi(bp+1);
	    *bp = '\0';
	  }
	  if (!standard_pam(smstr,ppst)) {
	    initpam (smstr, ppst);
	  }
	  ppst->pam_set=1;
	  break;
	case 'S':	/* turn on extended alphabet for seg */
	  ppst->ext_sq_set = 1;
	  break;
	case 't':
	  sscanf (optarg, "%d", &ppst->tr_type);
	  break;
	case 'T':
#ifdef PCOMPLIB
	  if (strchr(optarg,'-') != NULL) {
	    sscanf(optarg,"%d-%d",&worker_1,&worker_n);
	    if (worker_1 > worker_n) {
	      worker_1 = worker_n = 0;
	    }
	  }
	  else 
#endif
	    sscanf (optarg, "%d", &max_workers);
	  if (max_workers < 0) max_workers=1;
	  break;
	case 'v':
	  sscanf (optarg,"%d",&ppst->zs_win);
	  break;
	case 'w':
	  sscanf (optarg,"%d",&m_msg->aln.llen);
	  if (m_msg->aln.llen < 10) m_msg->aln.llen = 10;
	  if (m_msg->aln.llen > 200) m_msg->aln.llen = 200;
	  if (!m_msg->aln.llcntx_flg) m_msg->aln.llcntx = m_msg->aln.llen/2;
	  break;
	case 'W':
	  sscanf (optarg,"%d",&m_msg->aln.llcntx);
	  m_msg->aln.llcntx_flg = 1;
	  break;
	case 'X':
	  sscanf (optarg,"%ld %ld",&m_msg->sq0off,&m_msg->sq1off); break;
	case 'x':
	  sscanf (optarg,"%d",&ppst->pam_x);
	  break;
	case 'z':
	  sscanf(optarg,"%d",&ppst->zsflag);
	  break;
	case 'Z':
	  sscanf(optarg,"%ld",&ppst->zdb_size);
	  break;
	}
      }
      else if (strchr (s_optstr, copt))
	 f_getopt (copt, optarg, m_msg, ppst);
   }
   optind--;

   f_lastenv (m_msg, ppst);

   if (argc - optind < 3) return;
   if (argc - optind > 1) strcpy (m_msg->tname, argv[optind + 1]);
   if (argc - optind > 2) { strcpy(m_msg->lname, argv[optind + 2]); }
   m_msg->tnamesize = sizeof (m_msg->tname);
   f_getarg (argc, argv, optind, m_msg, ppst);
}
