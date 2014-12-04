/* Concurrent read version */

/* $Name: fa34t20b3 $ - $Id: defs.h,v 1.14 2002/08/04 22:46:23 wrp Exp $ */

#ifdef SUNOS
#include <sys/stdtypes.h>
#endif

#if !defined(MAX_WORKERS) && !defined(PCOMPLIB)
#define MAX_WORKERS 1
#endif

#define SEQTYPE_DNA 1
#define SEQTYPE_PROT 0
#define SEQTYPE_UNK -1

#ifndef DEF_NMLEN
#define DEF_NMLEN 6
#endif

/* unfortunately, there is an important relationship between MAXTRN and
   MAXTST+MAXLIB embedded here.  MAXTRN must be >= (MAXTST+MAXLIB)/3
   or it will be possible for a translated DNA sequence to be longer
   than the translation space available */

#define MAX_STR	512 /* standard label/message buffer */
#define MAX_FN  120 /* maximum size of a file name */
#define MAX_CH	40 /* maximum number of library choices */
#ifdef BIGMEM
#define MAX_LF  500 /* maximum numer of library files */
#else
#define MAX_LF  80 /* maximum numer of library files */
#endif

#define MAX_UID 20 /* length of libstr, used for character keys with SQL */

#define AVE_AA_LEN 400
#define AVE_NT_LEN 5000
#define MAX_AA_BUF 5000		/* 5000 later */
#define MAX_NT_BUF 1000		/* 2000 later */

#ifdef BIGMEM
#define MAXTST	20000		/* longest query */
#define MAXLIB	60000		/* longest library */
#define MAXPLIB	600000		/* longest library with p_comp* */
#define MIN_RES 2000		/* minimum amount allocated for alignment */
#ifndef TFASTX
#define MAXTRN  40000		/* buffer for fastx translation */
#else
#define MAXTRN 80000		/* buffer for tfastx translation */
#endif
#define SEQDUP	1200		/* future - overlap */
#ifndef PCOMPLIB
#ifndef MAXBEST
#define MAXBEST	60000	/* max number of best scores */
#endif
#define MAXSTATS 60000
#else
#ifndef MAXBEST
#define MAXBEST	60000	/* max number of best scores */
#endif
#define MAXSTATS 60000
#endif
#define BIGNUM  1000000000
#ifndef MAXINT
#define MAXINT 2147483647
#endif
#define MAXLN	120	/* size of a library name */
#else
#define MAXTST	1500
#define MAXLIB	10000
#define MAXPLIB	100000		/* longest library with p_comp* */
#define MIN_RES 1000
#ifndef TFASTX
#define MAXTRN  4000
#else
#define MAXTRN 11500
#endif
#define SEQDUP	300
#define MAXBEST 2000
#define MAXSTATS 20000
#define BIGNUM  32767
#define MAXINT  32767
#define MAXLN	40	/* size of a library name */
#endif
#if !defined(TFASTA) && !defined(TFASTX)
#define MAXTOT (MAXTST+MAXLIB)
#define MAXDIAG	(MAXTST+MAXLIB)
#else
#define MAXTOT (MAXTST+MAXTRN)
#define MAXDIAG	(MAXTST+MAXTRN)
#endif

#define MAXPAM	600	/* maximum allowable size of the pam matrix */
#define PROF_MAX 500
#define ALF_MAX 30

#ifdef SUPERFAMNUM
#define NSFCHAR '!'
#endif

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
