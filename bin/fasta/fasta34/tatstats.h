#ifndef MAXSQ
#include "param.h"
#endif

#ifndef MAXSAV
#define MAXSAV 30
#endif
#define MAXSEG 30

struct savestr {
  int     score;		/* pam score with segment optimization */
  int     score0;		/* pam score of best single segment */
  int     start0;		/* score from global match */
  int     dp;			/* diagonal of match */
  int     start;		/* start of match in lib seq */
  int     stop;                 /* end of match in lib seq */
#if defined(FASTF) || defined(TFASTF)
  int     *used;                /* array of positions in aa0 that were used */
#endif
};

struct dstruct	{	/* diagonal structure for saving current run */
   int     score;	/* hash score of current match */
   int     start;	/* start of current match */
   int     stop;	/* end of current match */
   struct savestr *dmax;   /* location in vmax[] where best score data saved */
};

struct tat_str {
  double *probs;
  int lowscore;
  int highscore;
};

struct f_struct {
  struct dstruct *diag;
  struct savestr vmax[MAXSAV];	/* best matches saved for one sequence */
  struct savestr *vptr[MAXSAV];
  struct savestr *lowmax;
  int shuff_cnt;
  int nsave;
  int ndo;
  int noff;
  int nm0;
#if defined(FASTS) || defined(TFASTS)
  int *nmoff;		/* offset number, start */
  int *nm_u;
  int *aa0b, *aa0e, *aa0i;
#else
  int nmoff;		/* offset number, start */
  unsigned char *aa0;
  int aa0ix;
#endif
  unsigned char *aa0t;
  int hmask;			/* hash constants */
  int *pamh1;			/* pam based array */
  int *pamh2;			/* pam based kfact array */
#if defined(FASTS) || defined(TFASTS)
  int *link, *harr;		/* hash arrays */
#else
  struct hlstr *link, *harr;            /* hash arrays */
#endif
  int kshft;			/* shift width */
  int nsav, lowscor;		/* number of saved runs, worst saved run */
#if defined(TFASTS) || defined(TFASTF)
  unsigned char *aa1x;
  int n10;
#endif
  struct bdstr *bss;
  struct swstr *ss;
  struct swstr *r_ss;
  int *waa;
  int *res;
  int max_res;
  double *priors;
#if defined(FASTS) || defined(TFASTS)
  struct tat_str **tatprobs;          /* array of pointers to tat structs */
  double **intprobs;                  /* array of integrated tatprobs */
#endif
  int dotat;
  double spacefactor;
};

struct slink {
  int     score;
  double  tatprob;
  struct tat_str *tat;
  struct tat_str *newtat;
  struct savestr *vp;
  struct slink *next;
  struct slink *prev;
};

struct segstr {
  double tatprob;
  int length;
};

void generate_tatprobs(const unsigned char *query,
		       int begin,
		       int end,
		       double *priors,
		       int **pam2,
		       int nsq,
		       struct tat_str **tatarg, struct tat_str *oldtat);

double
calc_tatusov ( struct slink *last,
	       struct slink *this,
	       const unsigned char *aa0, const int n0,
	       unsigned char *aa1, int n1,
	       int **pam2, int nsq,
	       struct f_struct *f_str,
	       int pseudocts,
	       int do_opt,
	       int zsflag
	       );

double seg_tatprob(struct slink *start,
		   unsigned char *aa0,
		   int n0,
		   unsigned char *aa1,
		   int n1,
		   struct f_struct *f_str,
		   struct pstruct *ppst,
		   int do_opt);

void calc_priors(double *priors,
		 struct pstruct *ppst,
		 struct f_struct *f_str,
		 unsigned char *aa1,
		 int n1, int pseudocts);

double factorial (int a, int b);

int max_score(int *scores, int nsq);

int min_score(int *scores, int nsq);

double calc_spacefactor(struct f_struct *f_str);

void linreg(double *lnx, double *x, double *lny,
	    int n,
	    double *a, double *b, double *c, int start);
