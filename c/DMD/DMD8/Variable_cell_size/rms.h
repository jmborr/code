#ifndef _RMS__
#define _RMS_

extern double get_rms(double energy);
extern int save_rms(void);
extern int init_rms(void);
extern void close_rms(void);
extern void rms(void);
extern double get_T(void);
extern int is_nucleus(void);
extern int is_rms(void);
extern int n_rms(void);
extern int n0_rms(void);
extern int nonnative(int p, int q, int ct);
#endif
