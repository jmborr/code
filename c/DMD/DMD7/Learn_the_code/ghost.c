#include <math.h>
#include "bcp.h"

#define mult  314159269   //multiplier
#define add   907633385   //aditive constant for rng
#define big  4294967296.0
#define PI  3.141592653589793
#define MIN 1.0e-100

external CollisionData *coll;
unsigned seed;
double T, sigma=sqrt(T) ; //T stands for temperature


/*uniform distribution*/
double u_d(){
  static unsigned rn= seed; 
  rn= rn*mult + add;
  return (double)rn/big;
}/*===========================================matches double u_d()

/*gaussian distribution, not normalized, but lacking 1/sigma */
double g_d(){
  double r= u_d();
  if(r<MIN) return sigma*sqrt(-2*log(MIN)) * cos( PI*2*u_d() );
  return sigma*sqrt(-2*log(r)) * cos( PI*2*u_d() );
}/*============================================matches double g_d()

/*exponentially distributed, positive, random deviate of unit mean */
double expdev(long *idum){  
  double r=u_d();
  if(r<MIN) return -log(MIN);
  return -log(r);
}/*=============================================matches double expdev(...)
   
/*compute new velocities for the collision of real particle with ghost
 particle */
long hit_by_ghost( atom *a, long i1, long ct1 ){
  long i ;
  long k, ct = ct1 ;
  atom *a1 = a+i ;
  moveatom( a1 ) ;
  a1->v.x = g_d() ;
  a1->v.y = g_d() ;
  a1->v.z = g_d() ;
  k = ct ;
  ct = coll[k].next ;
  return ct ;
}/*=============================================matches hit_by_ghost(...)
   
