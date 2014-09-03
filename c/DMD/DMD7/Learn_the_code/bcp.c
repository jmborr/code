#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "bcp.h"
#include "controls.h"
#include "movie.h"
#include "make_system.h"
#include "cluster.h"
#include "search.h"
#include "corr_func.h" 
#include "rms.h"
#include "bonds.h"
#include "fs_output.h"
#include "schedule.h"
/*#include<ghost.h>*/

char* text_name="junk_bcp.txt";
int text_no=10000;
FILE  * echo_path , * movie_path ;
FILE * text_path /*input configuration*/;
int done,is_open,m_is_open,t_is_open;
int mticks=20;
int p1,q1,ct1;
int n_p_mes;
double maxtime;

int xyz[4];
well_type ntot;/* wells indeces larger or equal to ntot, correspond to gaps 
inroduced to prevent overlap of hardspheres and breaking permanent bond
well indeces <nen corespond to wells with finite energies;
well indeces >=nen1 corespond to wells with infinite energies;*/

double gap;
double L_limit[3];
double press_limit;
double press_coeff;
double press_curr;
static double max_int_dist;
int n_gap_mes;
int max_gap_mes;
int new_density=0;
int var_density=0;
static int num_gaps;

int n,n1,n2,n3,nen,nen1,ll,llp,coll_type=0,deltall;
int nrt;
double timec;/*elapsed time from the beginning of the simulation. This is 
	       useful when the initial time is not zero*/
double delta,delta1,delta2,delta3,delta4,delta5,delta6,dticks,deltap,timed;
double timep;/*time elapsed since the last rescaling of "temperature"*/
double vvm,corr,corr_2,virial;
double ts;/*number of computed collisions*/
double timea;/*scheduled time for the collision about to be carried out*/
double timeb;/*it can be the time when last computed collision happened or the same as timea, depending on the context*/
double pressure,temperature,volume,dim;
double vvmtime;/*time integral of the kinetic energy*/
double mes_time,avpres,avtemp;
double temp_limit;/*heat bath temperature when using the thermal coefficient for heating rate*/
double coeff;/*thermal coefficient of heating rate. Initialized to zero*/
static double dblarg1=DBL1,dblarg2=DBL2; 

*double tballtime=0;
double totaltime=0;*/
double potential;
double avePot, avpot; 
double potTime;/*time integral of the potential energy*/
atom *a;
well_type ** ecoll;
well_type ** icoll;

well_type * collp; /*gloval variabele (scope is bcp.c file) to store coll.  */
well_type * collq; /*types of atoms with previous scheduled collisions with */
                   /*"p" and "q".Initialized in "make_key_system"           */
CollisionData *coll;
ReactionData *react;
dimensions bound[3];
crd o={0.0,0.0,0.0};
static int px;    /*maximum cell number along the x-axis. We have*/
static int py;    /*px = bound[0].period - 1 */
static int pz;
static double lx; /*system size along the x-axis. W have*/
static double ly; /*lx=bound[0].length*/
static double lz;

/* inpt - the pointer to the record with the next collision time.
*p1,q1,timea are the numbers of the particles next to collide and their
collision time, they are determined by the squeezetable which kill all
the records that contain the particles from the table. If q1>=n1, number of
atoms, (n is the number of atoms -1) 
it denotes the collision with the wall: q1-n is the direction
of the wall: x=1,y=2,z=3. twall defines with which wall it collides first,
newloc defines to which box the particle will go by means of reading
a triad r,v,i from a correct place of the atom structure. i is the
box number in each direction which is integer, but define as union with
double in order to maintaine the same length as velocity and position in
order to be able to read the triad by shifting. */

/*next_double returns the input double x plus a little bit. Let
X be the real number that x wants to represent. Then X=x+d, 
where d is the cut-off due to rounding errors in double. Thus,
next_double adds to x one digit in the last significative 
figure of the double precission*/
double next_double(double x)
{
  int y;
  double z=x;  
  z=frexp(z,&y)+ldexp(DBL_EPSILON,-1);
  z=ldexp(z,y); 
  return z;
}

double prev_double(double x)
{
  int y;
  double z=x;
  z=frexp(z,&y)-ldexp(DBL_EPSILON,-1);
  z=ldexp(z,y); 
  return z;
}


void set_new_bounds(double *L, double maxrb,int ndim)
{     
  double ccell=maxrb;
  set_search_volume(maxrb);
  int i,y;
  max_int_dist=maxrb;
  ccell=next_double(ccell);/*add one digit to the
last significative figure in the double precission*/
  dim=ndim;
  for(i=0;i<ndim;i++)
    {
      bound[i].period=L[i]/ccell;/*take only the integer part of this division to define bound[i].period*/
      
      if(bound[i].period<3)/*we need at least three cells
			     along each dimension*/
	{
	  bound[i].period=3;
	  bound[i].dl=ccell;
	  bound[i].length=ccell*3;
	}
      else
	{
	  bound[i].dl=L[i]/bound[i].period;/* in general, bound[i].dl will be bigger than maxrb, because we only took the integer part of L[i]/ccell in order to define bound[i].period. Thus, the "real" cell size is not maxrb, but slightly bigger*/
	  bound[i].length=bound[i].dl*bound[i].period;
          if((bound[i].dl<=maxrb)||(bound[i].length<L[i]))
	    {
	      bound[i].dl=next_double(bound[i].dl);
	    }
	  bound[i].length=bound[i].dl*bound[i].period;
	}
      printf("%d %lf %d %lf\n",i,bound[i].length,bound[i].period,bound[i].dl);
    }
  if (ndim==2)
    bound[2].length=bound[2].dl=bound[2].period=1;/*no need for
three cells along the third dimensions, since no calculations
will be done*/
}

void init_parameters(void)
{      
  corr_2=0;
  corr=0; /*square root of corr_2*/
  timea=0; 
  timeb=0;
  timec=0;
  timed=0;
  ts=0;
  ll=0;
  maxtime=0;
  m_is_open=0;
  is_open=0;
  t_is_open=0;
  return;
}
 
void init_update_param(int * is_x)
{  
  int i;
/*momemtum of the center of mass and mass of the whole system*/
  double mvx=0,mvy=0,mvz=0,mass=0;
  px=bound[0].period-1; /*maximum cell number along x-axis,*/
  py=bound[1].period-1; /*cell numbers begin with zero*/
  pz=bound[2].period-1;
  lx=bound[0].length; /*syze of system along x-axis*/
  ly=bound[1].length;
  lz=bound[2].length;

      
  n2=n+2; /*n1 in wall index along x-axis, n2 is wall index*/
  n3=n+3; /*along y-axis, and n3 for z-axis*/
  vvm=0.0; /* vvm will store kinetic energy */
   for(i=0;i<=n;i++)
    {

      if(a[i].r.x>=lx) a[i].r.x-=lx;/*deal the periodic*/
      if(a[i].r.x<0) a[i].r.x+=lx;  /*boundary conditions*/
      if(a[i].r.y>=ly) a[i].r.y-=ly;
      if(a[i].r.y<0) a[i].r.y+=ly;
      if(a[i].r.z>=lz) a[i].r.z-=lz;
      if(a[i].r.z<0) a[i].r.z+=lz;

      a[i].i.x.i = a[i].r.x / bound[0].dl; /*initialize the cell*/
      a[i].i.y.i=a[i].r.y/bound[1].dl; /*indexes for every atom*/
      a[i].i.z.i=a[i].r.z/bound[2].dl; 
      /*add is a single digit address for the three cell 
	index numbers*/
      a[i].add=a[i].i.x.i+(a[i].i.y.i+a[i].i.z.i*bound[1].period)*bound[0].period;
      a[i].t=timeb;/*initialized to zero*/
      mvx+=a[i].m*a[i].v.x;
      mvy+=a[i].m*a[i].v.y;
      mvz+=a[i].m*a[i].v.z;
      mass+=a[i].m;
    }
  mvx/=mass;
  mvy/=mass;
  mvz/=mass;
  for(i=0;i<=n;i++)
    {
      a[i].v.x-=mvx;/*set total momemtum to zero*/
      a[i].v.y-=mvy;
      a[i].v.z-=mvz;
      vvm+=a[i].m*dist(a[i].v,o);/* o={0.0,0.0,0.0} */
    }
  /*printf("%lf %lf %lf \n",mvx,mvy,mvz);*/
  vvm*=0.5; /*kinetic energy*/
  if(!corr_2){/*corr_2==0 at the beginning of the simulation.*/
    corr_2=1; 
    corr=sqrt(corr_2);
  }
  delta2=0;
  delta4=0;
  delta6=0;
  delta1=1000;
  delta3=1000;
  delta5=1000;
  dticks=60;
  xyz[0]=1;
  xyz[2]=2;
  xyz[3]=3;
  xyz[1]=1;
  deltall=n1;/*number of time steps between two updates of the temperature of the system by means of the thermal coefficient. By default set to number of particles*/
  if(deltall<DELTALL)deltall=DELTALL;
  potential=0;
  num_gaps=0;
  new_density=0;
  var_density=0;
  pressure=dblarg1;
  temperature=dblarg1;
  avePot=dblarg1;
  mes_time=dblarg1;
  coeff=0;
  temp_limit=0;
  timep=0;
  potTime=0; 
  vvmtime=0;
  virial=0;
  volume=1;
  dim=0;
  for(i=0;i<3;i++)
    if(is_x[i])
      {
	L_limit[i]=bound[i].length;
	volume*=bound[i].length;
	dim++;
      }
  if(!is_x[2])bound[2].period=1;/*bidimensional system*/
  return; 
 } 

int cleanup(void){  
  int i;
  double mvx=0,mvy=0,mvz=0,mass=0;/*center of mass, and totall mass*/
  double vvm1;
  double corr1=1/corr;
  px=bound[0].period-1;
  py=bound[1].period-1;
  pz=bound[2].period-1;
  lx=bound[0].length;
  ly=bound[1].length;
  lz=bound[2].length;
  /*update atoms positions for time===timeb, but do not check per.bound.cond*/
  update_atoms();
  /*store positions of atoms for time=timeb in "crd q" */
  moveatoms();
  /*scale velocities to match the temperature. This means collision
    times should be recomputed!*/
  for( i=0;i<n1;i++)
    {
      a[i].v.x*=corr1;
      a[i].v.y*=corr1;
      a[i].v.z*=corr1;
    }
  reset_colldata();/*set potential walls to "real" values*/
  vvm/=corr_2;  /*scale the kinetic energy*/
  corr_2=corr=1; /*because system T == real T now*/
  timed+=timec; /*store the elapsed time*/
  timec=timeb=0;/*initialize system time to zero*/
  for(i=0;i<=n;i++)
    {
      if(a[i].r.x>=lx) a[i].r.x-=lx;/*implement periodic boundary cond.*/
      if(a[i].r.x<0) a[i].r.x+=lx;
      if(a[i].r.y>=ly) a[i].r.y-=ly;
      if(a[i].r.y<0) a[i].r.y+=ly;
      if(a[i].r.z>=lz) a[i].r.z-=lz;
      if(a[i].r.z<0) a[i].r.z+=lz;
      a[i].i.x.i=a[i].r.x/bound[0].dl;
      a[i].i.y.i=a[i].r.y/bound[1].dl;
      a[i].i.z.i=a[i].r.z/bound[2].dl;            
      a[i].add=a[i].i.x.i+(a[i].i.y.i+a[i].i.z.i*bound[1].period)*bound[0].period;
      a[i].t=timeb;/*which by the way is zero*/
      mvx+=a[i].m*a[i].v.x;/*calculating center of mass*/
      mvy+=a[i].m*a[i].v.y;
      mvz+=a[i].m*a[i].v.z;
      mass+=a[i].m;
    }
  mvx/=mass;
  mvy/=mass;
  mvz/=mass;
  vvm1=0.0;/*Initialize the kinetic energy to zero for future update*/
  /*Calculate velocities so that momentum and angular momentum of the
    whole system are both zero, store velocities in "cdr u" of each
    atom. The resulting kinetic energy will be smaller*/
  stop_atoms((moved_iatom*)a,n1); 
  for(i=0;i<=n;i++)/*assign the velocities from "cdr u" to "cdr v"*/
    {
      a[i].v.x=a[i].u.x;
      a[i].v.y=a[i].u.y;
      a[i].v.z=a[i].u.z;
      vvm1+=a[i].m*dist(a[i].v,o);/*add to the kinetic energy*/
    }
  /*
  for(i=0;i<=n;i++)
    {
      a[i].v.x-=mvx;
      a[i].v.y-=mvy;
      a[i].v.z-=mvz;
      vvm1+=a[i].m*dist(a[i].v,o);
    }
  */
  printf("%lf %lf %lf \n",mvx,mvy,mvz);
  vvm1*=0.5;/*final kinetic energy*/
  printf("%lf %lf\n",vvm1,vvm);
  if(!vvm1)return 0;
  vvm1=sqrt(vvm/vvm1);
  /*rescale the velocities so that the potential energy will be again
    vvm. We are adding vvm-vvm1 energy to the system, but in such a
    way that both the momentum and angular momentum are zero*/
  for(i=0;i<=n;i++)
    {
      a[i].v.x*=vvm1;
      a[i].v.y*=vvm1;
      a[i].v.z*=vvm1;
    }
  potential=0;
  num_gaps=0;
  pressure=dblarg1;
  temperature=dblarg1;
  avePot=dblarg1;
  mes_time=dblarg1;
  potTime=0;
  timep=0;
  vvmtime=0;
  virial=0;
  i=init_tables();/*reset all collision times in the system*/
  if(i!=1)return 0;
  return 1; 
 } 


int change_density(void)
{  
  int i,j;
  double mvx=0,mvy=0,mvz=0,mass=0;
  double vvm1;
  double factor[3];
  double L[3]; 
  int ndim=dim;


  double corr1=1/corr;
  iatom *b =(iatom *)a;

  if(var_density)
    {
      for(i=0;i<ndim;i++)
	if (L_limit[i]>bound[i].length)
	  { 
	    if(L_limit[i]<bound[i].length*(1+gap))
	      L[i]=L_limit[i];
	    else
	      L[i]=bound[i].length*(1+gap*0.99);
	  }
	else	    
	  { 
	    if(L_limit[i]>bound[i].length/(1+gap))
	      L[i]=L_limit[i];
	    else
	      L[i]=bound[i].length/(1+gap*0.99);
	  }
      
      var_density=0;
      for(i=0;i<ndim;i++)
	if(L[i]!=L_limit[i])
	  var_density=1;
    }
  else
  {
    double dv=1+(press_curr/n_gap_mes-press_limit)*press_coeff;
    if(dv>=1+gap)dv=1+gap*0.99;
    if(dv*(1+gap)<1)dv=1/(1+gap*0.99);
    for(i=0;i<ndim;i++)
      L[i]=bound[i].length*dv;
  }
  for(i=0;i<ndim;i++)
    factor[i]=L[i]/bound[i].length;
  n_gap_mes=0;
  press_curr=0;

  update_atoms();
  moveatoms();

  for(j=0;j<3;j++)
    for(i=0;i<=n;i++)
    {
      b[i].r[j]*=factor[j];
    }

  set_new_bounds(L,max_int_dist,(int)dim);
  px=bound[0].period-1;
  py=bound[1].period-1;
  pz=bound[2].period-1;

  volume=1;
  for(i=0;i<ndim;i++)
    volume*=bound[i].length;

  lx=bound[0].length;
  ly=bound[1].length;
  lz=bound[2].length;


  for( i=0;i<n1;i++)
    {
      a[i].v.x*=corr1;
      a[i].v.y*=corr1;
      a[i].v.z*=corr1;
    }
  reset_colldata();
  vvm/=corr_2;    
  corr_2=corr=1;
  timed+=timec;
  timec=timeb=0;
   for(i=0;i<=n;i++)
    {
      if(a[i].r.x>=lx) a[i].r.x-=lx;
      if(a[i].r.x<0) a[i].r.x+=lx;
      if(a[i].r.y>=ly) a[i].r.y-=ly;
      if(a[i].r.y<0) a[i].r.y+=ly;
      if(a[i].r.z>=lz) a[i].r.z-=lz;
      if(a[i].r.z<0) a[i].r.z+=lz;
      a[i].i.x.i=a[i].r.x/bound[0].dl;
      a[i].i.y.i=a[i].r.y/bound[1].dl;
      a[i].i.z.i=a[i].r.z/bound[2].dl;            
      a[i].add=a[i].i.x.i+(a[i].i.y.i+a[i].i.z.i*bound[1].period)*bound[0].period;
      a[i].t=timeb;
      mvx+=a[i].m*a[i].v.x;
      mvy+=a[i].m*a[i].v.y;
      mvz+=a[i].m*a[i].v.z;
      mass+=a[i].m;
    }
  mvx/=mass;
  mvy/=mass;
  mvz/=mass;
  vvm1=0.0;

  stop_atoms((moved_iatom*)a,n1); 
  for(i=0;i<=n;i++)
    {
      a[i].v.x=a[i].u.x;
      a[i].v.y=a[i].u.y;
      a[i].v.z=a[i].u.z;
      vvm1+=a[i].m*dist(a[i].v,o);
    }
  printf("%lf %lf %lf \n",mvx,mvy,mvz);
  vvm1*=0.5;
  printf("%lf %lf\n",vvm1,vvm);
  if(!vvm1)return 0;
  vvm1=sqrt(vvm/vvm1);
  for(i=0;i<=n;i++)
    {
      a[i].v.x*=vvm1;
      a[i].v.y*=vvm1;
      a[i].v.z*=vvm1;
    }
  realloc_search();
  potential=0;
  num_gaps=0;
  pressure=dblarg1;
  temperature=dblarg1;
  avePot=dblarg1;
  mes_time=dblarg1;
  potTime=0;
  timep=0;
  vvmtime=0;
  virial=0;
  new_density=0;
  i=init_tables();
  if(i!=1)return 0;
  return 1; 
 } 


void add_potential(int ct)
{
  int k=ct;
  if((k>=ntot)||(coll[k].prev>=ntot))
    num_gaps++;             /*"etot" is the the difence in potential energy  */
  else                      /*from an infinite distance (zero potential) to  */
    potential+=coll[k].etot;/*the present distance                           */
}

int get_delta_ll(void)
{
  return deltall;
}

void set_delta_ll(int new_deltall)
{
  deltall=new_deltall;
}

double get_mes_time(void)
{
  if (mes_time!=dblarg1)
    return mes_time+timed;
  else 
    return timec+timed;
}

double get_temperature(void)
{
  /*at the beginning of the simulation , temperature=dblarg1*/
  if(temperature != dblarg1)
    return temperature;
  else if(timep)
    /*we average the kinetic energy in a "timep" lapse of time. This
  lapse of time corresponds to the time it takes "deltall" collisions
  to happen*/
    return 2*(vvmtime/timep)/(n1*dim*corr_2);
  else 
    /*if there is no average of the kinetic energy in time */
    return 2*vvm/(n1*dim*corr_2);
}

double get_avePot(void)
{
  if (avePot!=dblarg1)
  return avePot;
  else if(timep)
    return potTime/timep;
  else 
    return potential;
}

void set_temp(double temp)/*temp is the temperature at the beginning of the simulation, or every time we re-schedule*/
{ 
  int i;
  double vvmo;
  vvmo=(temp*n1*dim)/2;/*kinetic energy at equilibrium for T=temp, like 3NkT/2*/
  corr_2=vvm/vvmo;/*ratio of the old kinetic energy to the new kinetic energy*/
  for (i=0;i<nen;i++)
    {
      coll[i].e=coll[i].eo*corr_2;/*modify the potential enery steps inversely to the change of temperature*/
      coll[i].edm=coll[i].edmo*corr_2;
    }
  corr=sqrt(corr_2);
  llp=0;
  virial=0;
  vvmtime=0;
  timep=0;
}

double get_temp(void)
{return 2*vvm/(n1*dim*corr_2);} /*vvm is the total kinetic energy*/


void rescale(void)/*is called regularly, every "deltall" collisions*/
{ int i;
  double temp0=get_temperature();
  if (coeff)
    {
      double coeff1=coeff*(timep*corr);/*remember that time is also changed*/
      if(coeff1>1)coeff1=1;
      corr_2*=temp0/(temp0*(1.0-coeff1)+temp_limit*coeff1);  
      for (i=0;i<nen;i++)
	{
	  coll[i].e=coll[i].eo*corr_2;
	  coll[i].edm=coll[i].edmo*corr_2;
	}
      corr=sqrt(corr_2);
    }
  llp=0;
  virial=0;
  vvmtime=0;
  potTime=0;
  timep=0;
}

void set_temp_limit(double t)
{
  temp_limit=t;
}

double get_temp_limit(void)
{
  return temp_limit;
}

void set_coeff(double c)
{
  coeff=c;
}

double get_coeff(void)
{
  return coeff;
}

double get_rate(void)
{return delta5;}

void set_rate(double rate)
{delta5=rate;delta6=0;}
/*
double get_rate(void)
{return 3600/(double)dticks;}

void set_rate(double rate)
{ 
 if (rate<3600/2000000000.0) dticks=2000000000;
 else dticks=3600/rate;
}
*/
int open_echo_file(int is_open,  char * fname)
{
  int nbyte,i,lmax;
  unsigned char s[512]; 
  int fErr=noErr;


  if(is_open)
    fErr=fclose(echo_path);
  if(fErr!=noErr)return 0;

  do
    {
      printf("open echo file? y/n\n");
      scanf("%s",fname);
      if(!strcmp(fname,"n"))return 0;
    }
  while(strcmp(fname,"y"));
  printf("what is echo file name ?\n");
  scanf("%s",fname);
  echo_path=fopen(fname,"wb");
  if(!echo_path)return 0;
  n_p_mes=0;
  avpres=0;
  avtemp=0;
  avpot=0;
  nbyte=sprintf(s,
		"     time \t      temperature \t  energy \t       radius \t     pressure \n");
  if(nbyte<=0){ fclose(echo_path);return 0;}
  if(fwrite(&s[0],1,nbyte,echo_path)!=nbyte){fclose(echo_path);return 0;}
  else return 1;
  
}

int set_text_name(int is_open,  char * fname)
{
  int nbyte;
  unsigned char s[512]; 
  int fErr=noErr;

  do
    {
      printf("open text file? y/n\n");
      scanf("%s",fname);
      if(!strcmp(fname,"n"))return 0;
    }
  while(strcmp(fname,"y"));
  printf("what is text file name ?\n");
  scanf("%s",fname);

  text_name=fname;
  return 1;

}


int open_movie_file(int is_open, char * fname)
{
 int fErr=noErr;
 
 if(is_open)
   fErr=closemovie(movie_path);
 if(fErr!=noErr)return 0;
 do
   {
     printf("open movie file? y/n\n");
     scanf("%s",fname);
     if(!strcmp(fname,"n"))return 0;
   }
 while(strcmp(fname,"y"));
 printf("what is movie file name ?\n");
 scanf("%s",fname);
 movie_path=fopen(fname,"wb");
 if(!movie_path)return 0;
 if(write_movie_header(movie_path))return 1;
 return 0;			
}




int write_echo(void)
{ int nbyte;
  unsigned char s[512];
  double time1=get_mes_time();
  double energy=countenergy();
  double temp=n_p_mes ? avtemp/n_p_mes: get_temperature();
  double pot=n_p_mes ? avpot/n_p_mes: get_avePot();
  double gr=get_gr();
  double pressure=n_p_mes ? avpres/n_p_mes: get_pressure();
  /*  printf("%ld\n",n_p_mes);*/
  n_p_mes=0;
  avtemp=0;
  avpres=0;
  avpot=0;
  pot=-pot;
  printf("%lf\n",time1); 
  nbyte=sprintf(&s[0],"%12.3lf\11%12.3lf\11%12.3lf\11%17.10lf\11%17.10lf\n"
		,time1,temp,pot,gr,pressure);
  if(nbyte<=0){ fclose(echo_path);return 0;}
  if(fwrite(&s[0],1,nbyte,echo_path)!=nbyte){fclose(echo_path);return 0;}
  else 
    {
      fflush(echo_path);/*standar C funtion, which forces to write the location where echo_path is pointing of all buffered information*/
      return 1;
    }
}



void set_time(double time1)
{
  timed=time1;
}

double get_time(void)
{
  return timec+timed;
}

void set_frate(double frate)
{
  delta1=frate;
  delta2=0;
}
double get_frate(void)
{return delta1;}
void set_mfrate(double frate)
{
delta3=frate;
delta4=0;
}
double get_mfrate(void)
{return delta3;}


void vp(crd * a, crd * b, crd * c)
{
  c->x=a->y*b->z-a->z*b->y;
  c->y=a->z*b->x-a->x*b->z;
  c->z=a->x*b->y-a->y*b->x;
}

double dist(crd r,crd s)/*return square of the cartesian distance*/
{
  double x1,y1,z1;
  x1=r.x-s.x;
  y1=r.y-s.y;
  z1=r.z-s.z;
  return(y1*y1+x1*x1+z1*z1);
}


/*In reaction, we asses whether reaction successful. If so, then break the 
  existing bond or set a new bond. Also, since atoms change their types, we 
  need to update the collision types between "p" and the neighbors of "p", and
  between "q" and the neighbors of "q". If say, "p" and neighbor "ip" were
  bonded before p-q reaction, it may be that they can not be bonded after p-q
  reaction. In that case we may break the bond*/
int reaction (atom *a,int i1,int i2, int ct1,double sc,
	      double x,double y,double z)                 /*208 lines of code*/
{  
  int ct=ct1;
  int rtype;
  int revers;
  if(sc<0)   /*atoms are approaching to each other, the bond is to be created*/
    {
      rtype=coll[ct1].react;  /*Has potential barrier an associated reaction?*/
      if(rtype<=0)return -1;                             /*Return no reaction*/
      revers=0;        /*There's reaction. revers==0 means bond to be created*/
    }
  else/*atoms are approaching to each other,the bond is to be broken. We take*/
    {/*the external pot. barrier,which is "prev" barrier to the one stored   */
      ct=coll[ct1].prev;
      rtype=~coll[ct].react; 
      if(rtype<=0)return -1; 
      revers=1;         /*There's reaction. revers==1 means bond to be broken*/
    }
  {
    int out,i,j,ix,iy,iz;
    double ab1,ab2,vx,vy,vz,ab,di,ed;
    atom *a1,*a2;
    double old_pot=coll[ct1].etot;                   /*Energy before reaction*/
    double du,duc,new_pot=0;                 /*new_pot, energy after reaction*/
    int ct_new;
    int np=get_np();  /*number atoms with previous scheduled collision with p*/
    int nq=get_nq();  /*number atoms with previous scheduled collision with p*/
    int * ap=get_atomp();    /*list of atom numbers with previous p-collision*/
    int * aq=get_atomq();    /*list of atom numbers with previous q-collision*/
    int * cp=get_collp();    /*collision types of atoms with prev p-collision*/
    int * cq=get_collq();    /*collision types of atoms with prev q-collision*/
    int bond,iq,ip;
    int old1,old2,new1,new2;  
    double m1,m2;
    a1=a+i1;                                               /*points to atop p*/
    a2=a+i2;                                               /*points to atop q*/
    old1=a1->c;                       /*atom type of atom "p" before reaction*/
    old2=a2->c;                       /*atom type of atom "q" before reaction*/
    m1=a1->m;                                              /*mass of atom "p"*/
    m2=a2->m;                                              /*mass of atom "q"*/
    if(!revers)           /*The bond will be created: old1+old2 --> new1+new2*/
      {
	if(old1==react[rtype].old1)     /*old1 is type of "p" before bonding*/
	  {
	    new1=react[rtype].new1;      /*new1 is type of "p" after bonding*/ 
	    new2=react[rtype].new2;      /*new2 is type of "q" after bonding*/ 
	  }
	else /*If type of p listed as react[rtype].old2,then new1 still will*/
	  {  /*be type of p after bonding, but refers to react[rtype].new2 */
	    new1=react[rtype].new2;     /*new1 is type of "p" after bonding*/
	    new2=react[rtype].new1;     /*new2 is type of "q" after bonding*/
	  }
      }
    else                   /*The bond will be broken: old1+old2 <-- new1+new2*/
      {/*By convention, react[rtype].new1 and react[rtype].new2 are types of */
	if(old1==react[rtype].new1)          /*atoms when in the bonded state*/
	  {
	    new1=react[rtype].old1; /*new1, type of p after breaking the bond*/
	    new2=react[rtype].old2; /*new2, type of q after breaking the bond*/
	  }
	else 
	  {
	    new1=react[rtype].old2; /*new1, type of p after breaking the bond*/
	    new2=react[rtype].old1; /*new2, type of q after breaking the bond*/
	  }
      }


    if(revers)                                      /*The bond will be broken*/
      ct_new=react[rtype].out;/*react[rtype].out, coll. type of unbound pair */
    else                                           /*The bond will be created*/
      ct_new=react[rtype].in;  /*react[rtype].out, coll. type of bonded pair */


    new_pot+=coll[ct_new].etot;                       /*Energy after reaction*/
    if(old1!=new1)         /*the atomic type of "p" has changed upon reaction*/
      {
	a1->c=new1;                             /*Update the type of atom "p"*/
	for(i=0;i<np;i++)/*Go through atoms with previous scheduled collision*/
	  {
	    ip=ap[i];                    /*Atom number of one the "neighbors"*/
	    if(ip!=i2)                    /*Make sure the neighbor is not "q"*/
	      { /*cp[ip] is collision type of previously scheduled p-ip coll.*/
		old_pot+=coll[cp[ip]].etot;          /*energy before reaction*/
		moveatom(a+ip);     /*update position of "ip" to system time */
		bond=is_bond(cp[ip]);        /*if "p" and "ip" bonded, bond=1*/
		collp[ip]=after_type(i1,ip,&bond,cp[ip]);     /*new coll.type*/
		if(collp[ip]<0)return -1;
		new_pot+=coll[collp[ip]].etot;        /*Energy after reaction*/
		if(bond)collp[ip] = ~(collp[ip]);/*If "p" and "ip"were bonded*/
	      } /*before "p" underwent reaction,but then the bond was broken */
	  } /*in after_tupe(..),then bond==1 and we take the bitwise         */
      }/*complement of collp[ip], which is -collp[ip]-1, a negative number   */
    if(old2!=new2)         /*The atomic type of "q" has changed upon reaction*/
      {
	a2->c=new2;                             /*Update the type of atom "q"*/
	for(i=0;i<nq;i++)/*go through atoms with previous scheduled collision*/
	  {
	    iq=aq[i];             /*Atom number of one the "neighbors" of "q"*/
	    if(iq!=i1)                    /*make sure the neighbor is not "p"*/
	      {
		old_pot+=coll[cq[iq]].etot;          /*Energy before reaction*/
		moveatom(a+iq);     /*Update position of "iq" to system time */
		bond=is_bond(cq[iq]);        /*If "q" and "iq" bonded, bond=1*/
		collq[iq]=after_type(i2,iq,&bond,cq[iq]);     /*new coll.type*/
		if(collq[iq]<0)return -1;
		new_pot+=coll[collq[iq]].etot;        /*energy after reaction*/
		if(bond)collq[iq]=~(collq[iq]); /*We remember bond iis broken*/
	      }              /*by storing a negative number as collision type*/
	  }
      }
    du=new_pot-old_pot;            /*Change in total energy upon p-q reaction*/
    if(!react[rtype].bond)du+=react[rtype].eo;/*If bonding,release latentheat*/
    duc=du*corr_2;           /*Adjusted change in energy for the thermal bath*/
    if(m1==DBL1)
      ed=2*duc*coll[ct].dd/m2;
    else if(m2==DBL1)
      ed=2*duc*coll[ct].dd/m1;
    else
      ed=2*duc/(m1*m2*coll[ct].dm);
    
    di=1.0+ed/(sc*sc);
    if(di<=0)  /*When ed is large negative, unsuccessfull attempt to escape: */
      {	       /*reaction does not happen                                    */
	ab=-2.0*sc/coll[ct].dd;
	a1->c=old1;  /*Return the old types of a1 and a2 to the state before */
	a2->c=old2;  /*We evaluated reaction(..)                             */
	if(!revers)return -1;       /*revers==0 if the bond was to be created*/
	ct_new=ct1;   /*Return the old collision type between "i" and "k"    */
      }
    else                      /*di>0, thus reaction is energetically possible*/
      {
	ab=sc*(sqrt(di)-1.0)/coll[ct].dd;
	vvm+=duc;
	potential+=du;
	if((react[rtype].old1!=react[rtype].new1) ||
	   (react[rtype].old2!=react[rtype].new2)    )
	  setNewTypes(1);                 /*Set flag "newTypes" to value of 1*/
	if(revers)                            /*The reaction breaks the bond */
	  breakBond(i1,i2);
        else if(react[rtype].bond)             /*The reaction produces a bond*/
	  setBond(i1,i2);
	if(new2!=old2)       /*If the type of "q" has changed in the reaction*/
	  for(i=0;i<nq;i++)               /*Go through the "neighbors" of "q"*/
	    {
	      iq=aq[i];                 /*Atom number of one of the neighbors*/
	      if(iq!=i1)                  /*Make sure the neighbor is not "p"*/
		{
		  if(collq[iq]<0)/*Negative if the bond q-iq broke because of*/
		    {                                      /*the reaction p-q*/
		      breakBond(i2,iq);
		      cq[iq]=~collq[iq];/*Insert the correct (positive) type */
		    }          /*between "q" and "iq" into "search.collq[ip]"*/
		  else
		    cq[iq]=collq[iq];   /*Either the bond q-iq was conserved */
		}  /*after reaction or there was no bond between "q" and "iq"*/
	    }
	if(new1!=old1)       /*If the type of "p" has changed in the reaction*/
	  for(i=0;i<np;i++)               /*Go through the "neighbors" of "q"*/
	    {
	      ip=ap[i];                 /*Atom number of one of the neighbors*/
	      if(ip!=i2)                  /*Make sure the neighbor is not "q"*/
		{
		  if(collp[ip]<0)/*Negative if the bond p-ip broke because of*/
		    {                                      /*the reaction p-q*/
		      breakBond(i1,ip);
		      cp[ip]=~collp[ip];/*Insert the correct (positive) type */
		    }          /*between "p" and "ip" into "search.collp[ip]"*/
		  else
		    cp[ip]=collp[ip];   /*Either the bond p-ip was conserved */
		}  /*after reaction or there was no bond between "q" and "iq"*/
	    }
	cp[i2]=ct_new;  /*insert the new collisiont type between "p" and "q" */
      }                                                      /*into collp[i2]*/
    if(coll[ct].mij!=DBL1){
      ab1=ab*m2/coll[ct].mij;
      ab2=-ab*m1/coll[ct].mij;
      virial+=ab1*m1*coll[ct].dd;
      a1->v.x+=x*ab1;
      a1->v.y+=y*ab1;
      a1->v.z+=z*ab1;
      a2->v.x+=x*ab2;
      a2->v.y+=y*ab2;
      a2->v.z+=z*ab2;
    }
    else{
      if(m1==DBL1){
	virial+=ab*m2*coll[ct].dd;
	a2->v.x-=x*ab;
	a2->v.y-=y*ab;
	a2->v.z-=z*ab;
      }
      else if(m2==DBL1){
	virial+=ab*m1*coll[ct].dd;
	a1->v.x+=x*ab;
	a1->v.y+=y*ab;
	a1->v.z+=z*ab;
      }
    }
    return ct_new; 
  }
}/*Matches int reaction(...)*/

int newvel (atom *a,int i,int j, int ct1)
{  
  int out,ix,iy,iz;
  double ab1,ab2,vx,vy,vz,x,y,z,ab,sc,di,ed;
  atom *a1,*a2;
  int k, ct=ct1;/*collision type between "i" and "j" before reaction*/
  a1=a+i;
  a2=a+j;
  /*  if(ct1==18)
       if(((i==346)&&(j==347))||((i==347)&&(j==346)))
       printf("error\n"); */
  moveatom(a1);
  moveatom(a2);
  vx=a1->v.x-a2->v.x;
  vy=a1->v.y-a2->v.y;
  vz=a1->v.z-a2->v.z;
  x=a1->r.x-a2->r.x;
  y=a1->r.y-a2->r.y;
  z=a1->r.z-a2->r.z;
  ix=a1->i.x.i-a2->i.x.i;
  iy=a1->i.y.i-a2->i.y.i;
  iz=a1->i.z.i-a2->i.z.i;
  if (ix>1)x-=bound[0].length;
  else if (ix<-1)x+=bound[0].length;
  if (iy>1)y-=bound[1].length;
  else if (iy<-1)y+=bound[1].length;
  if (iz>1)z-=bound[2].length;
  else if (iz<-1)z+=bound[2].length;
  sc=vx*x+vy*y+vz*z;
  
  if(nrt)/*if reactions are possible in the system*/
    {
      int ctr;      
      ctr=reaction(a,i,j,ct1,sc,x,y,z);
      if(ctr>-1){return ctr;}/*no reaction means ctr==-1*/      
    }
  
  k=ct;/*k is collision type between "i" and "j" before reaction*/
  ed=coll[k].edm;
  /*  k=ct&UNSTABLE;  type of collision */
  
  if(sc>=0)
    {
      /*if sc>=0 the atoms are moving away, we should take outer parameters of the outer well */
      if(k>=ntot)
	{
	  num_gaps--;
	  return coll[k].prev;
	}
      k=coll[k].prev;
      ed=-coll[k].edm;  /* depth of potential well */
      /* attempt to escape from the well; external collisions are 
	 sometimes with finite energy; then the energy is taken with
	 negative sign in oppose to the case when the attoms are jumping 
	 in the well. coll[k].e is the depth of well is positive for
	 attraction */ 
    }
  else if(coll[k].prev>=ntot)
    {
      num_gaps--;
      return coll[k].next;
    }
  
  if(coll[k].e==-dblarg1)
    {
      /* energy as dblarg1 always means ellastic repulsion */
      ab=-2.0*sc/coll[k].dd;
      /*      if(sc>=0.0) 
	      ct|=STABLE; */
    }
  else 
    {
      di=1.0+ed/(sc*sc);
      if(di<=0)
	{
	  
	  /*	when ed is large negative, it is
		unsuccessfull attempt to escape: 
		ellastic collision */   
	  ab=-2.0*sc/coll[k].dd;
	  /*	  if(sc>=0.0) 
		  ct|=STABLE; */
	}
      else
	{
	  ab=sc*(sqrt(di)-1.0)/coll[k].dd;
	  if (sc>=0.0) 
	    {
	      /*	   atoms jumped out of the well and move to a previous well*/  
	      ct=k;
	      vvm-=coll[k].e;
	      potential-=coll[k].eo;
	    }
	  else 
	    {
	      /* atoms jumed into the next well */
	      vvm+=coll[k].e;
	      potential+=coll[k].eo;
	      ct=coll[k].next;
	    }
	}
    }
  if(coll[k].mij!=DBL1){
    ab1=ab*a2->m/coll[k].mij;
    ab2=-ab*a1->m/coll[k].mij;
    virial+=ab1*a1->m*coll[k].dd;
    a1->v.x+=x*ab1;
    a1->v.y+=y*ab1;
    a1->v.z+=z*ab1;
    a2->v.x+=x*ab2;
    a2->v.y+=y*ab2;
    a2->v.z+=z*ab2;
  }
  else{
    if(a1->m==DBL1){
      virial+=ab*a2->m*coll[k].dd;
      a2->v.x-=x*ab;
      a2->v.y-=y*ab;
      a2->v.z-=z*ab;
    }
    else if(a2->m==DBL1){
      virial+=ab*a1->m*coll[k].dd;
      a1->v.x+=x*ab;
      a1->v.y+=y*ab;
      a1->v.z+=z*ab;
    }
  }
  return ct; 
  
}

/* i1 atom number
   j1 is wall number :
      n+1 for x
      n+2 for y
      n+3 for z   */
newloc (atom *a,int i1,int j1)
{  
  int xy,i;
  atom *a1;
  triad *b;
  double *aa;
  dimensions *bound1;
  int address,step,period;
  a1=a+i1; /*pointer to atom i1*/
  moveatom(a1); /*update position of atom i1, thus it crosses to a new cell*/
  aa=(double *)a1;
  xy=j1-n1; /*xy=0 means X-axis, xy=1 means Y-axis, xy=2 means Z-axis*/
  aa+=xy;
  b=(triad *)aa;/*if xy=1, then b stores Y-coord., Y-vel. and Y-cell index*/
  bound1=&bound[xy];
/* take the old box number of the atom 
and decrease or increase it accordingly */
  i=b->i.i;/*if xy=1, then "i" is old cell index of the atom*/
  if(b->v>0)
    {
      i++;
      if (i==bound1->period)/*check for periodic boundary conditions*/
	{
	  i=0;
	  b->r-=bound1->length;
	}
    }
  else
    {      
      i--;
      if (i==-1)
	{
	  i+=bound1->period;
	  b->r+=bound1->length;
	}	
    }
  b->i.i=i;/*if xy=1, assing new cell index along Y-axis*/
  /*address now is the address of the cell where the particle is located,    */
  /*however, we have not updated this quantity in the atom structure of i1,  */
  /*i.e., a[i1]->add still contains the address of the old cell              */
  address=a1->i.x.i+bound[0].period*(a1->i.y.i+bound[1].period*a1->i.z.i);
  find_atoms(i1,address);
 }

int twall(int i, double * t1)
{
  double s, x ,d ,y,z,rx ,ry,rz,vx,vy,vz,hry,wrx,drz,vv;
  atom *pt;
  double t;
  int q;
  pt=a+i;/*pt points to info on atom number i*/
  q=n3;/*walls on z-axis*/
  x=dblarg1;/*x is time will the particle to leave the cell if
its velocity would be (vx,0,0). But we first initialize to dblarg1*/
  y=dblarg1;
  z=dblarg1;

  rx=pt->r.x;/*position of atom*/
  ry=pt->r.y;
  rz=pt->r.z;

  vx=pt->v.x;/*velocity of atom*/
  vy=pt->v.y;
  vz=pt->v.z;

  wrx= pt->i.x.i * bound[0].dl;/*pt->i.x.i is cell index of atom
i along x-axis, thus wrx is distance to the origin on the x-axis
of one of the two walls enclosing atom i in the x-axis. The wall
is actually the one which is closer to the origin, since indexes
begin with zero
------------------------>(this is rx vector)
|-------------|---------X----|-------------|---(these are the cells)
-------------->(this is wrx vector)
*/
  hry=pt->i.y.i*bound[1].dl;
  drz=pt->i.z.i*bound[2].dl;

  rx-=wrx;
  ry-=hry;
  rz-=drz;

  wrx=bound[0].dl-rx;  /*(----------)this is rx       */
  hry=bound[1].dl-ry;  /*|----------X----|            */
  drz=bound[2].dl-rz;  /*           (----)this is wrx */
  
  if (vx<0) 
    x=-rx/vx; /*time to hit along the x-axis the closer wall to the origin*/
  else if (vx>0)
    x=wrx/vx;/*time to hit along the x-axis the further wall to
		the origin*/
  if (vy<0)
    y=-ry/vy;
  if (vy>0)
    y=hry/vy;
  if (vz<0)
    z=-rz/vz;
  if (vz>0)
    z=drz/vz;

  t=z;

  if ((x<z)||(y<z))
    {
      t=y;
      q=n2;
      if(x<y)
	{
	  t=x;
	  q=n1;
	}
    }
  if(t<dblarg1) t+=pt->t;/*when initialized, pt->t is zero(in 
make_system.c, line sam[i].t=0.0; in make_key_system(...) )*/
  else t=dblarg1;

  pt->w=t;/*in pt->w we store the time to hit the closest wall
	    of the enclosing cell*/
  *t1=t;
  return q;
}


int tball(int i,int j,int ct, double * t1)
{
  int k,ix,iy,iz;

  atom *a1,*a2;
  double t=dblarg1;/*assume that particles will not collide*/
  int q=0;
  double ab,sc,di,de,x,y,z,u,dd,v,w,dt,delta;
  a1=a+i;
  a2=a+j;
  k=ct;
  u=a2->v.x-a1->v.x;
  v=a2->v.y-a1->v.y;
  w=a2->v.z-a1->v.z;
  ab=u*u+v*v+w*w;
  if(ab) /*two static particles have ab==0*/
    {
      x=a2->r.x-a1->r.x;
      y=a2->r.y-a1->r.y;
      z=a2->r.z-a1->r.z;
      delta=a1->t-a2->t;/*the last times when the position of the atoms where updated will differ*/
      if(delta>0)/*particle a1 was updated more recently than particle a2*/
	{
	  x+=delta*a2->v.x;
	  y+=delta*a2->v.y;
	  z+=delta*a2->v.z;
	}
      else
	{
	  x+=delta*a1->v.x;
	  y+=delta*a1->v.y;
	  z+=delta*a1->v.z;
	}
      ix=a2->i.x.i-a1->i.x.i;
      iy=a2->i.y.i-a1->i.y.i;
      iz=a2->i.z.i-a1->i.z.i;
      if (ix>1)x-=bound[0].length;
      if (ix<-1)x+=bound[0].length;
      if (iy>1)y-=bound[1].length;
      if (iy<-1)y+=bound[1].length;
      if (iz>1)z-=bound[2].length;
      if (iz<-1)z+=bound[2].length;

      sc=u*x+v*y+w*z;
      de=(x*x+y*y+z*z)*ab-sc*sc;
        
      if (sc<0.0)
	{
	  di=ab*coll[k].dd-de;
	  if (di>0.0)
	    t=(-sc-sqrt(di))/ab;
	}
      if((t==dblarg1)&&(coll[k].prev>-1))
	{ 
	  t=(-sc+sqrt(ab*coll[coll[k].prev].dd-de))/ab;
	}   
 
    }
  if (t<dblarg1)
    {
      if(delta>0)
	t+=a1->t;
      else
	t+=a2->t;
    }
  if((coll[k].prev>-1)||((t<=a1->w)&&(t<=a2->w))) q=1;
  *t1=t;
  return q;
}



int collision_type(int i, int k)
{  
  int ic=a[i].c;/*type of atom i                                             */
  double rx,ry,rz,dr;
  int ia,ky,kz;
  int link_err=isFriend(k,i);/*return one if i and k are linked              */
  int kc=a[k].c;       /*type of atom k                                      */
  int ie=ecoll[ic][kc];/*Index on array coll where info on non-ellastic      */
  int ix=a[i].i.x.i;   /*collision for atom types ic and kc begins. ie is the*/
  int iy=a[i].i.y.i;   /*collision type for the most external collision      */
  int iz=a[i].i.z.i;
  int kx=ix-a[k].i.x.i;/*Difference in cell index along x-axis               */
  int ct=ie;           /*Initialize the collision type as the most exterior  */
  ia=abs(kx);          /*collision                                           */
  if((ia>1)&&(ia!=px))goto far_away;/*if the difference in cell              */
             /*indexes along x-axis between cell where atom i is and cell    */
             /*where atom j is, is bigger than one (not neighboring cell) and*/
             /*different than the maximum cell number along the x-axis, then */
             /*goto far_away                                                 */
  ky=iy-a[k].i.y.i;
  ia=abs(ky);
  if((ia>1)&&(ia!=py))goto far_away;
  kz=iz-a[k].i.z.i;
  ia=abs(kz);
  if((ia>1)&&(ia!=pz))goto far_away;
  rx=a[i].r.x-a[k].r.x;
  ry=a[i].r.y-a[k].r.y;
  rz=a[i].r.z-a[k].r.z;
  if(kx<-1)rx+=lx;     /*Periodic boundary conditions                        */
  if(kx>1)rx-=lx;
  if(ky<-1)ry+=ly;
  if(ky>1)ry-=ly;
  if(kz<-1)rz+=lz;
  if(kz>1)rz-=lz;
  dr=rx*rx+ry*ry+rz*rz;/*Distance between the two atoms                      */

  if(link_err)         /*Atoms i and k are linked,thus ie isn't suitable type*/
    {
      int ii=icoll[ic][kc];/*icoll for linked atoms, ii is the index in array*/
                       /*coll where info for linked types ic and kc begins   */
      if(ii>-1)        /*-1 means we hit the interior hard core repulsion    */
	{
	  if(dr<coll[ii].dd)
	    { 
	      for(ii=coll[ii].next; ii>-1; ii=coll[ii].next)
		{     /*Go to the next interior collision                    */
		  if(dr>coll[ii].dd)
		    {
		      link_err=0;
		      ct=ii; /*We found the type, now return it              */
		      goto far_away;
		    }
		}
	    }
	}
    }
/*beginning with the most external collision, we are going towards deeper collisions, in case particles are very close. When we reach the interior
hard-core repulsion, ie==-1*/
  while(ie>-1)           /*Beginning with the most external collision, we are*/
    {			 /*going towards deeper collisions,in case particles */
      if(dr>=coll[ie].dd)/*are very close. When we reach interior hard-core  */
	{		 /*repulsion, ie==-1, we stop                        */
	  ct=ie;         /*We found the type, now return it                  */
          goto far_away;
	}
      ie=coll[ie].next; /*Go to the next interior collision                  */
    }
  too_close_dialog(k,i,sqrt(dr));/*atoms k and i are too close               */
  return -3;
 far_away:
  if(link_err){         /*if linked and the cells are not neighbors,then link*/
                        /*must be broken!                                    */
   bond_error_dialog(k,i,sqrt(dr));breakBond(i,k);
 }
  return ct;
}	


/*atoms "i" and "j" are undergoing a reaction. How are third party interactions
  affected? Here, "i" and "k" had a previously scheduled collision with 
  collision type old_ct. "link_err" signals wheter they were bonded before "i"
  underwent reaction with "j". After i-j reaction, what is the new collision
  type between "i" and "k" ?
  When we invoke after_type, "i" has already changed type due to i-j reaction.
  Also, the position of "i" and "k" have been updated to system type, so that 
  their local times coincide with the system time                           */
int after_type(int i, int k, int * link_err,int old_ct)  /*72 lines of code */
{          /*a[i].c is type of "i" after reaction (we already updated it in */
  int ic=a[i].c;                                      /*"reaction" function */
  int prev=coll[old_ct].prev;    /*Collision type of inner potential barrier*/
  double rx,ry,rz,dr;
  int ia,ky,kz;
  int kc=a[k].c;      /*Type of "k". Remember "k" is not undergoing reaction*/
  int ie=ecoll[ic][kc];/*Most external collision type between unbound "i"   */
  int ix=a[i].i.x.i;   /*and "k", if "i" and "k" were unbound (link_err==0) */
  int iy=a[i].i.y.i;
  int iz=a[i].i.z.i;
  int kx=ix-a[k].i.x.i;
  int ct=ie;
  ia=abs(kx);
  if((ia>1)&&(ia!=px))return ct; /*"i" and "k" were not bonded because they */
  ky=iy-a[k].i.y.i;     /*were separated by more than one cell in the x-axis*/
  ia=abs(ky);
  if((ia>1)&&(ia!=py))return ct;
  kz=iz-a[k].i.z.i;
  ia=abs(kz);
  if((ia>1)&&(ia!=pz))return ct;
  rx=a[i].r.x-a[k].r.x;
  ry=a[i].r.y-a[k].r.y;
  rz=a[i].r.z-a[k].r.z;
  if(kx<-1)rx+=lx;                  /*deal with periodic boundary conditions*/
  if(kx>1)rx-=lx;
  if(ky<-1)ry+=ly;
  if(ky>1)ry-=ly;
  if(kz<-1)rz+=lz;
  if(kz>1)rz-=lz;
  dr=rx*rx+ry*ry+rz*rz;                        /*distance between "i" and "k"*/
  /*because of floating point errors, it could be that "dr" smaller than
    coll[old_ct].dd. Ideally, dr>coll[old_ct].dd allways because "i" and "k" 
    have not collided yet at current system time. If it is the hard-core
    potential barrier (prev==-1), then "dr" should be bigger than
    coll[old_ct].dd. If there is a more internal potential barrier (prev>-1),
    then let "dr" be smaller than coll[old_ct].dd */
  if(dr<=coll[old_ct].dd)
    {
      dr=next_double(coll[old_ct].dd);
     }
  if(prev>-1)/*there is a more internal collision between "i" and "k" (with  */
    {        /*the type of "i" that one before "i" underwent collision       */
      if(dr>=coll[prev].dd)
	dr=prev_double(coll[old_ct].dd);
     }

  /*link_err==1 if "i" and "k" were bonded before "i" underwent reaction, but
    in this context we will use it to asses wheter the bond survives after
    "i" undergoes reaction                                                   */
  if(*link_err)    
    {
      int ii=icoll[ic][kc];/*Most external collision type between "i" and "k"*/
      /*If there's no bonding defined between the new type of "i" and that of
        "k", then ii==-1 */
      if(ii>=ntot)ii=coll[ii].next;/*due to variable volume "gap"*/
      if(ii>-1)     /*bond between type of "k" and new type of "i" is defined*/
	{
	  if(dr<coll[ii].dd)/*"i" and "k" distance bigger than distance corr-*/
	    { /*-esponding to most stretched bond.                           */
	      for(ii=coll[ii].next;ii>-1;ii=coll[ii].next)/*navigate towards */
		{                   /*more internal potential energy barriers*/
		  if(dr>coll[ii].dd)/*Found an potential barrier with smaller*/
		    {                                    /*distance than dr. */
		      *link_err=0;/*Bond survives after "i" undergoes        */
		      ct=ii;      /*reaction,   no link error                */
		      return ct;
		    }
		}/*if dr bigger than distance corresponding to most compres- */
	    }    /*-sed bond, then i-k bond will not survive (link_error==1) */
	}        /* and we turn now to the collision type between "i" and "k"*/
    }            /* when "i" and "k" are not bound                           */

  while(ie>-1)   /*if "i" and "k" were unbound before "i" underwent reaction,*/
    {            /*find unbound collision type between "i" and "k" after "i" */
      if(dr>=coll[ie].dd)        /*undergoes reaction (and changed atom type)*/
	{
	  ct=ie;
          return ct;
	}
      ie=coll[ie].next;
    }
  return -1;/*if dr below the hard-core distance between type of "k" and new*/
            /*type of "i"                                                   */
}/*Matches int after_type(int i, int k, int * link_err,int old_ct)*/



moveatom( atom *pt)
{ double delta=timeb-pt->t;
  if(delta)
  {
      pt->r.x=pt->r.x+pt->v.x*delta;
      pt->r.y=pt->r.y+pt->v.y*delta;
      pt->r.z=pt->r.z+pt->v.z*delta;
      pt->t=timeb;
   }
}
/*store positions of all atoms for time==timeb in structure pt->q for
  future use*/
void moveatoms(void)
{ double delta;
  atom *pt;
  for (pt=a;pt->c!=0;pt++)
    { 
      delta=timeb-pt->t;    
      pt->q.x=pt->r.x+pt->v.x*delta;
      pt->q.y=pt->r.y+pt->v.y*delta;
      pt->q.z=pt->r.z+pt->v.z*delta;
    }
}

/*move positions of all atoms to positions corresponding to time==timeb*/
void update_atoms(void)
{ double delta;
  atom *pt;
  for (pt=a;pt->c!=0;pt++)
    { 
      delta=timeb-pt->t;    
      pt->r.x=pt->r.x+pt->v.x*delta;
      pt->r.y=pt->r.y+pt->v.y*delta;
      pt->r.z=pt->r.z+pt->v.z*delta;
      pt->t=timeb;
    }
}
/*scale all potential walls*/
void reset_colldata(void)
{
  int i;
      for (i=0;i<nen;i++)
	{
	  coll[i].e=coll[i].eo;
	  coll[i].edm=coll[i].edmo;
	}
      return;
}

void corr_vel(void)
{
  int i;
  double corr1=1/corr;
  for( i=0;i<n1;i++)
    {
      a[i].u.x=a[i].v.x*corr1;
      a[i].u.y=a[i].v.y*corr1;
      a[i].u.z=a[i].v.z*corr1;
    }
}

double countenergy(void)
{
  return -potential;
}

double get_pressure(void)
{ 
  if (pressure!=dblarg1)
  return pressure;
  else if(timep)
  return (virial+vvmtime+vvmtime)/(volume*dim*timep*corr_2);
  else 
  return dblarg1;
}

int collision(){/*76 lines of code*/
  int i,k;
  int coll_type=0;
  double vvm0,t2,corrt1;
  t2=timea-timeb;/*increase in time from the last computed collision*/
  if(t2<0)
    printf("error %lf %d %d\n",timeb,p1,q1);
  timeb=timea;
  ts++;
  corrt1=t2*corr;
  timec+=corrt1;
  timep+=t2;
  vvmtime+=t2*vvm;
  potTime+=t2*potential;
  delta2+=corrt1;
  delta4+=corrt1;
  if(t_is_open)delta6+=corrt1;
  if (q1>=n1) newloc(a,p1,q1);
  else 
    {
      coll_type=1;
      ll++;
      llp++;
      ct1=newvel (a,p1,q1,ct1); /*virial is computed inside newvel*/
      if(llp==deltall)
	{
	  pressure=(virial+vvmtime+vvmtime)/(volume*dim*timep*corr_2);
	  temperature=2*(vvmtime/timep)/(n1*dim*corr_2);/*since "vvmtime/timep" is the average kinetic energy during a "timep" lapse of time, temperature is the average temperature during this lapse. Coefficient corr_2 rescales the kinetic energy of the system that is being simulated to the kinetic energy of the system that would result if we really changed the temperature of the particles. */
          avePot=potTime/timep;
	  mes_time=timec;
	  if(is_open)
	    {
	      n_p_mes++;
	      avpres+=pressure;
	      avtemp+=temperature;
              avpot+=avePot;
	    }
	  if(m_is_open) add_movie_param(temperature,-avePot,temperature*0.5*dim*n1-avePot,pressure);
	  if(fs_ok()) fs_add(temperature,-avePot,temperature*0.5*dim*n1-avePot,pressure); 
	  rescale();
	  if(gap&&(var_density||press_coeff))
	    {
	      n_gap_mes++;
	      press_curr+=pressure;
	      if((n_gap_mes>=max_gap_mes)&&(!num_gaps))
		new_density=1;
	    }
	  if(schedule_ok())read_schedule();

	}
    }
  
  if (delta2>delta1)
    {
      delta2-=delta1;
      if(is_open)
	{
	  is_open=write_echo();
	}     
    }
  if (delta4>delta3)
    {
      delta4-=delta3;
      if(m_is_open)m_is_open=write_movie_frame();
    }
  if (delta6>delta5)
    {
      if(!coll_type)
	{
	  delta6-=delta5;
	  writetext(text_name);
	}
    }

  corr_func(corrt1);
  return coll_type;
}/*Matches int collision()*/


int writetext (char* fname)
{
  char * dig="0123456789";
  char * newname;
  int fErr=0;
  FILE * path;
  int i,k,j=text_no;
  int name_length=strlen(fname);
  newname=malloc(name_length+5);
  for(i=0;i<name_length;i++)
    newname[i]=fname[i];
  name_length+=4;
  newname[name_length]=(char)0;
  for(i=name_length-1;i>=name_length-4;i--)
    {
      k=j % 10;
      newname[i]=dig[k];
      j=j/10;
    }
  path=fopen(newname,"wb");
  free(newname);
  if(!path)return 1;
  fErr=write_key_coord(path);
  if(fErr==noErr)
    {
      fflush(path);
      fclose(path);
      text_no++;
      return 0;
    }
  else 
    return 1;			
}

/*funcion stop_atoms put the center of mass in (0,0,0), eliminates to total momentum of the system*/
void stop_atoms(moved_iatom * a, int n)
{    
  double  rx,sx,dx,sm;
  int i,j,k;
  double A[3][3];
  double W[3],W1[3],M[3];
  double I,smax,norm;
  double corr1=1/corr;
  sm=0;           /*sm, system mass*/
  for(i=0;i<n;i++)
    sm+=a[i].m;
  
  for(j=0;j<3;j++)              
    {
      rx=a[0].r[j];/*store in rx the "j" componento of first atom position*/
      a[0].r[j]=0;/*put first atom at origin of coordinates*/
      sx=0;       /*Current component (given by "j") of the center of mass*/
      for(i=1;i<n;i++){ 
	dx=a[i].r[j]-rx;
	/*next a trick: dx+dx<-L cheaper to compute than dx<-L/2 */
	if(dx+dx<-bound[j].length)dx+=bound[j].length;
	if(dx+dx>bound[j].length)dx-=bound[j].length;
	rx=a[i].r[j];
	a[i].r[j]=a[i-1].r[j]+dx;  /*Eliminate per.bound.cond. of coordinates*/
	sx+=a[i].r[j]*a[i].m;
      }
      sx/=sm;
      for(i=0;i<n;i++)
	a[i].r[j]-=sx;         /*We will place the center of mass in (0,0,0)*/
    }

    for(j=0;j<3;j++)
      {
	sx=0;      /*Current component (given by "j") of the system momemtum*/
	for(i=0;i<n;i++){ 
	  sx+=a[i].v[j]*a[i].m;
	}
	sx/=sm;
	for(i=0;i<n;i++)
	  a[i].u[j]=(a[i].v[j]-sx)*corr1;   /*Total momentum will be (0,0,0)*/
      }

  I=0;
  for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
	{
	  double  rij=0;
	  for(k=0;k<n;k++)
	    {
	      rij+=a[k].r[i]*a[k].r[j]*a[k].m;
	    }
	  A[i][j]=rij;                                /*A[][], inertia tensor*/
	}
      I+=A[i][i];
      M[i]=0;
    }
  if(!I) return;
  for(k=0;k<n;k++)                              /*M[], total angular momentum*/
    {
    M[2]+=(a[k].r[0]*a[k].u[1]-a[k].r[1]*a[k].u[0])*a[k].m;
    M[0]+=(a[k].r[1]*a[k].u[2]-a[k].r[2]*a[k].u[1])*a[k].m;
    M[1]+=(a[k].r[2]*a[k].u[0]-a[k].r[0]*a[k].u[2])*a[k].m;
  }
    printf("%lf %lf %lf\n",M[0],M[1],M[2]);
norm=0;
  for(i=0;i<3;i++)
    {
      for(j=0;j<3;j++)
	{
	  A[i][j]/=I;
	}
      M[i]/=I;
      W[i]=M[i];
      norm+=M[i]*M[i];
    }
  k=0;
  norm*=1.0e-32;
  do{
    smax=0;
    for(i=0;i<3;i++)
      {
	double s=0;
	for(j=0;j<3;j++)
	  s+=A[i][j]*W[j];
	W1[i]=s;
      }
    for(i=0;i<3;i++)
      {double s1=W1[i]+M[i];
       double s=s1-W[i];
       smax+=s*s; 
      W[i]=s1;
     }
    /*    printf("%le\n",smax);*/
    k++;
if(k==1000)break;
  }while(smax>norm);


    for(k=0;k<n;k++)
      {
	a[k].u[2]-=W[0]*a[k].r[1]-W[1]*a[k].r[0];
	a[k].u[0]-=W[1]*a[k].r[2]-W[2]*a[k].r[1];
	a[k].u[1]-=W[2]*a[k].r[0]-W[0]*a[k].r[2];
      }

    /*
  for(j=0;j<3;j++)
    M[j]=0;
  for(k=0;k<n;k++)
    {
      M[2]+=(a[k].r[0]*a[k].u[1]-a[k].r[1]*a[k].u[0])*a[k].m;
      M[0]+=(a[k].r[1]*a[k].u[2]-a[k].r[2]*a[k].u[1])*a[k].m;
      M[1]+=(a[k].r[2]*a[k].u[0]-a[k].r[0]*a[k].u[2])*a[k].m;
    }
   printf("%le %le %le\n",M[0],M[1],M[2]);
*/

  for(j=0;j<3;j++)
    {
      sx=bound[j].length*0.5;
      for(i=0;i<n;i++) 
	a[i].r[j]+=sx;
    }
  /*
	printf("%lf %lf %lf\n%d\n",
bound[0].length,bound[1].length,bound[2].length,n);
  for(i=0;i<n;i++)
    {
	printf("%d %d ",i+1,a[i].c);
      for(j=0;j<3;j++)
	printf("%lf ",a[i].r[j]);
      for(j=0;j<3;j++)
	printf("%lf ",a[i].u[j]);
      printf("\n");
    }
  */

  return;

}



main()
{ int ares;
/*  tail=fopen("junk","w");*/
  if((ares=startup())<1){if(!ares)StopAlert(FILE_ALRT);return;}
  else if(ares>MOVIE){ares-=MOVIE;text_error_dialog(ares);return;}
  if(ares!=MOVIE)
  ct1=squeeze_table(&p1,&q1,&timea); 
  else
  return;
  schedule_init();
  options_dialog();
  event_loop();
/*printf("%lf %lf\n",totaltime,tballtime);*/
   if(!fs_ok())writetext(text_name);
   else fs_close();
   
  if(is_open){fclose(echo_path);}
  if(m_is_open)closemovie(movie_path);
  close_corr_func();
  close_rms();
/*fclose(tail);*/

}
void event_loop(void)
/* Wait for events from the user, and respond to them.  Exit if the functions
   handling the events set the done variable to true. */
{
  while ((get_time()<maxtime)||coll_type) 
    {
      fs_output();
      /*      int ticks=clock(); */


      coll_type=collision();
      update_table(p1,q1,ct1);
      rms();
      if((corr_2>1.0e1)||(corr_2<1.0e-1))
	if(!coll_type)
	  if(!cleanup())break;

      if(new_density)
	  if(!change_density())break;


      ct1=squeeze_table(&p1,&q1,&timea);
      /*    totaltime+=clock()-ticks;*/

    }
}

int readfile (void)
{
  
  int fErr,ares;
  int nbyte;
  int filetype;
  FILE  *path;
  unsigned char * s;
  char fname[80];
  int nn;

  printf("what is file name ?\n");
  scanf("%s",fname);
  path=fopen(fname,"rb");
  if(!path)return 0;
  text_path=path;return TEXT;
}  


double atom_dist(int i,int j)
{
  double dr,dd=0;
  int k;
  iatom * a1=(iatom *)(a+i);
  iatom * a2=(iatom *)(a+j);
  for(k=0;k<dim;k++)
    {
      dr=a1->r[k]-a2->r[k];
      if(dr+dr>bound[k].length)dr-=bound[k].length;
      if(dr+dr<-bound[k].length)dr+=bound[k].length;
      dd+=dr*dr;
   }
 return sqrt(dd);
}



atom * get_atom(void){return a;}
int get_atom_number(void){return n1;}
int get_dimension(void){return (int)dim;}
dimensions * get_bounds(void){return &bound[0];}
double get_movie_dt(void){return delta3;}
int is_reaction(well_type k){int i=coll[k].react;return i>-1?i:~i;}
int is_bond(well_type k){int i=coll[k].react;return (int)(i<0);} 
int is_internal(well_type k){return coll[k].prev==-1?0:1;} 
double etot(well_type k){return coll[k].etot;}
double get_corr(void){return 1/corr_2;}
int get_ll(void){return ll;}
void set_timeb(double t){timeb=t;}
void advance_timeb(){if(timea>timeb)timeb=(timea+timeb)*0.5;}
double get_timeb(void){return timeb;}
double get_timec(void){return timec;}
char * get_text_name(void){return text_name;}
well_type ** get_ecoll(void){return ecoll;}

void set_press_limit(double P0){press_limit=P0;}
void set_press_coeff(double dv_dp){press_coeff=dv_dp;}
void set_L_limit(double * L)
{
  int i=0;
  int ndim=dim; 
  for(i=0;i<ndim;i++)
    if(L[i]>0)
    {
      if(L[i]!=bound[i].length)
	var_density=1;
      L_limit[i]=L[i];
    }
  press_curr=0;
  n_gap_mes=0;
}
 
void set_max_gap_mes(int n_coll_V){max_gap_mes=n_coll_V;}
void set_max_time(double newmaxtime){maxtime=newmaxtime;}
double get_max_time(void){return maxtime;}
