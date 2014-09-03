#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "bcp.h"
#include "search.h"

atom *a ;
dimensions * bound;
double mytimeb;
int N, coll_count, selected_atom_ID;
FILE *stats,*stats2;
int srate,count,sbyte;
unsigned char buf[12900];
long npart,ncell;
double displ2,max_displ2, tot_disp[3];

void init_stats(void){
  a=get_atom();
  N=get_atom_number();
  stats=fopen("stats.txt","w");
  stats2=fopen("stats2.txt","w");
  mytimeb=get_time();
  bound=get_bounds();
  srate=100;
  count=0;
  sbyte=0;
  npart=0;/*number of particle collisions*/
  ncell=0;/*number of cell crossings*/
  sbyte+=sprintf(&buf[sbyte],"# dt\tdr1\tdr2\tdr\td\tneigh1\tneigh2\n");
  coll_count=0;
  selected_atom_ID=1;/*let's pick atom with atom number selected_atom_ID*/
  max_displ2=1.0*1.0; /*let's update every 0.25Angstrom*/
  tot_disp[0]=tot_disp[1]=tot_disp[2]=0.0;
}

void update_stats(int p,int q, double timea){
  
  count++;
  atom *ap,*aq;
  double dt,dx1,dy1,dz1,dx2,dy2,dz2,dx,dy,dz,vx,vy,vz,x,y,z,ix,iy,iz;
  double dr1,dr2,dr,d,mytimea ;

  /*mytimea=get_time();  printf("mytimea=%lf\n",mytimea);*/
  mytimea=timea;
  if(p>=N||q>=N){ncell++;}  
  else{
    npart++;
    ap=a+p;
    aq=a+q;
  
    dt=mytimea-ap->t;
    dx1=ap->v.x*dt;  dy1=ap->v.y*dt;  dz1=ap->v.z*dt;
    dr1=sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
    dt=mytimea-aq->t;
    dx2=aq->v.x*dt;  dy2=aq->v.y*dt;  dz2=aq->v.z*dt;
    dr2=sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
    dx=dx2-dx1; dy=dy2-dy1; dz=dz2-dz1;
    dr=sqrt(dx*dx+dy*dy+dz*dz);    
    
    x=(ap->r.x-aq->r.x)+(ap->v.x-aq->v.x)*dt;
    y=(ap->r.y-aq->r.y)+(ap->v.y-aq->v.y)*dt;
    z=(ap->r.z-aq->r.z)+(ap->v.z-aq->v.z)*dt;
    ix=ap->i.x.i-aq->i.x.i;
    iy=ap->i.y.i-aq->i.y.i;
    iz=ap->i.z.i-aq->i.z.i;
    if (ix>1)x-=bound[0].length;
    else if (ix<-1)x+=bound[0].length;
    if (iy>1)y-=bound[1].length;
    else if (iy<-1)y+=bound[1].length;
    if (iz>1)z-=bound[2].length;
    else if (iz<-1)z+=bound[2].length;
    d=sqrt(x*x+y*y+z*z);

    sbyte+=sprintf(&buf[sbyte],"%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",
		   dt,dr1,dr2,dr,d,get_np(),get_nq());
    if(count>srate){
      fwrite(&buf[0],1,sbyte,stats);  
      count=sbyte=0;
    }
    
    if(p==selected_atom_ID){
      coll_count++;
      tot_disp[0]+=dx1;tot_disp[1]+=dy1;tot_disp[2]+=dz1;
      displ2=tot_disp[0]*tot_disp[0]+tot_disp[1]*tot_disp[1]+tot_disp[2]*tot_disp[2];
      printf("%d %lf %lf\n",coll_count, dr1, sqrt(displ2));
    }
    else if(q==selected_atom_ID){
      coll_count++; 
      tot_disp[0]+=dx2;tot_disp[1]+=dy2;tot_disp[2]+=dz2;
      displ2=tot_disp[0]*tot_disp[0]+tot_disp[1]*tot_disp[1]+tot_disp[2]*tot_disp[2];
      printf("%d %lf %lf\n",coll_count, dr2, sqrt(displ2));
    }
    if(displ2>max_displ2){
      fprintf(stats2, "%d\n",coll_count);
      tot_disp[0]=tot_disp[1]=tot_disp[2]=0.0;
      coll_count=0;
    }
  }
  mytimeb=mytimea;
}


void close_stats(void){
  fwrite(&buf[0],1,sbyte,stats);
  fclose(stats);
}
