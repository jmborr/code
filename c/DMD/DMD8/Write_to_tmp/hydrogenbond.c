#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "bcp.h"
#include "make_system.h"
#include "search.h"
#include "bonds.h"
#include "controls.h"
#include "hydrogenbond.h"

extern deltaGID;
extern dimensions bound[3];
extern CollisionData* coll;
extern well_type ** ecoll;
extern well_type ** icoll;
extern ReactionData *react;
extern double vvm,corr,corr_2,timea,ts,virial;
extern well_type * collp;
extern well_type * collq;
extern double potential;

int isHBRelated(atom* a1, atom* a2){

  if(a1->hb && a2->hb){
    int delg=fabs(a2->gid-a1->gid);
    if(delg<deltaGID)return 0;
    return 1;
  }
  return 0;
}

int process_hb(atom* a, int i1, int i2, int ct1, double sc, double x, double y, double z, 
	       int ct, int rtype, int revers)
{
  int out,i,j,ix,iy,iz;
  double ab1,ab2,vx,vy,vz,ab,di,ed;
  atom *a1,*a2;
  double old_pot=coll[ct1].etot;
  double du,duc,new_pot=0;
  int ct_new;
  int np=get_np();  
  int nq=get_nq();  
  int * ap=get_atomp();  
  int * aq=get_atomq();  
  int * cp=get_collp();  
  int * cq=get_collq();
  int bond,iq,ip;
  int old1,old2,new1,new2;  
  double m1,m2;
  int ia11, ia12=-1, ia21, ia22=-1;
  atom *a11, *a12=NULL, *a21, *a22=NULL;
  int a11_c, a12_c, a21_c, a22_c;
  int a2_a11_ct;
  int a2_a12_ct;
  int a1_a21_ct;
  int a1_a22_ct;
  a1=a+i1;
  a2=a+i2;

  ia11=a1->hb_a1;
  a11=a + a1->hb_a1;
  a11_c=a11->c;
  if(a1->hb==2){ia12=a1->hb_a2; a12=a + a1->hb_a2; a12_c=a12->c;}
  ia21=a2->hb_a1;
  a21=a + a2->hb_a1;
  a21_c=a21->c;
  if(a2->hb==2){ia22=a2->hb_a2; a22=a + a2->hb_a2; a22_c=a22->c;}
  
  old1=a1->c;
  old2=a2->c;
  m1=a1->m;
  m2=a2->m;
  if(!revers)
    {
      if(old1==react[rtype].old1)
	{
	  new1=react[rtype].new1;
	  new2=react[rtype].new2;
	}
      else 
	{
	  new1=react[rtype].new2;
	  new2=react[rtype].new1;
	}

      a2_a11_ct = checkAssociatedBond(a2, new2, a11, a11_c);
      if(a2_a11_ct==-1)return -1;
      if(a12){
	a2_a12_ct = checkAssociatedBond(a2, new2, a12, a12_c); 
	if(a2_a12_ct==-1)return -1;
      }
      
      a1_a21_ct = checkAssociatedBond(a1, new1, a21, a21_c);
      if(a1_a21_ct==-1)return -1;
      if(a22){
	a1_a22_ct = checkAssociatedBond(a1, new1, a22, a22_c); 
	if(a1_a22_ct==-1)return -1;
      }
    }
  else
    {
      if(old1==react[rtype].new1)
	{
	  new1=react[rtype].old1;
	  new2=react[rtype].old2;
	}
      else 
	{
	  new1=react[rtype].old2;
	  new2=react[rtype].old1;
	}

      a2_a11_ct = checkDessociatedBond(a2, new2, a11, a11_c);
      if(a2_a11_ct==-1)return -1;
      if(a12){
	a2_a12_ct = checkDessociatedBond(a2, new2, a12, a12_c); 
	if(a2_a12_ct==-1)return -1;
      }
      
      a1_a21_ct = checkDessociatedBond(a1, new1, a21, a21_c);
      if(a1_a21_ct==-1)return -1;
      if(a22){
	a1_a22_ct = checkDessociatedBond(a1, new1, a22, a22_c); 
	if(a1_a22_ct==-1)return -1;
      }
      
    }
  
  if(revers)
    ct_new=react[rtype].out;
  else
    ct_new=react[rtype].in;
  
  new_pot+=coll[ct_new].etot;
  
  if(old1!=new1)
    {
      a1->c=new1;
      for(i=0;i<np;i++)
	{
	  ip=ap[i];
	  if(ip==ia21){
	    old_pot+=coll[cp[ip]].etot;
	    collp[ip]=a1_a21_ct;
	    new_pot+=coll[collp[ip]].etot;
	  }
	  else if(ip==ia22){
	    old_pot+=coll[cp[ip]].etot;
	    collp[ip]=a1_a22_ct;
	    new_pot+=coll[collp[ip]].etot;
	  }
	  else if(ip!=i2)
	    {
	      old_pot+=coll[cp[ip]].etot;
	      moveatom(a+ip);
	      bond=is_bond(cp[ip]);
	      collp[ip]=after_type(i1,ip,&bond,cp[ip]);
	      if(collp[ip]<0)return -1;
	      new_pot+=coll[collp[ip]].etot;
	      if(bond)collp[ip]=~(collp[ip]);
	    }
	}
    }
  
  if(old2!=new2)
    {
      a2->c=new2;
      for(i=0;i<nq;i++)
	{
	  iq=aq[i];
	  if(iq==ia11){
	    old_pot+=coll[cq[iq]].etot;
	    collq[iq]=a2_a11_ct;
	    new_pot+=coll[collq[iq]].etot;
	  }
	  else if(iq==ia12){
	    old_pot+=coll[cq[iq]].etot;
	    collq[iq]=a2_a12_ct;
	    new_pot+=coll[collq[iq]].etot;
	  }
	  else if(iq!=i1)
	    {
	      old_pot+=coll[cq[iq]].etot;
	      moveatom(a+iq);
	      bond=is_bond(cq[iq]);
	      collq[iq]=after_type(i2,iq,&bond,cq[iq]);
	      if(collq[iq]<0)return -1;
	      new_pot+=coll[collq[iq]].etot;
	      /* we remember that bonds was broken storing negatives in collq */ 
	      if(bond)collq[iq]=~(collq[iq]);
	    }
	}
    }
  
  du=new_pot-old_pot;
  if(!react[rtype].bond)du+=react[rtype].eo;
  duc=du*corr_2;
  ed=2*duc/(m1*m2*coll[ct].dm);
  di=1.0+ed/(sc*sc);
  if(di<=0)
    {
      
      /*	when ed is large negative, it is
		unsuccessfull attempt to escape: 
		reaction do not happen */   
      ab=-2.0*sc*coll[ct].dm;
      a1->c=old1;
      a2->c=old2;
      if(!revers)return -1;           
      ct_new=ct1;
    }
  else
    {
      ab=sc*coll[ct].dm*(sqrt(di)-1.0);
      vvm+=duc;
      potential+=du;
      if((react[rtype].old1!=react[rtype].new1)||(react[rtype].old2!=react[rtype].new2))
	setNewTypes(1);
      if(revers){
	breakBond(i1,i2);
	breakBond(i1,ia21);
	if(ia22>=0) breakBond(i1, ia22);
	breakBond(i2,ia11);
	if(ia12>=0) breakBond(i2, ia12);
      }
      else if(react[rtype].bond){
	setBond(i1,i2);
	setBond(i1, ia21);
	if(ia22>=0) setBond(i1, ia22);
	setBond(i2, ia11);
	if(ia12>=0) setBond(i2, ia12);
	/*printf("Successful hbonding!\n");*/
	  
      }
      if(new2!=old2)
	for(i=0;i<nq;i++)
	  {
	    iq=aq[i];
	    if(iq!=i1)
	      {
		if(collq[iq]<0)
		  {
		    breakBond(i2,iq);
		    cq[iq]=~collq[iq];
		  }
		else
		  cq[iq]=collq[iq];
	      }
	  }
      if(new1!=old1)
	for(i=0;i<np;i++)
	  {
	    ip=ap[i];
	    if(ip!=i2)
	      {
		if(collp[ip]<0)
		  {
		    breakBond(i1,ip);
		    cp[ip]=~collp[ip];
		  }
		else
		  cp[ip]=collp[ip];
	      }
	  }
      cp[i2]=ct_new;
      
    }
  
  ab1=ab*m2;
  ab2=-ab*m1;
  virial+=ab1*m1*coll[ct].dd;
  a1->v.x+=x*ab1;
  a1->v.y+=y*ab1;
  a1->v.z+=z*ab1;
  a2->v.x+=x*ab2;
  a2->v.y+=y*ab2;
  a2->v.z+=z*ab2;
  
  return ct_new; 
  
}

int checkAssociatedBond(atom* a, int newType, atom* asso, int assoType){
  int nct;
  double rx, ry, rz;
  int irx, iry, irz;
  double dd;
  moveatom(asso);
  rx=a->r.x-asso->r.x;
  ry=a->r.y-asso->r.y;
  rz=a->r.z-asso->r.z;
  irx=a->i.x.i-asso->i.x.i;
  iry=a->i.y.i-asso->i.y.i;
  irz=a->i.z.i-asso->i.z.i;
  if (irx>1)rx-=bound[0].length;
  else if (irx<-1)rx+=bound[0].length;
  if (iry>1)ry-=bound[1].length;
  else if (iry<-1)ry+=bound[1].length;
  if (irz>1)rz-=bound[2].length;
  else if (irz<-1)rz+=bound[2].length;   
  dd=rx*rx+ry*ry+rz*rz;
  nct=icoll[newType][assoType];
  if(nct<0)return nct;
  if(dd>coll[nct].dd || coll[nct].prev!=-1)return -1;
  nct=coll[nct].next;
  while(nct>-1){
    if(dd>coll[nct].dd)return nct;
    nct = coll[nct].next;
  }
  return -1;
}
  
  
int checkDessociatedBond(atom* a, int newType, atom* asso, int assoType){
  int nct;
  double rx, ry, rz;
  int irx, iry, irz;
  double dd;
  moveatom(asso);
  rx=a->r.x-asso->r.x;
  ry=a->r.y-asso->r.y;
  rz=a->r.z-asso->r.z;
  irx=a->i.x.i-asso->i.x.i;
  iry=a->i.y.i-asso->i.y.i;
  irz=a->i.z.i-asso->i.z.i;
  if (irx>1)rx-=bound[0].length;
  else if (irx<-1)rx+=bound[0].length;
  if (iry>1)ry-=bound[1].length;
  else if (iry<-1)ry+=bound[1].length;
  if (irz>1)rz-=bound[2].length;
  else if (irz<-1)rz+=bound[2].length;   
  dd=rx*rx+ry*ry+rz*rz;
  nct = ecoll[newType][assoType];
  while(nct>-1){
    if(dd>=coll[nct].dd)return nct;
    nct=coll[nct].next;
  }
  return -1;
}
