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

extern deltaGID;//minimal amino acid separation to allow a HB 
extern dimensions bound[3];
extern CollisionData* coll;
extern well_type ** ecoll;
extern well_type ** icoll;
extern ReactionData *react;
extern double vvm,corr,corr_2,timea,ts,virial;
extern well_type * collp;
extern well_type * collq;
extern double potential;

int isHBRelated(atom* a1, atom* a2){  /*check whether a1 and a2 can form a HB*/
  if(a1->hb && a2->hb){             /*atoms can react to form a hydrogen bond*/
    int delg=fabs(a2->gid-a1->gid);/*delg,index diff. between two amino acids*/
    if(delg<deltaGID)return 0;/*deltaGID minimum separation to allow to amino*/
    return 1;                                 /*acids to form a hydrogen bond*/
  }
  return 0;
}

int process_hb(atom* a, int i1, int i2, int ct1, double sc, double x, double y,
	       double z, int ct, int rtype, int revers )  /*270 lines of code*/
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

  /*find association atoms for a1*/
  ia11=a1->hb_a1;        /*atom index of previous association atom of atom a1*/
  a11=a + a1->hb_a1;                   /*previous association atom of atom a1*/
  a11_c=a11->c;        /*current type of previous association atom of atom a1*/
  if(a1->hb==2){ /*if the atom needs of two helper association atoms, then...*/
    ia12=a1->hb_a2;          /*atom index of next association atom of atom a1*/
    a12=a + a1->hb_a2;                     /*next association atom of atom a1*/
    a12_c=a12->c;          /*current type of next association atom of atom a1*/
  }

  /*find association atoms for a2*/
  ia21=a2->hb_a1;
  a21=a + a2->hb_a1;
  a21_c=a21->c;                        /*previous association atom of atom a2*/
  if(a2->hb==2){
    ia22=a2->hb_a2;
    a22=a + a2->hb_a2;                    /*next  association atom of atom a2*/
    a22_c=a22->c;
  }
  
  old1=a1->c;                                /*backup the current types of a1*/
  old2=a2->c;
  m1=a1->m;                                                      /*mass of a1*/
  m2=a2->m;

  if(!revers)  /*revers==0 if particles approaching each other. The bond will*/
    {                                   /*be created: old1+old2 --> new1+new2*/
      if(old1==react[rtype].old1)        /*old1 is type of "p" before bonding*/
	{
	  new1=react[rtype].new1;         /*new1 is type of "p" after bonding*/
	  new2=react[rtype].new2;         /*new2 is type of "q" after bonding*/
	}
      else    /*If type of p listed as react[rtype].old2,then new1 still will*/
	{      /*be type of p after bonding, but refers to react[rtype].new2 */
	  new1=react[rtype].new2;         /*new1 is type of "p" after bonding*/
	  new2=react[rtype].new1;         /*new2 is type of "q" after bonding*/
	}
      /*Now check whether auxiliary bond a1-a11 is feasible.If so,return coll*/
      a2_a11_ct = checkAssociatedBond(a2, new2, a11, a11_c);  /*type of bond*/
      if(a2_a11_ct==-1)return -1;         /*auxiliary bond a11-a2 is feasible*/
      if(a12){              /*Check whether auxiliary bond a12-a2 is feasible*/
	a2_a12_ct = checkAssociatedBond(a2, new2, a12, a12_c); 
	if(a2_a12_ct==-1)return -1;   /*a12-a2 not possible, thus no reaction*/
      }
      /*Now check whether auxiliary bond a1-a21 is feasible.If so,return coll*/
      a1_a21_ct = checkAssociatedBond(a1, new1, a21, a21_c);   /*type of bond*/
      if(a1_a21_ct==-1)return -1;
      if(a22){              /*Check whether auxiliary bond a1-a22 is feasible*/
	a1_a22_ct = checkAssociatedBond(a1, new1, a22, a22_c); 
	if(a1_a22_ct==-1)return -1;   /*a12-a2 not possible, thus no reaction*/
      }
    }
  else /*revers==1 if particles moving away from each other. The bond will be*/
    {                                      /* broken: old1+old2 <-- new1+new2*/
      if(old1==react[rtype].new1)            /*atoms when in the bonded state*/
	{
	  new1=react[rtype].old1;   /*new1, type of p after breaking the bond*/
	  new2=react[rtype].old2;   /*new2, type of q after breaking the bond*/
	}
      else 
	{
	  new1=react[rtype].old2;   /*new1, type of p after breaking the bond*/
	  new2=react[rtype].old1;   /*new2, type of q after breaking the bond*/
	}
      
      a2_a11_ct = checkDessociatedBond(a2, new2, a11, a11_c); /*Check whether*/
      if(a2_a11_ct==-1)return -1;        /*auxiliary bond a11-a2 is breakable*/
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
  
  if(revers)                                        /*The bond will be broken*/
    ct_new=react[rtype].out;  /*react[rtype].out, coll. type of unbound pair */
  else                                             /*The bond will be created*/
    ct_new=react[rtype].in;    /*react[rtype].out, coll. type of bonded pair */
  
  new_pot+=coll[ct_new].etot;                         /*Energy after reaction*/
  
  if(old1!=new1)           /*The atomic type of "p" has changed upon reaction*/
    {
      a1->c=new1;                               /*Update the type of atom "p"*/
      for(i=0;i<np;i++)  /*Go through atoms with previous scheduled collision*/
	{
	  ip=ap[i];               /*Atom number of one the "neighbors" of "p"*/
	  if(ip==ia21){        /*If neighbor is first association atom of "q"*/
	    old_pot+=coll[cp[ip]].etot;       /*update energy before reaction*/
	    collp[ip]=a1_a21_ct;  /*store a1-a21 coll. type in collp[]. Mind */
	     /*collp is not search.collp,but an extra array declared in bcp.c*/
	    new_pot+=coll[collp[ip]].etot;       /*total energy after bonding*/
	  }
	  else if(ip==ia22){  /*If neighbor is second association atom of "q"*/
	    old_pot+=coll[cp[ip]].etot;
	    collp[ip]=a1_a22_ct;
	    new_pot+=coll[collp[ip]].etot;
	  }
	  else if(ip!=i2)                            /*If neighbor is not "q"*/
	    {   /*cp[ip] is collision type of previously scheduled p-ip coll.*/
	      old_pot+=coll[cp[ip]].etot;
	      moveatom(a+ip);       /*Update position of "ip" to system time */
	      bond=is_bond(cp[ip]);          /*if "p" and "ip" bonded, bond=1*/
	      collp[ip]=after_type(i1,ip,&bond,cp[ip]);       /*new coll.type*/
	      if(collp[ip]<0)return -1;
	      new_pot+=coll[collp[ip]].etot;
	      if(bond)collp[ip]=~(collp[ip]);   /*If "p" and "ip" were bonded*/
	    }   /*before "p" underwent reaction,but then the bond was broken */
	}   /*in after_type(..),then bond==1 and we take the bitwise         */
    }  /*complement of collp[ip], which is -collp[ip]-1, a negative number   */
  
  if(old2!=new2)           /*The atomic type of "p" has changed upon reaction*/
    {
      a2->c=new2;                               /*Update the type of atom "p"*/
      for(i=0;i<nq;i++)/*Go through atoms with previous schedld. coll. with q*/
	{
	  iq=aq[i];               /*Atom number of one the "neighbors" of "q"*/
	  if(iq==ia11){        /*If neighbor is first association atom of "p"*/
	    old_pot+=coll[cq[iq]].etot;       /*update energy before reaction*/
	    collq[iq]=a2_a11_ct;  /*store a1-a21 coll. type in search.collq[]*/
	    new_pot+=coll[collq[iq]].etot;       /*total energy after bonding*/
	  }
	  else if(iq==ia12){  /*If neighbor is second association atom of "q"*/
	    old_pot+=coll[cq[iq]].etot;
	    collq[iq]=a2_a12_ct;
	    new_pot+=coll[collq[iq]].etot;
	  }
	  else if(iq!=i1)                            /*If neighbor is not "p"*/
	    {   /*cq[ip] is collision type of previously scheduled q-iq coll.*/
	      old_pot+=coll[cq[iq]].etot;
	      moveatom(a+iq);
	      bond=is_bond(cq[iq]);          /*if "q" and "iq" bonded, bond=1*/
	      collq[iq]=after_type(i2,iq,&bond,cq[iq]);       /*new coll.type*/
	      if(collq[iq]<0)return -1;
	      new_pot+=coll[collq[iq]].etot; 
	      if(bond)collq[iq]=~(collq[iq]);  /*If "q" and "iq" were bonded*/
	    }   /*before "q" underwent reaction,but then the bond was broken */
	}   /*in after_type(..),then bond==1 and we take the bitwise         */
    }  /*complement of collq[iq], which is -collp[iq]-1, a negative number   */
  
  du=new_pot-old_pot;              /*Change in total energy upon p-q reaction*/
  if(!react[rtype].bond)du+=react[rtype].eo; /*If bonding,release latent heat*/
  duc=du*corr_2;             /*Adjusted change in energy for the thermal bath*/
  ed=2*duc/(m1*m2*coll[ct].dm);
  di=1.0+ed/(sc*sc);
  if(di<=0)    /*When ed is large negative, unsuccessfull attempt to escape: */
    {	       /*reaction does not happen                                    */
      ab=-2.0*sc*coll[ct].dm;
      a1->c=old1;    /*Return the old types of a1 and a2 to the state before */
      a2->c=old2;    /*We evaluated reaction(..)                             */
      if(!revers)return -1;         /*revers==0 if the bond was to be created*/
      ct_new=ct1; ,   /*Return the old collision type between "i" and "k"    */
    }
  else                        /*di>0, thus reaction is energetically possible*/
    {
      ab=sc*coll[ct].dm*(sqrt(di)-1.0);
      vvm+=duc;
      potential+=du;
      if( (react[rtype].old1!=react[rtype].new1) ||
	  (react[rtype].old2!=react[rtype].new2)    )
	setNewTypes(1);                   /*Set flag "newTypes" to value of 1*/
      if(revers){                             /*The reaction breaks the bond */
	breakBond(i1,i2);                                   /*break p-q* bond*/
	breakBond(i1,ia21);                            /*break auxiliary bond*/
	if(ia22>=0) breakBond(i1, ia22);
	breakBond(i2,ia11);
	if(ia12>=0) breakBond(i2, ia12);
      }
      else if(react[rtype].bond){              /*The reaction produces a bond*/
	setBond(i1,i2);
	setBond(i1, ia21);
	if(ia22>=0) setBond(i1, ia22);
	setBond(i2, ia11);
	if(ia12>=0) setBond(i2, ia12);
      }
      if(new2!=old2)         /*If the type of "q" has changed in the reaction*/
	for(i=0;i<nq;i++)                 /*Go through the "neighbors" of "q"*/
	  {
	    iq=aq[i];                   /*Atom number of one of the neighbors*/
	    if(iq!=i1)                    /*Make sure the neighbor is not "p"*/
	      {
		if(collq[iq]<0)  /*Negative if the bond q-iq broke because of*/
		  {                                        /*the reaction p-q*/
		    breakBond(i2,iq);                      
		    cq[iq]=~collq[iq];  /*Insert the correct (positive) type */
		  }            /*between "q" and "iq" into "search.collp[ip]"*/
		else
		  cq[iq]=collq[iq];     /*Either the bond q-iq was conserved */
	      }    /*after reaction or there was no bond between "q" and "iq"*/
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
      cp[i2]=ct_new;    /*insert the new collisiont type between "p" and "q" */
      
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
}/*Matches int process_hb(..)*/

/*Check whether there is a bond defined between the bonded type of "a"(newType,
  new because "a" is not yet participating in  a hydrogen bond) and the type of
  "asso" (assoType). If so, check interparticle distance smaller that most 
  external potential barrier of the a-asso bond*/
int checkAssociatedBond(atom* a, int newType, atom* asso, int assoType){
  int nct;
  double rx, ry, rz;
  int irx, iry, irz;
  double dd;
  moveatom(asso);

  rx=a->r.x-asso->r.x;               /*relative vector between "a" and "asso"*/
  ry=a->r.y-asso->r.y;
  rz=a->r.z-asso->r.z;

  irx=a->i.x.i-asso->i.x.i;     /*relative cell-vector between "a" and "asso"*/
  iry=a->i.y.i-asso->i.y.i;
  irz=a->i.z.i-asso->i.z.i;

  if (irx>1)rx-=bound[0].length;         /*check periodic boundary conditions*/
  else if (irx<-1)rx+=bound[0].length;
  if (iry>1)ry-=bound[1].length;
  else if (iry<-1)ry+=bound[1].length;
  if (irz>1)rz-=bound[2].length;
  else if (irz<-1)rz+=bound[2].length;   

  dd=rx*rx+ry*ry+rz*rz;              /*Square distance between "a" and "asso"*/
  nct=icoll[newType][assoType];          /*a-asso collission type when bonded*/
  if(nct<0)return nct;             /*No bond possible, return negative number*/
  if(dd>coll[nct].dd || coll[nct].prev!=-1)return -1;/*"a" and "asso" farther*/
  nct=coll[nct].next;        /*than most external pot. barrier of a-asso bond*/
  while(nct>-1){     /*return appropiate potential barrier of the a-asso bond*/
    if(dd>coll[nct].dd)return nct;
    nct = coll[nct].next;
  }
  return -1;
}
  
/*Check if auxiliary bond a-asso can be broken. For this, the distance between
  "a" and "asso" must be bigger than hard-core distance between type of "a" 
  after breaking the hydrogen bond(newType) and the type of "asso"(assoType)*/ 
int checkDessociatedBond(atom* a, int newType, atom* asso, int assoType){
  int nct;
  double rx, ry, rz;
  int irx, iry, irz;
  double dd;
  moveatom(asso);
  rx=a->r.x-asso->r.x;               /*relative vector between "a" and "asso"*/
  ry=a->r.y-asso->r.y;
  rz=a->r.z-asso->r.z;
  irx=a->i.x.i-asso->i.x.i;     /*relative cell-vector between "a" and "asso"*/
  iry=a->i.y.i-asso->i.y.i;
  irz=a->i.z.i-asso->i.z.i;
  if (irx>1)rx-=bound[0].length;         /*check periodic boundary conditions*/
  else if (irx<-1)rx+=bound[0].length;
  if (iry>1)ry-=bound[1].length;
  else if (iry<-1)ry+=bound[1].length;
  if (irz>1)rz-=bound[2].length;
  else if (irz<-1)rz+=bound[2].length;   
  dd=rx*rx+ry*ry+rz*rz;              /*Square distance between "a" and "asso"*/
  nct = ecoll[newType][assoType];       /*a-asso collission type when unbound*/
  while(nct>-1){
    if(dd>=coll[nct].dd)return nct;
    nct=coll[nct].next;
  }
  return -1;
}
