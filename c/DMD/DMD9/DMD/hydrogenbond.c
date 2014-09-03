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
extern well_type * collr;
extern double potential;

int isHBRelated(atom* a1, atom* a2){

  if(a1->hb>0 && a2->hb>0){
    int delg=fabs(a2->gid-a1->gid);
    if(delg<deltaGID)return -1;
    return 1;
  }
  else if((a1->hb>0 && a2->hb_asso) || (a2->hb>0 && a1->hb_asso)){
    return -2;
  }
  else if(a1->hb>0 || a2->hb>0){
    int delg=fabs(a2->gid-a1->gid);
    if(a1->gid<0 || a2->gid<0 || delg<deltaGID)return -1;
    return 2;
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
  int ia11=-1, ia12=-1, ia21=-1, ia22=-1;
  atom *a11=NULL, *a12=NULL, *a21=NULL, *a22=NULL;
  int a11_c, a12_c, a21_c, a22_c;
  int a2_a11_ct;
  int a2_a12_ct;
  int a1_a21_ct;
  int a1_a22_ct;
  a1=a+i1;
  a2=a+i2;
  if(a1->hb==1){
    ia11=a1->hb_a1;
    a11=a + a1->hb_a1;
    a11_c=a11->c;
  }
  else if(a1->hb==2){
    ia11=a1->hb_a1;
    a11=a + a1->hb_a1;
    a11_c=a11->c;

    ia12=a1->hb_a2; 
    a12=a + a1->hb_a2; 
    a12_c=a12->c;
  }
  if(a2->hb==1){
    ia21=a2->hb_a1;
    a21=a + a2->hb_a1;
    a21_c=a21->c;
  }
  else if(a2->hb==2){
    ia21=a2->hb_a1;
    a21=a + a2->hb_a1;
    a21_c=a21->c;

    ia22=a2->hb_a2; 
    a22=a + a2->hb_a2; 
    a22_c=a22->c;
  }
  
  old1=a1->c;
  old2=a2->c;
  m1=a1->m;
  m2=a2->m;

  //printf("%ld %ld\n", i1, i2);
  if(!revers)
    {
      /*printf("Try to form the hbonds\n");*/
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

      if(a11){
	a2_a11_ct = checkAssociatedBond(a2, new2, a11, a11_c);
	if(a2_a11_ct==-1)return -1;
      }
      if(a12){
	a2_a12_ct = checkAssociatedBond(a2, new2, a12, a12_c); 
	if(a2_a12_ct==-1)return -1;
      }
      
      if(a21){
	a1_a21_ct = checkAssociatedBond(a1, new1, a21, a21_c);
	if(a1_a21_ct==-1)return -1;
      }
      if(a22){
	a1_a22_ct = checkAssociatedBond(a1, new1, a22, a22_c); 
	if(a1_a22_ct==-1)return -1;
      }
    }
  else
    {
      /*printf("Try to break the hbonds\n");*/
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

      if(a11){
	a2_a11_ct = checkDessociatedBond(a2, new2, a11, a11_c);
	if(a2_a11_ct==-1)return -1;
      }
      if(a12){
	a2_a12_ct = checkDessociatedBond(a2, new2, a12, a12_c); 
	if(a2_a12_ct==-1)return -1;
      }
      
      if(a21){
	a1_a21_ct = checkDessociatedBond(a1, new1, a21, a21_c);
	if(a1_a21_ct==-1)return -1;
      }
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
	      if(collp[ip]<0){
		printf("here %ld %ld\n", i1, ip);
		exit(2);
	      }
	      new_pot+=coll[collp[ip]].etot;
	      //printf("old_pot: %lf new_pot: %lf\n", old_pot, new_pot);
	      //printf("%ld %ld %ld %ld %lf %ld %lf\n",i1+1, ip+1, cp[ip], coll[cp[ip]].etot, collp[ip], coll[collp[ip]].etot);
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
	      if(collq[iq]<0){
		printf("here2\n");
		exit(2);
	      }
	      new_pot+=coll[collq[iq]].etot;
	      //printf("old_pot: %lf new_pot: %lf\n", old_pot, new_pot);
	      //printf("%ld %ld %ld %ld %lf %ld %lf\n",i2+1, iq+1, cq[iq], coll[cq[iq]].etot, collq[iq], coll[collq[iq]].etot);
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
	/**///break_hb(i1,i2);
	if(ia21>=0) breakBond(i1,ia21);
	if(ia22>=0) breakBond(i1, ia22);
	if(ia11>=0) breakBond(i2,ia11);
	if(ia12>=0) breakBond(i2, ia12);
	//printf("%ld %ld break HB\n", i1, i2);
      }
      else if(react[rtype].bond){
	setBond(i1,i2);
	/**///set_hb(i1,i2);
	if(ia21>=0) setBond(i1, ia21);
	if(ia22>=0) setBond(i1, ia22);
	if(ia11>=0) setBond(i2, ia11);
	if(ia12>=0) setBond(i2, ia12);
	//printf("%ld %ld form HB\n", i1, i2);
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

/*(A[B]+C-->A'+B'+C):: atom A or atom B collide with C and change their
  type to A' or B' in the mean time B/A also change their type to B'/A'. 
  The reaction condition have similar citeria as process_hb:
  A's partners---C2, C3
  B's partbers---C,  C1
  checking A-C1, C2-B, C3-B, and A-B distances to see wheather the 
  reaction/[reverse] happens*/
int process_3hb(atom* a, int i1, int i2, int ct1, double sc, double x, double y, double z,
		int ct, int rtype, int revers){
  int A, B, C, C1=-1, C2=-1, C3=-1;
  int i, iA, iB, bond, ip, iq;
  int A_B, A_C1, C2_B, C3_B;
  int nA, nB, nC;
  int * aA=NULL;
  int * aB=NULL;
  int * aC=NULL;
  int * cA=NULL;
  int * cB=NULL;
  int * cC=NULL;
  double sign=1.0;

  double old_pot=coll[ct1].etot;
  double du,duc,ed,di;
  double abA,abC,vx,vy,vz,ab;
  double new_pot=0;
  double mA, mC;
  int ct_new,next,prev,rct;
  int oldA,oldB,oldC,newA,newB=-1,newC; 
  /*find A, B, C repectively*/
  if(a[i1].hb_asso){ 
    A = i2;
    C = i1;
    aC = get_atomp();
    cC = get_collp();
    nC = get_np();

    aA = get_atomq();
    cA = get_collq();
    nA = get_nq();
    sign=-1.0;
    
  }
  else{
    C=i2;
    A=i1;
    aA = get_atomp();
    cA = get_collp();
    nA = get_np();

    aC = get_atomq();
    cC = get_collq();
    nC = get_nq();
  }
  /*and also the partners-C1,C2,C3*/
  if(a[C].hb_p1>=0 &&  a[a[C].hb_p1].origc!=a[A].origc){
    B=a[C].hb_p1;
    if(a[B].hb==0){
      printf("fatal error\n");
      exit(2);
    }
    else if(a[B].hb==2){
      if(C==a[B].hb_a1) C1 = a[B].hb_a2;
      else C1=a[B].hb_a1;
    }
  }
  else if(a[C].hb_p2>=0 && a[a[C].hb_p2].origc!=a[A].origc){
    B=a[C].hb_p2;
    if(a[B].hb==0){
      printf("fatal error\n");
      exit(2);
    }
    else if(a[B].hb==2){
      if(C==a[B].hb_a1) C1 = a[B].hb_a2;
      else C1=a[B].hb_a1;
    }
  }
  else{
    printf("fatal error\n");
    exit(2);
  }
  C2 = a[A].hb_a1;
  C3 = a[A].hb_a2;
  /*checking the seperations*/
  if(fabs(a[A].gid-a[B].gid)<deltaGID) return -1;
  
  oldA = a[A].c; oldB=a[B].c; oldC=a[C].c;
  moveatom(a+B);
  /*find out the new and old types*/
  if(!revers){
    if(oldA==react[rtype].old1){
      newA=react[rtype].new1;
      newC=react[rtype].new2;
    }
    else{
      newA=react[rtype].new2;
      newC=react[rtype].new1;
    }
    //printf("hi %ld %ld\n", oldA, oldB);
    next=ecoll[oldA][oldB];
    while(next>-1){
      rct=coll[next].react;
      if(rct){
	if(oldA==react[rct].old1&&
	   newA==react[rct].new1&&
	   oldB==react[rct].old2){
	  newB=react[rct].new2;
	  break;
	}
	else if(oldA==react[rct].old2&&
		newA==react[rct].new2&&
		oldB==react[rct].old1){
	  newB=react[rct].new1;
	  break;
	}
      }
      next=coll[next].next;
    }
    if(newB<0)return -1;
    //printf("form:");
    ct_new = react[rtype].in;
    /*A-B*/
    A_B = checkAssociatedBond(a+A, newA, a+B, newB);
    if(A_B==-1) return -1;
    /*A-C1*/
    if(C1>-1){
      A_C1 = checkAssociatedBond(a+A, newA, a+C1, a[C1].c);
      if(A_C1==-1) return -1;
    }
    /*C2-B*/
    if(C2>-1){
     C2_B = checkAssociatedBond(a+B, newB, a+C2, a[C2].c);
     if(C2_B==-1) return -1;
   }
   /*C3-B*/
   if(C3>-1){
     C3_B = checkAssociatedBond(a+B, newB, a+C3, a[C3].c);
     if(C3_B==-1) return -1;
   }
  }
  else{
    if(oldA==react[rtype].new1){
      newA = react[rtype].old1;
      newC = react[rtype].old2;
    }
    else{
      newA = react[rtype].old2;
      newC = react[rtype].old1;
    }
    prev = cA[B];
    while(prev>-1){
      prev=coll[prev].prev;
      rct=~coll[prev].react;
      if(rct){
	if(oldA==react[rct].new1&&
	   newA==react[rct].old1&&
	   oldB==react[rct].new2){
	  newB=react[rct].old2;
	  break;
	}
	else if(oldA==react[rct].new2&&
		newA==react[rct].old2&&
		oldB==react[rct].new1){
	  newB=react[rct].old1;
	  break;
	}
      }
    }
    if(newB<0)return -1;
    //printf("%ld %ld %ld %ld %ld %ld %ld %ld %ld \n", A, oldA, newA, B, oldB, newB, C, oldC, newC);
    //printf("break:");
    ct_new = react[rtype].out;
    /*A-B*/
    //printf("%ld %lf %lf\n", cA[B], dist(a[A].r, a[B].r), coll[icoll[7][8]].dd);
    A_B = checkDessociatedBond(a+A, newA, a+B, newB);
    if(A_B==-1) return -1;
    /*A-C1*/
    if(C1>-1){
      A_C1 = checkDessociatedBond(a+A, newA, a+C1, a[C1].c);
      if(A_C1==-1) return -1;
    }
    /*C2-B*/
    if(C2>-1){
      C2_B = checkDessociatedBond(a+B, newB, a+C2, a[C2].c);
      if(C2_B==-1) return -1;
    }
    /*C3-B*/
    if(C3>-1){
      C3_B = checkDessociatedBond(a+B, newB, a+C3, a[C3].c);
      if(C3_B==-1) return -1;
    }
  }
  if(oldC!=newC){
    printf("fatal error\n");
    exit(2);
  }
  
  new_pot += coll[ct_new].etot;
  
  aB = get_atomr();
  cB = get_collr();
  nB=collect_neighbour(B, cB, aB);
  
  if(newA!=oldA){
    a[A].c = newA;
    for(i=0; i<nA; i++){
      iA=aA[i];
      if(iA==B){
	old_pot+=coll[cA[iA]].etot;
	collp[iA]=A_B;
	new_pot+=coll[A_B].etot;
      }
      else if(iA==C1){
	old_pot+=coll[cA[iA]].etot;
	collp[iA]=A_C1;
	new_pot+=coll[A_C1].etot;
      }
      else if(iA!=C){
	old_pot+=coll[cA[iA]].etot;
	moveatom(a+iA);
	bond=is_bond(cA[iA]);
	collp[iA]=after_type(A, iA, &bond, cA[iA]);
	//printf("%ld %ld %ld %ld\n", A, iA, cA[iA], collp[iA]);
	if(collp[iA]<0){
	  printf("here A %ld %ld %lf %lf \n", A, iA, a[A].t, a[iA].t);
	  a[A].c=oldA;
	  clear_r(nB);
	  exit(2);
	}
	new_pot+=coll[collp[iA]].etot;
	if(bond)collp[iA]=~(collp[iA]);
      }
    }
  }

  if(newB!=oldB){
    moveatom(a+B);
    a[B].c=newB;
    for(i=0; i<nB; i++){
      iB=aB[i];
      //printf(" %ld %ld %ld \n", B, iB, cB[iB]);
      if(iB==C2){
	old_pot+=coll[cB[iB]].etot;
	collr[iB]=C2_B;
	new_pot+=coll[C2_B].etot;
      }
      else if(iB==C3){
	old_pot+=coll[cB[iB]].etot;
	collr[iB]=C3_B;
	new_pot+=coll[C3_B].etot;
      }
      else{
	old_pot+=coll[cB[iB]].etot;
	moveatom(a+iB);
	bond=is_bond(cB[iB]);
	collr[iB]=after_type(B, iB, &bond, cB[iB]);
	if(collr[iB]<0){
	  printf("here B %ld %ld %ld %ld %ld %ld\n", 
		 B, iB, oldB, newB, a[iB].c, is_bond(cB[iB]));
	  a[A].c=oldA;
	  a[B].c=oldB;
	  clear_r(nB);
	  exit(2);
	}
	new_pot+=coll[collr[iB]].etot;
	if(bond) collr[iB]=~(collr[iB]);
      }
    }
    /*C-B*/
    moveatom(a+C);
    bond=is_bond(cC[B]);
    collq[B]=after_type(C, B, &bond, cC[B]);
    if(collq[C]<0){
      printf("here BC %ld %ld\n", B, iB);
      a[A].c=oldA;
      a[B].c=oldB;
      clear_r(nB);
      exit(2);
    }
    if(bond) collq[B]=~(collq[B]);
  }
  du=new_pot-old_pot;
  //printf("du=%lf\n", du);

  if(!react[rtype].bond)du+=react[rtype].eo;
  duc=du*corr_2;
  mA=a[A].m; mC=a[C].m;
  ed=2*duc/(mA*mC*coll[ct].dm);
  di=1.0+ed/(sc*sc);
  if(di<=0){
    /*	when ed is large negative, it is
	unsuccessfull attempt to escape: 
	reaction do not happen */   
    ab=-2.0*sc*coll[ct].dm;
    a[A].c=oldA;
    a[B].c=oldB;
    clear_r(nB);
    if(!revers){
      return -1;           
    }
    ct_new=ct1;
  }
  else{
    ab=sc*coll[ct].dm*(sqrt(di)-1.0);
    vvm+=duc;
    potential+=du;
    if((react[rtype].old1!=react[rtype].new1)||(react[rtype].old2!=react[rtype].new2))
      setNewTypes(1);
    if(revers){
      breakBond(A,B);
      /**///break_hb(A,B);
      breakBond(A,C);
      if(C1>=0) breakBond(A,C1);
      if(C2>=0) breakBond(B,C2);
      if(C3>=0) breakBond(B,C3);
      //printf("Successful breaking!--3HB %ld %ld %ld @ %lf\n", A, B, C, timea);
    }
    else if(react[rtype].bond){
      setBond(A,B);
      /**///set_hb(A,B);
      setBond(A,C);
      if(C1>=0) setBond(A,C1);
      if(C2>=0) setBond(B,C2);
      if(C3>=0) setBond(B,C3);
      //printf("Successful hbonding!--3HB %ld %ld %ld @ %lf\n", A, B, C, timea);
    }
    if(newA!=oldA){
      for(i=0; i<nA; i++){
	ip=aA[i];
	if(ip!=C)
	  if(collp[ip]<0){
	    breakBond(A, ip);
	    cA[ip]=~collp[ip];
	  }
	  else
	    cA[ip]=collp[ip];
      }
      cA[C]=ct_new;
    }
    if(newB!=oldB){
      for(i=0; i<nB; i++){
	iq = aB[i];
	//printf("%ld\n", iq);
	if(iq!=C && iq!=A){
	  if(collr[iq]<0){
	    breakBond(iq,B);
	    cB[iq]=~collr[iq];
	  }
	  else
	    cB[iq]=collr[iq];
	}
      }
      /*C_B*/
      if(collq[B]<0){
	breakBond(C,B);
	cC[B]=~collq[B];
      }
      else{
	cC[B]=collq[B];
      }
      set_nr(B, nB);
    }
  }

  if(sign==-1){
    abA=-ab*mC;
    abC=ab*mA;
  }
  else{
    abA=ab*mC;
    abC=-ab*mA;
  }
  virial+=abA*mA*coll[ct].dd;
  a[A].v.x+=x*abA;
  a[A].v.y+=y*abA;
  a[A].v.z+=z*abA;
  a[C].v.x+=x*abC;
  a[C].v.y+=y*abC;
  a[C].v.z+=z*abC;

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
  //if(nct<0)printf("%ld\n", nct);
  if(nct<0)return nct;
  if(dd>coll[nct].dd || coll[nct].prev!=-1)return -1;
  nct=coll[nct].next;
  while(nct>-1){
    if(dd>coll[nct].dd)return nct;
    nct = coll[nct].next;
  }
  return -1;
}
  
  
int checkDessociatedBond(atom* a1, int newType, atom* asso, int assoType){
  int nct;
  double rx, ry, rz;
  int irx, iry, irz;
  double dd;
  moveatom(asso);
  rx=a1->r.x-asso->r.x;
  ry=a1->r.y-asso->r.y;
  rz=a1->r.z-asso->r.z;
  irx=a1->i.x.i-asso->i.x.i;
  iry=a1->i.y.i-asso->i.y.i;
  irz=a1->i.z.i-asso->i.z.i;
  if (irx>1)rx-=bound[0].length;
  else if (irx<-1)rx+=bound[0].length;
  if (iry>1)ry-=bound[1].length;
  else if (iry<-1)ry+=bound[1].length;
  if (irz>1)rz-=bound[2].length;
  else if (irz<-1)rz+=bound[2].length;   
  dd=rx*rx+ry*ry+rz*rz;
  nct = ecoll[newType][assoType];
  while(nct>-1){
    //if(timea>4198)printf("Desso: %ld %ld %lf, %lf, %ld\n", newType, assoType, sqrt(dd), sqrt(coll[nct].dd), nct);
    if(dd>=coll[nct].dd)return nct;
    nct=coll[nct].next;
  }
  return -1;
}
