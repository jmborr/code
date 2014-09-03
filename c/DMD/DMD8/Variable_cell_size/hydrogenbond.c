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
#include "cell_size.h"

extern atom* a;
extern tsearch search;
extern deltaGID;
extern dimensions bound[3];
extern CollisionData* coll;
extern well_type ** ecoll;
extern well_type ** icoll;
extern ReactionData *react;
extern double vvm,corr,corr_2,timea,ts,virial;
extern double potential;
extern tree_scheme calendar;
alist ***tlists;
enum list_type *winlist;
int amso;
void set_amso_for_HB(int value){
  amso=value;
  tlists=search.tlists;
  winlist=search.winlist;
}
/*=========================================================================*/
int isHBRelated(atom* a1, atom* a2){
  /*printf("isHBRelated(...)*/
  if(a1->hb && a2->hb){
    int delg=fabs(a2->gid-a1->gid);
    if(delg<deltaGID)return 0;
    return 1;
  }
  return 0;
}
/*=========================================================================*/
/*In process_hbII, we assess whether react. successful. If so,then break the 
  existing bonds or set new bonds. Also, since atoms change their types, we 
  need to update the collision types between "p" and the neighbors of "p", and
  between "q" and the neighbors of "q". If say, "p" and neighbor "ip" were
  bonded before p-q reaction, it may be that they can not be bonded after p-q
  reaction. In that case we may break the bond*/
/*ct1 is collision type of the "next" potential step between i1 and i2. ct is
  the collision type of the reaction. ct can be ct1 or coll[ct1].prev*/
/*IMPORTANT NOTE: We assume that interaction potentials between i1 and
  any other neighbor atom are the same either when free or hydrogen
  bonded (except for i2 and the two associations of i1, of
  course). The same runs for i2. If assumption is not true, then we
  would have to change each collision in addition to the collisions with the
  association atoms. We would also check whether the collision after reaction
  needs to change from {HOT,WARM,COLD} to other list type. We would
  also have to check for new neighbors and/or discard old neighbors,
  since the surface shell orders may also have changed.*/
/*IMPORTANT NOTE: assumed i1<i2 in the circular sense*/
coll_type process_hbII (atom* a,atom_number i1,atom_number i2,coll_type ct1,
			double sc,double x,double y,double z,coll_type ct,
			coll_type rtype,int revers){/*270 lines*/
  double m1,m2,ab,ab1,ab2,di,ed,du,duc,old_pot,new_pot,t;
  double a1_a21_dd,a1_a22_dd,a2_a11_dd,a2_a12_dd;
  atom *a1,*a2;
  int old1,old2,new1,new2,bond;  
  coll_type ct_new;
  coll_type a1_a21_ct,a1_a22_ct,a2_a11_ct,a2_a12_ct;
  atom_number ia11, ia12, ia21, ia22;
  atom_type a11_c,a12_c,a21_c,a22_c;
  atom *a11,*a12,*a21,*a22;
  alist *a1_a2_coll,*a1_a21_coll,*a1_a22_coll,*a2_a11_coll,*a2_a12_coll;
  enum list_type lt,new_lt,a1_a2_lt,a1_a21_lt,a1_a22_lt,a2_a11_lt,a2_a12_lt;
  shell_order mso,so;
  /*printf("process_hbII(%d,%d)\n",i1,i2);*/

  a1=a+i1;  a2=a+i2;
  /*assotiations of a1*/
  ia11=a1->hb_a1;  a11=a + a1->hb_a1; a11_c=a11->c;
  ia12=a1->hb_a2;  a12=a + a1->hb_a2; a12_c=a12->c;
  /*assotiations of a2*/
  ia21=a2->hb_a1;  a21=a + a2->hb_a1; a21_c=a21->c;
  ia22=a2->hb_a2;  a22=a + a2->hb_a2; a22_c=a22->c;

  old1=a1->c;  old2=a2->c;

  if(!revers){/*no reverse reaction, bond may form*/
    /*determine who is N, who is C*/
    if(old1==react[rtype].old1 && old2==react[rtype].old2){
      new1=react[rtype].new1; new2=react[rtype].new2;
    }
    else if(old1==react[rtype].old2 && old2==react[rtype].old1){
      new1=react[rtype].new2;/*type of p after bonding, but refers to */
      new2=react[rtype].new1;/*react[rtype].new2 */
    }
    else{/*types of atoms are not adequate to form a bond*/
      return -1;
    }
    /*check association distances for associations of a2*/
    a2_a11_ct = checkAssociatedBond(a2, new2, a11, a11_c);
    if(a2_a11_ct==-1)return -1;
    a2_a12_ct = checkAssociatedBond(a2, new2, a12, a12_c); 
    if(a2_a12_ct==-1)return -1;

    /*check association distances for associations of a1*/
    a1_a21_ct = checkAssociatedBond(a1, new1, a21, a21_c);
    if(a1_a21_ct==-1)return -1;
    a1_a22_ct = checkAssociatedBond(a1, new1, a22, a22_c); 
    if(a1_a22_ct==-1)return -1;
  }
  else{/*reverse reaction, bond may break*/
    if(old1==react[rtype].new1 && old2==react[rtype].new2){
      new1=react[rtype].old1;/*allways refer to bonded types*/
      new2=react[rtype].old2;
    }
    else if(old1==react[rtype].new2 && old2==react[rtype].new1){
      new1=react[rtype].old2;
      new2=react[rtype].old1;
    }
    else{/*types of atoms are not adequate to break a bond*/
      return -1;
    }
    a2_a11_ct = checkDessociatedBond(a2, new2, a11, a11_c);
    if(a2_a11_ct==-1)return -1;
    a2_a12_ct = checkDessociatedBond(a2, new2, a12, a12_c); 
    if(a2_a12_ct==-1)return -1;    
      
    a1_a21_ct = checkDessociatedBond(a1, new1, a21, a21_c);
    if(a1_a21_ct==-1)return -1;
    a1_a22_ct = checkDessociatedBond(a1, new1, a22, a22_c); 
    if(a1_a22_ct==-1)return -1;
  }

  /*collision type for (i1,i2) after reaction (forward or reverse) happens*/
  if(revers)/*reverse reaction, bond will break*/
    ct_new=react[rtype].out;/*coll type of unbound pair*/
  else/*forward reaction, bond will form*/
    ct_new=react[rtype].in;

  a1_a2_coll=find_coll(i1,i2,&a1_a2_lt); /*print_list_type(a1_a2_lt);*/
  old_pot=coll[ct1].etot;
  new_pot=coll[ct_new].etot; /*energy after reaction*/

  if(old1!=new1){
    a1->c=new1;/*change the atom type of i1*/
    if(a1_a21_coll=find_coll(i1,ia21,&a1_a21_lt)){
      old_pot+=coll[a1_a21_coll->c].etot;
    }
    new_pot+=coll[a1_a21_ct].etot;
    if(a1_a22_coll=find_coll(i1,ia22,&a1_a22_lt)){
      old_pot+=coll[a1_a22_coll->c].etot;
    }
    new_pot+=coll[a1_a22_ct].etot;
  }
  if(old2!=new2){
    a2->c=new2;/*change the atom type of i2*/
    if(a2_a11_coll=find_coll(i2,ia11,&a2_a11_lt)){
      old_pot+=coll[a2_a11_coll->c].etot;
    }
    new_pot+=coll[a2_a11_ct].etot;
    if(a2_a12_coll=find_coll(i2,ia12,&a2_a12_lt)){
      old_pot+=coll[a2_a12_coll->c].etot;
    }
    new_pot+=coll[a2_a12_ct].etot;
  }
  

  m1=a1->m;  m2=a2->m;
  du=new_pot-old_pot;
  if(!react[rtype].bond)du+=react[rtype].eo;
  duc=du*corr_2;
  ed=2*duc/(m1*m2*coll[ct].dm);
  di=1.0+ed/(sc*sc);
  if(di<=0){
    ab=-2.0*sc*coll[ct].dm;
    a1->c=old1;
    a2->c=old2;
    if(!revers)return -1;           
    ct_new=ct1;
  }
  else{
    ab=sc*coll[ct].dm*(sqrt(di)-1.0);
    vvm+=duc;
    potential+=du;
    if( (react[rtype].old1!=react[rtype].new1) ||
	(react[rtype].old2!=react[rtype].new2)   )
      setNewTypes(1);
    if(revers){
      /*printf("breaking HB %4d (%4d %4d) %4d (%4d %4d)\n",
	i1,ia11,ia12,i2,ia21,ia22);*/
      /*breakbonds. After breaking auxiliary bond (i,aux), it may be
	that i and aux are futher away than the maximum shell order
	for a non-bound pair of i and aux atoms. Thus, we have to
	remove the collision. If aux<i in the circular sense,
	collision (i,aux) could have been winning collision for plists
	of aux. In that case, we have to update the winning collision
	for the plists of aux. (i1,aux) cannot be the winning
	collision of plists of i1 because (i1,i2) is the winning
	collision. However, (i2,aux) can be the winning collision of
	the plist of i2 if i2<aux in the circular sense. If we do not
	have to remove collision (i,aux), then we don't look if
	(i,aux) was the winning coll. of the plist of aux, because
	we'll do that later, in update_tables, when we update all
	coll. times for collsions involving i1 or i2. However, we
	check whether we have to change the list type. For example,
	(i1,aux) may correspond to a HOT collision when bonded, and
	WARM or COLD when unbound*/
      breakBond(i1,i2); 
      mso=return_mso(i1,i2);/*mso of unbound i1 and i2, since we broke bond*/
      /*We have return_so(i1,i2)==return_mso(i1,i2)*/
      a1_a2_coll->c=ct_new;a1_a2_coll->mso=mso; 
      /*Putative move to a different list type is done in 
	update_collision_times*/
      
      breakBond(i1,ia21); mso=return_mso(i1,ia21);
      if(return_so(i1,ia21)>mso){/*too far away, remove collision*/
	if(less_than(i1,ia21))
	  /*No need to check whether a1_a21_coll was the winning coll
	    of plists[a1_a21_lt][i1], because we'll find the winning
	    collision later in udpate_collision_times*/
	  remove_neighbor_noupdlist(i1,a1_a21_lt,a1_a21_coll);
	else if(remove_neighbor(ia21,a1_a21_lt,a1_a21_coll))
	  /*the collision was the winning collision of ia12. We have
	    to mark here an update of winning coll time for ia21
	    because ia21 is not a anymore a neighbor of i1. Thus, it
	    will not be checked in update_collision_times for i1. Note
	    that if a1_a21_lt==COLD, then remove_neighbor will return 0*/
	  calendar.update_node_for[calendar.n_updates++]=ia21;
      }
      else{/*We do not remove coll., but maybe change its list type*/
	a1_a21_coll->c=a1_a21_ct;a1_a21_coll->mso=mso;
	new_lt=find_lt(i1,ia21,a1_a21_ct);
	/*If the list type has changed after change of coll type, then
	  we have to move the collision to the new list type*/
	if(new_lt!=a1_a21_lt){
	  mv_coll(a1_a21_coll,a1_a21_lt,new_lt);
	  if(less_than(ia21,i1) && a1_a21_lt!=COLD)
	    /*a1_a21_coll may have been first coll of list a1_a21_lt
	      of ia21. In that case we have to find the collision in
	      this list with the smallest coll time. This will not be
	      done later in update_collision_times, because we moved
	      a1_a21_coll to list type new_lt, and it is this list
	      that update_collision_times will automatically update.
	      In the rare case that the bond (i1,ia21) was adscribed
	      to the COLD list, none of this is neccessary because
	      COLD list does not have a winning time defined*/
	    if(tlists[a1_a21_lt][ia21]==a1_a21_coll){
	      update_tlist(a1_a21_lt,ia21);
	      if(winlist[ia21]==a1_a21_lt){
		/*a1_a21_coll might have been winning coll of plists of ia21.
		  In that case, we have to forward this fact to function
		  update_collision_times. However, this function is going to 
		  look for a1_a21_coll in list type new_lt, and not in old 
		  list type a1_a21_lt. Thus, in order to force an update of
		  the winnnig collision, we make a1_a21_coll to be the winner
		  collision of all plists of ia21. This is wrong for now, but
		  will be corrected later in update_collision_times*/
		if(new_lt!=COLD){
		  winlist[ia21]=new_lt;
		  tlists[new_lt][ia21]=a1_a21_coll;
		}
		/*if collision (i1,ia21) ends up in the COLD list,
		  after the bond is broken, then the coll time will
		  not be calculated in update_collision_times (because
		  this function does not look at the COLD list), and
		  the previous forwarding scheme will not work.  We
		  thus signal right here that ia21 is to be updated*/
		else{
		  update_winlist(ia21);
		  calendar.update_node_for[calendar.n_updates++]=ia21;
		}
	      }
	    }
	}
      }

      breakBond(i1,ia22); mso=return_mso(i1,ia22);
      if(return_so(i1,ia22)>mso){
	if(less_than(i1,ia22))
	  remove_neighbor_noupdlist(i1,a1_a22_lt,a1_a22_coll);
	else if(remove_neighbor(ia22,a1_a22_lt,a1_a22_coll))
	  calendar.update_node_for[calendar.n_updates++]=ia22;
      }
      else{
	a1_a22_coll->c=a1_a22_ct;a1_a22_coll->mso=mso;
	new_lt=find_lt(i1,ia22,a1_a22_ct);
	if(new_lt!=a1_a22_lt){
	  mv_coll(a1_a22_coll,a1_a22_lt,new_lt);
	  if(less_than(ia22,i1) && a1_a22_lt!=COLD)
	    if(tlists[a1_a22_lt][ia22]==a1_a22_coll){
	      update_tlist(a1_a22_lt,ia22);
	      if(winlist[ia22]==a1_a22_lt){
		if(new_lt!=COLD){
		  winlist[ia22]=new_lt;
		  tlists[new_lt][ia22]=a1_a22_coll;
		}
		else{
		  update_winlist(ia22);
		  calendar.update_node_for[calendar.n_updates++]=ia22;
		}
	      }
	    }
	}
       }

      breakBond(i2,ia11); mso=return_mso(i2,ia11);
      if(return_so(i2,ia11)>mso){
	if(less_than(i2,ia11))
	  remove_neighbor_noupdlist(i2,a2_a11_lt,a2_a11_coll);
	else if(remove_neighbor(ia11,a2_a11_lt,a2_a11_coll))
	  calendar.update_node_for[calendar.n_updates++]=ia11;
      }
      else{
	a2_a11_coll->c=a2_a11_ct;a2_a11_coll->mso=mso;
	new_lt=find_lt(i2,ia11,a2_a11_ct);
	if(new_lt!=a2_a11_lt){
	  mv_coll(a2_a11_coll,a2_a11_lt,new_lt);
	  if(less_than(ia11,i2) && a2_a11_lt!=COLD)
	    if(tlists[a2_a11_lt][ia11]==a2_a11_coll){
	      update_tlist(a2_a11_lt,ia11);
	      if(winlist[ia11]==a2_a11_lt){
		if(new_lt!=COLD){
		  winlist[ia11]=new_lt;
		  tlists[new_lt][ia11]=a2_a11_coll;
		}
		else{
		  update_winlist(ia11);
		  calendar.update_node_for[calendar.n_updates++]=ia11;
		}
	      }
	    }
	}
      }

      breakBond(i2,ia12); mso=return_mso(i2,ia12);
      if(return_so(i2,ia12)>mso){
	if(less_than(i2,ia12))
	  remove_neighbor_noupdlist(i2,a2_a12_lt,a2_a12_coll);
	else if(remove_neighbor(ia12,a2_a12_lt,a2_a12_coll))
	    calendar.update_node_for[calendar.n_updates++]=ia12;
      }
      else{
	a2_a12_coll->c=a2_a12_ct;a2_a12_coll->mso=mso;
	new_lt=find_lt(i2,ia12,a2_a12_ct);
	if(new_lt!=a2_a12_lt){	  
	  mv_coll(a2_a12_coll,a2_a12_lt,new_lt);
	  if(less_than(ia12,i2) && a2_a12_lt!=COLD)
	    if(tlists[a2_a12_lt][ia12]==a2_a12_coll){
	      update_tlist(a2_a12_lt,ia12);
	      if(winlist[ia12]==a2_a12_lt){
		if(new_lt!=COLD){
		  winlist[ia12]=new_lt;
		  tlists[new_lt][ia12]=a2_a12_coll;
		}
		else{
		  update_winlist(ia12);
		  calendar.update_node_for[calendar.n_updates++]=ia12;
		}
	      }
	    }
	}
      }

    }/*Matches if(revers)*/
    else if(react[rtype].bond){
      /*printf("forming HB %4d (%4d %4d) %4d (%4d %4d)\n",
	i1,ia11,ia12,i2,ia21,ia22);*/
      /*setBonds. When setting (i1,i2) bond, we don't have to check if
	we have to change the list type because we do that later in
	update_collision_times. When setting an auxiliary bond
	(i,aux), we don't need neither to compute the collision time,
	neither to update the winning time for the plists of aux,
	since we'll check these two items in update_collision_times,
	unless we move the collision to the COLD list. We may have to
	update the list type. For example, when setting bond
	(i1,ia21), we may change from WARM or COLD list type to HOT
	list type.*/
      mso=INF_SO; t=DBL2;

      setBond(i1,i2);
      /*Putative move to a different list type is done in
	update_collision_times*/
      a1_a2_coll->c=ct_new; a1_a2_coll->mso=mso;

      setBond(i1,ia21); new_lt=find_lt(i1,ia21,a1_a21_ct);
      if(a1_a21_coll){/*there was a collision entry for (i1,ia21)*/
	a1_a21_coll->c=a1_a21_ct;a1_a21_coll->mso=mso;
	if(new_lt!=a1_a21_lt){
	  mv_coll(a1_a21_coll,a1_a21_lt,new_lt);
	  /*a1_a21_coll may have been first coll of list a1_a21_lt of
	    ia21. In that case we have to find the collision in this
	    list with the smallest coll time. This will not be done
	    later in update_collision_times, because we moved
	    a1_a21_coll to list type new_lt, and it is this list that
	    update_collision_times will automatically update.  In the
	    rare case that the bond (i1,ia21) was adscribed to the
	    COLD list, none of this is neccessary because COLD list
	    does not have a winning time defined*/
	  if(a1_a21_lt!=COLD && less_than(ia21,i1))
	    if(tlists[a1_a21_lt][ia21]==a1_a21_coll){
	      update_tlist(a1_a21_lt,ia21);
	      if(winlist[ia21]==a1_a21_lt){
		/*a1_a21_coll might have been winning coll of plists of ia21.
		  In that case, we have to forward this fact to function
		  update_collision_times. However, this function is going to 
		  look for a1_a21_coll in list type new_lt, and not in old 
		  list type a1_a21_lt. Thus, in order to force an update of
		  the winnnig collision, we make a1_a21_coll to be the winner
		  collision of all plists of ia21. This is wrong for now, but
		  will be corrected later in update_collision_times*/
		if(new_lt!=COLD){
		  winlist[ia21]=new_lt;
		  tlists[new_lt][ia21]=a1_a21_coll;
		}
		else{
		/*if collision (i1,ia21) ends up in the COLD list,
		  after the bond is formed (very rare case), then the
		  coll time will not be calculated in
		  update_collision_times (because this function does
		  not look at the COLD list), and the previous
		  forwarding scheme will not work.  We thus signal
		  right here that ia21 is to be updated*/
		  update_winlist(ia21);
		  calendar.update_node_for[calendar.n_updates++]=ia21;
		}
	      }
	    }
	}
      }
      /*i1 and ia21 were further away than mso for unbound i1 and
	ia21.  insert_collision2(..) takes care if i1<ia21 in the
	circular sense. If ia21<i1 in the circular sense, we do not
	have to find the winning collision of plists[new_lt][ia21],
	since we do that later in update_collision_times(..). If we
	insert the coll onto the COLD list, then
	update_collision_times will not look at this list, since it's
	not neccessary. Neigther we have to calculate collision time
	t. Again, it'll be done in update_collision_times*/
      else insert_collision2(new_lt,t,i1,ia21,a1_a21_ct,mso);

      setBond(i1,ia22); new_lt=find_lt(i1,ia22,a1_a22_ct);
      if(a1_a22_coll){
	a1_a22_coll->c=a1_a22_ct;a1_a22_coll->mso=mso;
	if(new_lt!=a1_a22_lt){
	  mv_coll(a1_a22_coll,a1_a22_lt,new_lt);
	  if(a1_a22_lt!=COLD && less_than(ia22,i1))
	    if(tlists[a1_a22_lt][ia22]==a1_a22_coll){
	      update_tlist(a1_a22_lt,ia22);
	      if(winlist[ia22]==a1_a22_lt){
		if(new_lt!=COLD){
		  winlist[ia22]=new_lt;
		  tlists[new_lt][ia22]=a1_a22_coll;
		}
		else{
		  update_winlist(ia22);
		  calendar.update_node_for[calendar.n_updates++]=ia22;
		}
	      }
	    }
	}
      }
      else insert_collision2(new_lt,t,i1,ia22,a1_a22_ct,mso);

      setBond(i2,ia11); new_lt=find_lt(i2,ia11,a2_a11_ct);
      if(a2_a11_coll){
	a2_a11_coll->c=a2_a11_ct;a2_a11_coll->mso=mso;
	if(new_lt!=a2_a11_lt){
	  mv_coll(a2_a11_coll,a2_a11_lt,new_lt);
	  if(a2_a11_lt!=COLD && less_than(ia11,i2))
	    if(tlists[a2_a11_lt][ia11]==a2_a11_coll){
	      update_tlist(a2_a11_lt,ia11);
	      if(winlist[ia11]==a2_a11_lt){
		if(new_lt!=COLD){
		  winlist[ia11]=new_lt;
		  tlists[new_lt][ia11]=a2_a11_coll;
		}
		else{
		  update_winlist(ia11);
		  calendar.update_node_for[calendar.n_updates++]=ia11;
		}
	      }
	    }
	}
      }
      else insert_collision2(new_lt,t,i2,ia11,a2_a11_ct,mso);

      setBond(i2,ia12); new_lt=find_lt(i2,ia12,a2_a12_ct);
      if(a2_a12_coll){
	a2_a12_coll->c=a2_a12_ct;a2_a12_coll->mso=mso;
	if(new_lt!=a2_a12_lt){
	  mv_coll(a2_a12_coll,a2_a12_lt,new_lt);
	  if(a2_a12_lt!=COLD && less_than(ia12,i2))
	    if(tlists[a2_a12_lt][ia12]==a2_a12_coll){
	      update_tlist(a2_a12_lt,ia12);
	      if(winlist[ia12]==a2_a12_lt){
		if(new_lt!=COLD){
		  winlist[ia12]=new_lt;
		  tlists[new_lt][ia12]=a2_a12_coll;
		}
		else{
		  update_winlist(ia12);
		  calendar.update_node_for[calendar.n_updates++]=ia12;
		}
	      }
	    }
	}
      }
      else insert_collision2(new_lt,t,i2,ia12,a2_a12_ct,mso);

    }
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
  /*printf("End of process_hbII(..)\n");*/
  return ct_new; 
}/*Matches process_hbII(..)*/
/*=======================================================================*/
/*DUPLICATE FOR DEBUGGING PURPOSES*/
coll_type process_hbII2 (atom* a,atom_number i1,atom_number i2,coll_type ct1,
			double sc,double x,double y,double z,coll_type ct,
			coll_type rtype,int revers){/*270 lines*/
  double m1,m2,ab,ab1,ab2,di,ed,du,duc,old_pot,new_pot,t;
  double a1_a21_dd,a1_a22_dd,a2_a11_dd,a2_a12_dd;
  atom *a1,*a2;
  int old1,old2,new1,new2,bond;  
  coll_type ct_new;
  coll_type a1_a21_ct,a1_a22_ct,a2_a11_ct,a2_a12_ct;
  atom_number ia11, ia12, ia21, ia22;
  atom_type a11_c,a12_c,a21_c,a22_c;
  atom *a11,*a12,*a21,*a22;
  alist *a1_a2_coll,*a1_a21_coll,*a1_a22_coll,*a2_a11_coll,*a2_a12_coll;
  enum list_type lt,new_lt,a1_a2_lt,a1_a21_lt,a1_a22_lt,a2_a11_lt,a2_a12_lt;
  shell_order mso,so;
  printf("process_hbII2(%d,%d)\n",i1,i2);

  a1=a+i1;  a2=a+i2;
  /*assotiations of a1*/
  ia11=a1->hb_a1;  a11=a + a1->hb_a1; a11_c=a11->c;
  ia12=a1->hb_a2;  a12=a + a1->hb_a2; a12_c=a12->c;
  /*assotiations of a2*/
  ia21=a2->hb_a1;  a21=a + a2->hb_a1; a21_c=a21->c;
  ia22=a2->hb_a2;  a22=a + a2->hb_a2; a22_c=a22->c;

  old1=a1->c;  old2=a2->c;

  if(!revers){/*no reverse reaction, bond may form*/
    /*determine who is N, who is C*/
    if(old1==react[rtype].old1 && old2==react[rtype].old2){
      new1=react[rtype].new1; new2=react[rtype].new2;
    }
    else if(old1==react[rtype].old2 && old2==react[rtype].old1){
      new1=react[rtype].new2;/*type of p after bonding, but refers to */
      new2=react[rtype].new1;/*react[rtype].new2 */
    }
    else{
      printf("types of atoms are not adequate to form a bond\n");
      printf("Returning ct_new=-1\n");
      return -1;
    }
    /*check association distances for associations of a2*/
    a2_a11_ct = checkAssociatedBond(a2, new2, a11, a11_c);
    if(a2_a11_ct==-1)return -1;
    a2_a12_ct = checkAssociatedBond(a2, new2, a12, a12_c); 
    if(a2_a12_ct==-1)return -1;

    /*check association distances for associations of a1*/
    a1_a21_ct = checkAssociatedBond(a1, new1, a21, a21_c);
    if(a1_a21_ct==-1)return -1;
    a1_a22_ct = checkAssociatedBond(a1, new1, a22, a22_c); 
    if(a1_a22_ct==-1)return -1;
  }
  else{/*reverse reaction, bond may break*/
    if(old1==react[rtype].new1 && old2==react[rtype].new2){
      new1=react[rtype].old1;/*allways refer to bonded types*/
      new2=react[rtype].old2;
    }
    else if(old1==react[rtype].new2 && old2==react[rtype].new1){
      new1=react[rtype].old2;
      new2=react[rtype].old1;
    }
    else{
      printf("types of atoms are not adequate to break a bond\n");
      printf("returning ct_new=-1\n");
      return -1;
    }
    a2_a11_ct = checkDessociatedBond(a2, new2, a11, a11_c);
    if(a2_a11_ct==-1)return -1;
    a2_a12_ct = checkDessociatedBond(a2, new2, a12, a12_c); 
    if(a2_a12_ct==-1)return -1;    
      
    a1_a21_ct = checkDessociatedBond(a1, new1, a21, a21_c);
    if(a1_a21_ct==-1)return -1;
    a1_a22_ct = checkDessociatedBond(a1, new1, a22, a22_c); 
    if(a1_a22_ct==-1)return -1;
  }

  /*collision type for (i1,i2) after reaction (forward or reverse) happens*/
  if(revers){
    ct_new=react[rtype].out;/*coll type of unbound pair*/
    printf("reverse reaction, bond will break. ct_new=%d\n",ct_new);
  }
  else/*forward reaction, bond will form*/
    ct_new=react[rtype].in;

  a1_a2_coll=find_coll(i1,i2,&a1_a2_lt); /*print_list_type(a1_a2_lt);*/
  old_pot=coll[ct1].etot;
  new_pot=coll[ct_new].etot; /*energy after reaction*/

  if(old1!=new1){
    a1->c=new1;/*change the atom type of i1*/
    if(a1_a21_coll=find_coll(i1,ia21,&a1_a21_lt)){
      old_pot+=coll[a1_a21_coll->c].etot;
    }
    new_pot+=coll[a1_a21_ct].etot;
    if(a1_a22_coll=find_coll(i1,ia22,&a1_a22_lt)){
      old_pot+=coll[a1_a22_coll->c].etot;
    }
    new_pot+=coll[a1_a22_ct].etot;
  }
  if(old2!=new2){
    a2->c=new2;/*change the atom type of i2*/
    if(a2_a11_coll=find_coll(i2,ia11,&a2_a11_lt)){
      old_pot+=coll[a2_a11_coll->c].etot;
    }
    new_pot+=coll[a2_a11_ct].etot;
    if(a2_a12_coll=find_coll(i2,ia12,&a2_a12_lt)){
      old_pot+=coll[a2_a12_coll->c].etot;
    }
    new_pot+=coll[a2_a12_ct].etot;
  }
  

  m1=a1->m;  m2=a2->m;
  du=new_pot-old_pot;
  if(!react[rtype].bond)du+=react[rtype].eo;
  duc=du*corr_2;
  ed=2*duc/(m1*m2*coll[ct].dm);
  di=1.0+ed/(sc*sc);
  if(di<=0){
    ab=-2.0*sc*coll[ct].dm;
    a1->c=old1;
    a2->c=old2;
    if(!revers){
      return -1;
    }
    ct_new=ct1;
  }
  else{
    ab=sc*coll[ct].dm*(sqrt(di)-1.0);
    vvm+=duc;
    potential+=du;
    if( (react[rtype].old1!=react[rtype].new1) ||
	(react[rtype].old2!=react[rtype].new2)   )
      setNewTypes(1);
    if(revers){
	printf("break i1=%4d(ia11=%4d ia12=%4d) i2=%4d(ia21=%4d ia22=%4d)\n",
	       i1,ia11,ia12,i2,ia21,ia22);
      /*breakbonds. After breaking auxiliary bond (i,aux), it may be
	that i and aux are futher away than the maximum shell order
	for a non-bound pair of i and aux atoms. Thus, we have to
	remove the collision. If aux<i in the circular sense,
	collision (i,aux) could have been winning collision for plists
	of aux. In that case, we have to update the winning collision
	for the plists of aux. (i1,aux) cannot be the winning
	collision of plists of i1 because (i1,i2) is the winning
	collision. However, (i2,aux) can be the winning collision of
	the plist of i2 if i2<aux in the circular sense. If we do not
	have to remove collision (i,aux), then we don't look if
	(i,aux) was the winning coll. of the plist of aux, because
	we'll do that later, in update_tables, when we update all
	coll. times for collsions involving i1 or i2. However, we
	check whether we have to change the list type. For example,
	(i1,aux) may correspond to a HOT collision when bonded, and
	WARM or COLD when unbound*/
      breakBond(i1,i2); 
      mso=return_mso(i1,i2);/*mso of unbound i1 and i2, since we broke bond*/
      /*We have return_so(i1,i2)==return_mso(i1,i2)*/
      a1_a2_coll->c=ct_new;a1_a2_coll->mso=mso; 
      /*Putative move to a different list type is done in 
	update_collision_times*/
      
      breakBond(i1,ia21); 
      mso=return_mso(i1,ia21); /*maximum shell order of unbound atoms*/
      if(return_so(i1,ia21)>mso){/*too far away, remove collision*/
	if(less_than(i1,ia21))
	  /*No need to check whether a1_a21_coll was the winning coll
	    of plists[a1_a21_lt][i1], because we'll find the winning
	    collision later in udpate_collision_times*/
	  remove_neighbor_noupdlist(i1,a1_a21_lt,a1_a21_coll);
	else if(remove_neighbor(ia21,a1_a21_lt,a1_a21_coll))
	  /*the collision was the winning collision of ia12. We have
	    to mark here an update of winning coll time for ia21
	    because ia21 is not a anymore a neighbor of i1. Thus, it
	    will not be checked in update_collision_times for i1. Note
	    that if a1_a21_lt==COLD, then remove_neighbor will return 0*/
	  calendar.update_node_for[calendar.n_updates++]=ia21;
      }
      else{/*We do not remove coll., but maybe change its list type*/
	a1_a21_coll->c=a1_a21_ct;a1_a21_coll->mso=mso;
	new_lt=find_lt(i1,ia21,a1_a21_ct);
	/*If the list type has changed after change of coll type, then
	  we have to move the collision to the new list type*/
	if(new_lt!=a1_a21_lt){
	  mv_coll(a1_a21_coll,a1_a21_lt,new_lt);
	  if(less_than(ia21,i1) && a1_a21_lt!=COLD)
	    /*a1_a21_coll may have been first coll of list a1_a21_lt
	      of ia21. In that case we have to find the collision in
	      this list with the smallest coll time. This will not be
	      done later in update_collision_times, because we moved
	      a1_a21_coll to list type new_lt, and it is this list
	      that update_collision_times will automatically update.
	      In the rare case that the bond (i1,ia21) was previouly adscribed
	      to the COLD list, none of this is neccessary because
	      COLD list does not have a winning time defined*/
	    if(tlists[a1_a21_lt][ia21]==a1_a21_coll){
	      update_tlist(a1_a21_lt,ia21);
	      if(winlist[ia21]==a1_a21_lt){
		/*a1_a21_coll might have been winning coll of plists of ia21.
		  In that case, we have to forward this fact to function
		  update_collision_times. However, this function is going to 
		  look for a1_a21_coll in list type new_lt, and not in old 
		  list type a1_a21_lt. Thus, in order to force an update of
		  the winnnig collision, we make a1_a21_coll to be the winner
		  collision of all plists of ia21. This is wrong for now, but
		  will be corrected later in update_collision_times*/
		if(new_lt!=COLD){
		  winlist[ia21]=new_lt;
		  tlists[new_lt][ia21]=a1_a21_coll;
		}
		/*if collision (i1,ia21) ends up in the COLD list,
		  after the bond is broken, then the coll time will
		  not be calculated in update_collision_times (because
		  this function does not look at the COLD list), and
		  the previous forwarding scheme will not work.  We
		  thus signal right here that ia21 is to be updated*/
		else{
		  update_winlist(ia21);
		  calendar.update_node_for[calendar.n_updates++]=ia21;
		}
	      }
	    }
	}
      }

      dump_info_of(i1); dump_info_of(ia22);
      breakBond2(i1,ia22);
      mso=return_mso(i1,ia22); 
      printf("mso=%d so=%d\n",mso,return_so(i1,ia22));
      if(return_so(i1,ia22)>mso){
	if(less_than(i1,ia22))
	  remove_neighbor_noupdlist(i1,a1_a22_lt,a1_a22_coll);
	else if(remove_neighbor(ia22,a1_a22_lt,a1_a22_coll))
	  calendar.update_node_for[calendar.n_updates++]=ia22;
      }
      else{
	printf("Before: a1_a22_coll->c=%d a1_a22_coll->mso=%d\n",
	       a1_a22_coll->c,a1_a22_coll->mso);
	a1_a22_coll->c=a1_a22_ct;a1_a22_coll->mso=mso;
	printf("After: a1_a22_coll->c=%d a1_a22_coll->mso=%d\n",
	       a1_a22_coll->c,a1_a22_coll->mso);
	new_lt=find_lt2(i1,ia22,a1_a22_ct); print_list_type(new_lt);
	if(new_lt!=a1_a22_lt){
	  mv_coll2(a1_a22_coll,a1_a22_lt,new_lt);
	  if(less_than(ia22,i1) && a1_a22_lt!=COLD)
	    if(tlists[a1_a22_lt][ia22]==a1_a22_coll){
	      printf("a1_a22_coll winner for list type %d of atom %d\n",
		     a1_a22_lt,ia22);
	      update_tlist(a1_a22_lt,ia22);
	      if(winlist[ia22]==a1_a22_lt){
		if(new_lt!=COLD){
		  winlist[ia22]=new_lt;
		  tlists[new_lt][ia22]=a1_a22_coll;
		}
		else{
		  update_winlist(ia22);
		  calendar.update_node_for[calendar.n_updates++]=ia22;
		}
	      }
	    }
	}
       }

      breakBond(i2,ia11); mso=return_mso(i2,ia11);
      if(return_so(i2,ia11)>mso){
	if(less_than(i2,ia11))
	  remove_neighbor_noupdlist(i2,a2_a11_lt,a2_a11_coll);
	else if(remove_neighbor(ia11,a2_a11_lt,a2_a11_coll))
	  calendar.update_node_for[calendar.n_updates++]=ia11;
      }
      else{
	a2_a11_coll->c=a2_a11_ct;a2_a11_coll->mso=mso;
	new_lt=find_lt(i2,ia11,a2_a11_ct);
	if(new_lt!=a2_a11_lt){
	  mv_coll(a2_a11_coll,a2_a11_lt,new_lt);
	  if(less_than(ia11,i2) && a2_a11_lt!=COLD)
	    if(tlists[a2_a11_lt][ia11]==a2_a11_coll){
	      update_tlist(a2_a11_lt,ia11);
	      if(winlist[ia11]==a2_a11_lt){
		if(new_lt!=COLD){
		  winlist[ia11]=new_lt;
		  tlists[new_lt][ia11]=a2_a11_coll;
		}
		else{
		  update_winlist(ia11);
		  calendar.update_node_for[calendar.n_updates++]=ia11;
		}
	      }
	    }
	}
      }

      breakBond(i2,ia12); mso=return_mso(i2,ia12);
      if(return_so(i2,ia12)>mso){
	if(less_than(i2,ia12))
	  remove_neighbor_noupdlist(i2,a2_a12_lt,a2_a12_coll);
	else if(remove_neighbor(ia12,a2_a12_lt,a2_a12_coll))
	    calendar.update_node_for[calendar.n_updates++]=ia12;
      }
      else{
	a2_a12_coll->c=a2_a12_ct;a2_a12_coll->mso=mso;
	new_lt=find_lt(i2,ia12,a2_a12_ct);
	if(new_lt!=a2_a12_lt){	  
	  mv_coll(a2_a12_coll,a2_a12_lt,new_lt);
	  if(less_than(ia12,i2) && a2_a12_lt!=COLD)
	    if(tlists[a2_a12_lt][ia12]==a2_a12_coll){
	      update_tlist(a2_a12_lt,ia12);
	      if(winlist[ia12]==a2_a12_lt){
		if(new_lt!=COLD){
		  winlist[ia12]=new_lt;
		  tlists[new_lt][ia12]=a2_a12_coll;
		}
		else{
		  update_winlist(ia12);
		  calendar.update_node_for[calendar.n_updates++]=ia12;
		}
	      }
	    }
	}
      }

    }/*Matches if(revers)*/
    else if(react[rtype].bond){
      /*printf("forming HB %4d (%4d %4d) %4d (%4d %4d)\n",
	i1,ia11,ia12,i2,ia21,ia22);*/
      /*setBonds. When setting (i1,i2) bond, we don't have to check if
	we have to change the list type because we do that later in
	update_collision_times. When setting an auxiliary bond
	(i,aux), we don't need neither to compute the collision time,
	neither to update the winning time for the plists of aux,
	since we'll check these two items in update_collision_times,
	unless we move the collision to the COLD list. We may have to
	update the list type. For example, when setting bond
	(i1,ia21), we may change from WARM or COLD list type to HOT
	list type.*/
      mso=INF_SO; t=DBL2;

      setBond(i1,i2);
      /*Putative move to a different list type is done in
	update_collision_times*/
      a1_a2_coll->c=ct_new; a1_a2_coll->mso=mso;

      setBond(i1,ia21); new_lt=find_lt(i1,ia21,a1_a21_ct);
      if(a1_a21_coll){/*there was a collision entry for (i1,ia21)*/
	a1_a21_coll->c=a1_a21_ct;a1_a21_coll->mso=mso;
	if(new_lt!=a1_a21_lt){
	  mv_coll(a1_a21_coll,a1_a21_lt,new_lt);
	  /*a1_a21_coll may have been first coll of list a1_a21_lt of
	    ia21. In that case we have to find the collision in this
	    list with the smallest coll time. This will not be done
	    later in update_collision_times, because we moved
	    a1_a21_coll to list type new_lt, and it is this list that
	    update_collision_times will automatically update.  In the
	    rare case that the bond (i1,ia21) was adscribed to the
	    COLD list, none of this is neccessary because COLD list
	    does not have a winning time defined*/
	  if(a1_a21_lt!=COLD && less_than(ia21,i1))
	    if(tlists[a1_a21_lt][ia21]==a1_a21_coll){
	      update_tlist(a1_a21_lt,ia21);
	      if(winlist[ia21]==a1_a21_lt){
		/*a1_a21_coll might have been winning coll of plists of ia21.
		  In that case, we have to forward this fact to function
		  update_collision_times. However, this function is going to 
		  look for a1_a21_coll in list type new_lt, and not in old 
		  list type a1_a21_lt. Thus, in order to force an update of
		  the winnnig collision, we make a1_a21_coll to be the winner
		  collision of all plists of ia21. This is wrong for now, but
		  will be corrected later in update_collision_times*/
		if(new_lt!=COLD){
		  winlist[ia21]=new_lt;
		  tlists[new_lt][ia21]=a1_a21_coll;
		}
		else{
		/*if collision (i1,ia21) ends up in the COLD list,
		  after the bond is formed (very rare case), then the
		  coll time will not be calculated in
		  update_collision_times (because this function does
		  not look at the COLD list), and the previous
		  forwarding scheme will not work.  We thus signal
		  right here that ia21 is to be updated*/
		  update_winlist(ia21);
		  calendar.update_node_for[calendar.n_updates++]=ia21;
		}
	      }
	    }
	}
      }
      /*i1 and ia21 were further away than mso for unbound i1 and
	ia21.  insert_collision2(..) takes care if i1<ia21 in the
	circular sense. If ia21<i1 in the circular sense, we do not
	have to find the winning collision of plists[new_lt][ia21],
	since we do that later in update_collision_times(..). If we
	insert the coll onto the COLD list, then
	update_collision_times will not look at this list, since it's
	not neccessary. Neigther we have to calculate collision time
	t. Again, it'll be done in update_collision_times*/
      else insert_collision2(new_lt,t,i1,ia21,a1_a21_ct,mso);

      setBond(i1,ia22); new_lt=find_lt(i1,ia22,a1_a22_ct);
      if(a1_a22_coll){
	a1_a22_coll->c=a1_a22_ct;a1_a22_coll->mso=mso;
	if(new_lt!=a1_a22_lt){
	  mv_coll(a1_a22_coll,a1_a22_lt,new_lt);
	  if(a1_a22_lt!=COLD && less_than(ia22,i1))
	    if(tlists[a1_a22_lt][ia22]==a1_a22_coll){
	      update_tlist(a1_a22_lt,ia22);
	      if(winlist[ia22]==a1_a22_lt){
		if(new_lt!=COLD){
		  winlist[ia22]=new_lt;
		  tlists[new_lt][ia22]=a1_a22_coll;
		}
		else{
		  update_winlist(ia22);
		  calendar.update_node_for[calendar.n_updates++]=ia22;
		}
	      }
	    }
	}
      }
      else insert_collision2(new_lt,t,i1,ia22,a1_a22_ct,mso);

      setBond(i2,ia11); new_lt=find_lt(i2,ia11,a2_a11_ct);
      if(a2_a11_coll){
	a2_a11_coll->c=a2_a11_ct;a2_a11_coll->mso=mso;
	if(new_lt!=a2_a11_lt){
	  mv_coll(a2_a11_coll,a2_a11_lt,new_lt);
	  if(a2_a11_lt!=COLD && less_than(ia11,i2))
	    if(tlists[a2_a11_lt][ia11]==a2_a11_coll){
	      update_tlist(a2_a11_lt,ia11);
	      if(winlist[ia11]==a2_a11_lt){
		if(new_lt!=COLD){
		  winlist[ia11]=new_lt;
		  tlists[new_lt][ia11]=a2_a11_coll;
		}
		else{
		  update_winlist(ia11);
		  calendar.update_node_for[calendar.n_updates++]=ia11;
		}
	      }
	    }
	}
      }
      else insert_collision2(new_lt,t,i2,ia11,a2_a11_ct,mso);

      setBond(i2,ia12); new_lt=find_lt(i2,ia12,a2_a12_ct);
      if(a2_a12_coll){
	a2_a12_coll->c=a2_a12_ct;a2_a12_coll->mso=mso;
	if(new_lt!=a2_a12_lt){
	  mv_coll(a2_a12_coll,a2_a12_lt,new_lt);
	  if(a2_a12_lt!=COLD && less_than(ia12,i2))
	    if(tlists[a2_a12_lt][ia12]==a2_a12_coll){
	      update_tlist(a2_a12_lt,ia12);
	      if(winlist[ia12]==a2_a12_lt){
		if(new_lt!=COLD){
		  winlist[ia12]=new_lt;
		  tlists[new_lt][ia12]=a2_a12_coll;
		}
		else{
		  update_winlist(ia12);
		  calendar.update_node_for[calendar.n_updates++]=ia12;
		}
	      }
	    }
	}
      }
      else insert_collision2(new_lt,t,i2,ia12,a2_a12_ct,mso);

    }
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
  printf("End of process_hbII2(..)\nReturn ct_new=%d\n",ct_new);
  return ct_new; 
}/*Matches process_hbII2(..)*/
/*=======================================================================*/
/*Returns collision type between a and asso after hydrogen bond has formed*/
/*if a and asso are too far away, then association is not possible*/
int checkAssociatedBond(atom* ap, int newType, atom* asso, int assoType){
  int nct;
  double rx, ry, rz;
  int irx, iry, irz;
  double dd;
  moveatom(asso);
  rx=ap->r.x-asso->r.x;
  ry=ap->r.y-asso->r.y;
  rz=ap->r.z-asso->r.z;
  irx=ap->i.x.i-asso->i.x.i;
  iry=ap->i.y.i-asso->i.y.i;
  irz=ap->i.z.i-asso->i.z.i;
  if      (irx> amso) rx-=bound[0].length;
  else if (irx<-amso) rx+=bound[0].length;
  if      (iry> amso) ry-=bound[1].length;
  else if (iry<-amso) ry+=bound[1].length;
  if      (irz> amso) rz-=bound[2].length;
  else if (irz<-amso) rz+=bound[2].length;   
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
}/*Matches checkAssociatedBond(..)*/
/*=======================================================================*/
int checkAssociatedBond2(atom* ap, int newType, atom* asso, int assoType){
  int nct;
  double rx, ry, rz;
  int irx, iry, irz;
  double dd;
  moveatom(asso);
  rx=ap->r.x-asso->r.x;
  ry=ap->r.y-asso->r.y;
  rz=ap->r.z-asso->r.z;
  irx=ap->i.x.i-asso->i.x.i;
  iry=ap->i.y.i-asso->i.y.i;
  irz=ap->i.z.i-asso->i.z.i;
  if      (irx> amso) rx-=bound[0].length;
  else if (irx<-amso) rx+=bound[0].length;
  if      (iry> amso) ry-=bound[1].length;
  else if (iry<-amso) ry+=bound[1].length;
  if      (irz> amso) rz-=bound[2].length;
  else if (irz<-amso) rz+=bound[2].length;   
  dd=rx*rx+ry*ry+rz*rz;
  nct=icoll[newType][assoType];
  if(nct<0)return nct;
  if(dd>coll[nct].dd || coll[nct].prev!=-1)return -1;
  nct=coll[nct].next;
  while(nct>-1){
    if(dd>coll[nct].dd){ 
      printf("d(%d,%d)=%lf\n",(int)(ap-a),(int)(asso-a),sqrt(dd));
      return nct;
    }
    nct = coll[nct].next;
  }
  return -1;
}/*Matches checkAssociatedBond2(..)*/
/*=========================================================================*/  
/*Return type between atom a and bonded assoc.n asso after reaction breaks*/
int checkDessociatedBond(atom* a,int newType,atom* asso,int assoType){
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
  if      (irx> amso) rx-=bound[0].length;
  else if (irx<-amso) rx+=bound[0].length;
  if      (iry> amso) ry-=bound[1].length;
  else if (iry<-amso) ry+=bound[1].length;
  if      (irz> amso) rz-=bound[2].length;
  else if (irz<-amso) rz+=bound[2].length;   
  dd=rx*rx+ry*ry+rz*rz;
  nct = ecoll[newType][assoType];
  while(nct>-1){
    if(dd>=coll[nct].dd)return nct;
    nct=coll[nct].next;
  }
  return -1;
}
/*=========================================================================*/
