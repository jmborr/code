#include <math.h>
#include <stdio.h>
#include <float.h>
#include "bcp.h"
#include "search.h"
#include "controls.h"

extern nevents;
extern shell_order_scheme soc;
extern CollisionData *coll;
tsearch search;
atom *a;
static well_type ** ecoll;
static FILE * fp;
alist ***tlists,***plists,***qlists;
enum list_type *winlist;
extern double Lh[3];
int n1;
extern int n2,n3,n4;
tree_scheme calendar;
int ncoll=0;
extern flag;

/*=========================================================================*/
/*Initialize several structures. Remember that the last collision of
  all p-lists and q-list is search.storage[0]. This collision is only
  used to signal the end of any particular p-lists or q-list that we
  may be navigating

 Each collision stores the interacting pair (p,q), the type "c" which
 defines the inner discontinuity of the potential barrier, the
 collision time "t", the maximum shell order, the pointer to the next
 collision in the p-list and in the 1-list

 search.begin is an array whose elements are pointers to the elements
 of search.storage, to keep track of which elements of search.storage
 can be reused. An element of search.storage contains a collision. If
 a collision does not exists anymore, then we have to signal that we
 can use that bit of memory for a new collision. The scheme
 implemented by search.begin and search.free keeps track of this.
*/
void initsearch (void){ 
  alist * pt;
  int i,j;
  int n=search.n;
  /*printf("initsearch()\n");*/
  for (i=0;i<n;i++){
    search.collp[i]=-1;
    search.collq[i]=-1;
    search.atom_storage[i].n=i;
    search.atom_storage[i].next=NULL;
  }
  search.np=0;
  search.nq=0;
  pt=search.storage;
  for(i=0;i<search.maxfree;i++){/*initialize elements of search.begin*/
    search.begin[i]=pt;
    pt++;
  }
  /*pt points to the beginning of storage. Then, initialize only the
    first element of storage*/
  pt=search.storage;
  pt->t=DBL2;          
  pt->p=-1;      
  pt->q=-1;      
  pt->c=-1;      
  pt->mso=-1;    
  pt->pnext=NULL;/*This will end the end of any plist*/
  pt->qnext=NULL;/*This will end the end of any qlist*/
  
/*free points to the second element search.begin, which in turn points
  to the second element of search.storage. We reserve the first
  element of search.storage for special purposes*/
  search.free=search.begin+1;
  search.begin[search.maxfree]=NULL;/*This signal we cannot allocate
				      more than maxfree collisions*/
  j=n*nlt;
  for(i=0;i<j;i++){          /*Initially, all lists contain only the*/
    search.plists[0][i]=pt;  /*end-of-list signaling collision*/
    search.qlists[0][i]=pt;    
  }
  j=n*nft;
  for(i=0;i<j;i++) 
    search.tlists[0][i]=pt;
  
  for(i=0;i<=search.z;i++){
    search.cell_atoms[i]=NULL;
  }
  /*printf("End of initsearch()\n");*/
}/*Matches void initsearch(void)*/

/*========================================================================*/
/*Allocate memory space for the several variables we will use*/
int allocsearch (int n){
 int i,j,k;
 int maxadd;
 dimensions * bound=get_bounds();
 /*printf("allocsearch()\n");*/
 a=get_atom();
 ecoll=get_ecoll();
 search.x=bound[0].period;
 search.y=search.x*bound[1].period;
 search.z=search.y*bound[2].period;
 maxadd=search.z;
 search.n=n; 
 search.nhalf =n>>1;
 if(maxadd>MAXADD){
   printf("ERROR: system size too big. Increase value of MAXADD in bcp.h\n");
   return 0;
 }
 
 if(!(search.collp=(int *)malloc(n*sizeof(int))))return 0;
 if(!(search.collq=(int *)malloc(n*sizeof(int))))return 0;
 
 if(!(search.atomp=(int *)malloc(n*sizeof(int))))return 0;    
 if(!(search.atomq=(int *)malloc(n*sizeof(int))))return 0;

 if(!(search.cell_atoms=(clist **)malloc((maxadd+1)*sizeof(clist *))))return 0;
 if(!(search.atom_storage=(clist *)malloc(n*sizeof(clist))))return 0;

 search.maxfree=NFREE;
 if(!(search.storage =(alist *)malloc(search.maxfree*sizeof(alist))))return 0;
 if(!(search.begin=(alist **)malloc((search.maxfree+1)*sizeof(alist *))))
   return 0;

 if(!(search.neighbors=(alist **)malloc((n)*sizeof(alist *))))return 0;
 for(k=0;k<n;k++) search.neighbors[k]=NULL;
 if(!(search.neig_lt=(enum list_type *)malloc(n*sizeof(enum list_type))))
   return 0;

 if(!(search.plists=(alist ***)malloc(nlt*sizeof(alist **))))return 0;
 if(!(search.plists[0]=(alist **)malloc(nlt*n*sizeof(alist *))))return 0;
 for(k=1;k<nlt;k++)
   search.plists[k]=search.plists[k-1]+n;
 if(!(search.qlists=(alist ***)malloc(nlt*sizeof(alist **))))return 0;
 if(!(search.qlists[0]=(alist **)malloc(nlt*n*sizeof(alist *))))return 0;
 for(k=1;k<nlt;k++)
   search.qlists[k]=search.qlists[k-1]+n;
 if(!(search.tlists=(alist ***)malloc(nft*sizeof(alist **))))return 0;
 if(!(search.tlists[0]=(alist **)malloc(nft*n*sizeof(alist *))))return 0;
 for(k=1;k<nft;k++)
   search.tlists[k]=search.tlists[k-1]+n;
 k=0;
 if(!(search.winlist=(enum list_type *)malloc(n*sizeof(enum list_type))))
   return 0;
 for(k=0;k<n;k++) /*initialize the winlist entries*/
   search.winlist[k]=WALL ;

 plists=search.plists;/*initialize these pointers*/
 qlists=search.qlists;
 tlists=search.tlists;
 winlist=search.winlist;
 /*printf("End of allocsearch()\n");*/
 return 1;
}/*Matches int allocsearch(int n)*/
/*=======================================================================*/
/* list atoms in a certain cell, returns number of atoms in a cell */
int list_atoms(size_al address,int *atomx){ 
  int i=0;
  clist *inpt1;
  for(inpt1=search.cell_atoms[address];inpt1;inpt1=inpt1->next)
    atomx[i++]=inpt1->n;
 return i;	
}
/*=========================================================================*/
/*search.cell_atoms contains all the linked lists that are associated with
  a particular cell address.*/
/*Insert at front new atom in list*/
void add_atom2cell (int i)
{
  size_al address=a[i].add;
  clist * pt=search.cell_atoms[address];
  search.cell_atoms[address]=&search.atom_storage[i];
  search.cell_atoms[address]->next=pt;
}
/*=========================================================================*/
void dump_cell_occupants(cell_index address){
  clist *pt;
  pt=search.cell_atoms[address];
  if(pt){
    printf("Occupants of cell %d are: ",address);
    while(pt){printf("%d ",pt->n);pt=pt->next;} 
    printf("\n");
  }
}
/*===============================================================*/
/*At the beginning of the simulation, find all collisions*/
int init_tables (void){
  cell_index i,i0,i1,i2,j,j0,j1,j2,k,k0,k1,k2,I,J,K;
  cell_index di,dj,dk,address,addressj,addressk;
  shell_order so,soi,soj,sok,sojk,mso,amso=soc.amso;
  dimensions *bound=get_bounds();
  atom_number p,q;
  int n=search.n-1,ip;
  atom_type p_type,q_type;
  coll_type ct;
  double t;
  enum list_type lt;
  alist *new_coll;
  /*printf("init_tables()\n");*/
  n1=search.n; /*initialization of global n1 (scope is this file)*/
  initsearch();
  I=bound[0].period; J=bound[1].period; K=bound[2].period;
  for(p=0;p<n1;p++){ 
    add_atom2cell(p); /*find the cell address to which atom p belongs*/
    q=twall(p,&t);/*q is the wall index, either X-axis, Y-axis or Z-axis*/
    insert_collision(WALL,t,p,q,-1,0);
    q=tsf(p,&t);/*q is the index of a safe-limit coll (n4)*/
    insert_collision(SAFE,t,p,q,-1.0,0);
  }

  /*find collisions for each atom*/
  for(p=0;p<n1;p++){
    /*printf("p=%d\n",p); dump_info_of(p);*/
    p_type=a[p].origc;
    i0=a[p].i.x.i; j0=a[p].i.y.i; k0=a[p].i.z.i;/*cell indexes where p is*/
    /*we will look (i2-i1+1)*(j2-j1+1)*(k2-k1+1) cells surrounding p*/
    i1=i0-amso; i2=i0+amso;
    j1=j0-amso; j2=j0+amso;
    k1=k0-amso; k2=k0+amso; 

    dk=-amso;
    for(k=k1;k<=k2;k++){
      if(k<0) addressk=I*J*(k+K); /*apply periodic boundary conditions*/
      else if(k>=K) addressk=I*J*(k-K);
      else addressk=I*J*k; /*no need to apply periodic bound. cond.*/
      sok=abs(dk); /*initial estimate of shell order*/
      dj=-amso;
      for(j=j1;j<=j2;j++){
	if(j<0) addressj=I*(j+J)+addressk;
	else if(j>=J) addressj=I*(j-J)+addressk;   
	else addressj=I*j+addressk;
	soj=abs(dj);
	sojk = sok>soj ? sok : soj; /*pick biggest as shell order*/
	di=-amso;
	for(i=i1;i<=i2;i++){
	  if(i<0) address=i+I+addressj;
	  else if(i>=I) address=i-I+addressj;
	  else address=i+addressj;
	  soi=abs(di);
	  so = sojk>soi ? sojk : soi;
	  /*single integer "address" is uniquely associated with
	    (i,j,k), ie, address=i+I*(j+J*k). Now we look for atoms
	    that are in cell "address"*/
	  search.np=list_atoms(address,search.atomp);/*atoms in the cell*/
	  for(ip=0;ip<search.np;ip++){
	    q=search.atomp[ip]; /*one of the atoms in cell "address"*/
	    if(less_than(p,q)){ /*p<q in the circularly linked atom index*/
	      /*if(p==42&&q==127){printf("%d %d\n",p,q);exit(1);}*/
	      q_type=a[q].origc;
	      mso=return_mso(p,q);/*Returns INF_SO if p,q bonded*/
	      if(so<=mso){ /*p and q are within interaction range*/
		/*ct is collision type that uniquely determines the
		  inner potential barrier of the particular potential
		  well corresponding to the actual distance between p
		  and q*/
		ct=collision_type(p,q);
		if(ct<0) return (int) ct; /*program will exit*/
		add_potential(ct); /*add potential energy*/
		/*find collision time and return list type (HOT,COLD,or WARM)*/
		lt=tball_sf(p,q,ct,&t);/*t=dblarg1 if no collision*/
		/*print_list_type(lt);*/
		new_coll=insert_collision(lt,t,p,q,ct,mso);
	      }
	    }
	  }
	  di++;
	}
	dj++;
      }
      dk++;
    }
  }

  for(p=0;p<n1;p++){
    for(lt=0;lt<nlt-1;lt++)/*ntl-1 because we don't update COLD list*/
      update_tlist(lt,p);/*update the first collision time of each list*/
    update_winlist(p); /*update the winner list*/
  }
  /*printf("End of init_tables()\n");exit(1);*/
  return 1;
}/*Matches init_tables()*/

/*=======================================================================*/
int less_than(atom_number i,atom_number j){
  if( (i<j && j<=i+search.nhalf) || (i>j && i>j+search.nhalf) ) return 1;
  return 0;
}
/*=======================================================================*/
/*it is assumed that p<q in the circular atom number list*/
alist * insert_collision (enum list_type lt, double t, atom_number p,
			 atom_number q, coll_type c, shell_order mso){
  alist *pt_list,*new_coll;
  char l[10];
  strcpy(l,print_char_list_type(lt));
  
  ncoll++;

  /*printf("insert_collision(..)\n");*/
  new_coll=*(search.free);
  if(!new_coll){/*last element of search.begin points to NULL,indicating 
		  no more space in search.storage to insert new collisions*/
    printf("Message from \"insert_collision\"...\n");
    printf("ERROR:no more space in search.storage for new collisions!\n");
    writetext(get_text_name());
    exit(0);
  }
  search.free++;/*search.free points to next component of search.begin,which in
		  turns points to the next empty memory space in storage */
  new_coll->t=t;   
  new_coll->p=p;
  new_coll->q=q;
  new_coll->c=c;
  new_coll->mso=mso;
  new_coll->pprev=NULL;
  pt_list=plists[lt][p];
  new_coll->pnext=pt_list;/*insert at front*/
  pt_list->pprev=new_coll;
  plists[lt][p]=new_coll;
  switch(lt){
  case HOT:
  case WARM:
  case COLD:
    {
      new_coll->qprev=NULL;
      pt_list=qlists[lt][q];
      new_coll->qnext=pt_list;/*insert at front*/
      pt_list->qprev=new_coll;
      qlists[lt][q]=new_coll;
      break;
    }
  case WALL:
  case SAFE:
    {
      break;
    }
  default :
    {
      printf("Message from \"insert_collision\"...\n");
      printf("ERROR: no match for specified collision list\n");
      exit(0);
    }
  }
  /*printf("End of insert_collision(..)\n");*/
  return new_coll;
}
/*=======================================================================*/
/*We do NOT assume that i<j in the circular atom number list*/
alist * insert_collision2(enum list_type lt, double t, atom_number i,
			 atom_number j, coll_type c, shell_order mso){
  if(less_than(i,j))
    return insert_collision(lt,t,i,j,c,mso);
  else
    return insert_collision(lt,t,j,i,c,mso);
}
/*=======================================================================*/
void remove_neighbor_noupdlist(atom_number p, enum list_type lt, alist *pt){
  alist *prev,*next;
  int f=0;
  char l[10];
  strcpy(l,print_char_list_type(lt));

  ncoll--;

  prev=pt->pprev; /*take the neighbor out the plist of pt->p*/
  next=pt->pnext;
  if(prev){ prev->pnext=next; next->pprev=prev; }
  else{plists[lt][pt->p]=next; next->pprev=NULL; }

  prev=pt->qprev;  /*take the neighbor out the qlist of q*/
  next=pt->qnext;
  if(prev){ prev->qnext=next; next->qprev=prev; }
  else{ qlists[lt][pt->q]=next; next->qprev=NULL; }
  
  search.free--;
  *(search.free)=pt;/*liberate the space in search.storage for later use */
}
/*=========================================================================*/
int remove_neighbor (atom_number p, enum list_type lt, alist *pt){

  remove_neighbor_noupdlist(p,lt,pt);

  if(lt!=COLD)/*cold list does not have coll times*/
    if(tlists[lt][p]==pt){/*removing first HOT/WARM coll of plist of q !*/
      update_tlist(lt,p);
      if(lt==winlist[p]){/*if removed coll was winning coll*/
	update_winlist(p);
	return 1;
      }
    }

  return 0;
}/*Matches remove_neighbor(..)*/
/*=========================================================================*/ 
enum list_type find_lt (atom_number p, atom_number q, coll_type ct){
  atom *ap,*aq;
  double x,y,z,dd;
  int ix,iy,iz;
  dimensions *bound=get_bounds();

  ap=a+p; aq=a+q;
  x=aq->sfr.x-ap->sfr.x;
  y=aq->sfr.y-ap->sfr.y;
  z=aq->sfr.z-ap->sfr.z;
  if      (x>= Lh[0]) x-=bound[0].length;
  else if (x< -Lh[0]) x+=bound[0].length;
  if      (y>= Lh[1]) y-=bound[1].length;
  else if (y< -Lh[1]) y+=bound[1].length;
  if      (z>= Lh[2]) z-=bound[2].length;
  else if (z< -Lh[2]) z+=bound[2].length;
  dd=x*x+y*y+z*z;
  return return_lt(dd,ct);
}
/*=========================================================================*/ 
/*DUPLICATE FOR DEBUGGING PURPOSES*/
enum list_type find_lt2 (atom_number p, atom_number q, coll_type ct){
  atom *ap,*aq;
  double x,y,z,dd;
  int ix,iy,iz;
  dimensions *bound=get_bounds();
  printf("Beginning of find_lt2(p=%d q=%d ct=%d)\n",p,q,ct);
  ap=a+p; aq=a+q;
  x=aq->sfr.x-ap->sfr.x;
  y=aq->sfr.y-ap->sfr.y;
  z=aq->sfr.z-ap->sfr.z;
  if      (x>= Lh[0]) x-=bound[0].length;
  else if (x< -Lh[0]) x+=bound[0].length;
  if      (y>= Lh[1]) y-=bound[1].length;
  else if (y< -Lh[1]) y+=bound[1].length;
  if      (z>= Lh[2]) z-=bound[2].length;
  else if (z< -Lh[2]) z+=bound[2].length;
  dd=x*x+y*y+z*z;
  printf("d(%d,%d)=%lf\n",p,q,sqrt(dd));
  printf("End of find_lt2(..)\n");
  return return_lt2(dd,ct);
}
/*=========================================================================*/ 
enum list_type return_lt (double dd,coll_type ct){
  enum list_type lt=coll[ct].lt;/*either HOT or WARM*/

  if(lt==WARM){ /*a WARM collision should be COLD if above the safe limit*/
    if(coll[ct].dpsf2<dd){/*lower potential step is COLD*/	  
      if(coll[ct].prev>-1){/*let's check the upper potential step*/
	if(coll[coll[ct].prev].dmsf2>dd)/*upper step is also COLD*/
	  return COLD;
      }
      else/*there's not upper potential step*/
	return COLD;
    }
  }
  return lt;
}
/*=========================================================================*/ 
/*DUPLICATE FOR DEBUGGING PURPOSES*/
enum list_type return_lt2 (double dd,coll_type ct){
  enum list_type lt=coll[ct].lt;/*either HOT or WARM*/

  if(lt==WARM){ /*a WARM collision should be COLD if above the safe limit*/
    printf("coll[%d].dpsf2=%lf dd=%lf\n",ct,coll[ct].dpsf2,dd);
    if(coll[ct].dpsf2<dd){/*lower potential step is COLD*/	  
      if(coll[ct].prev>-1){/*let's check the upper potential step*/
	if(coll[coll[ct].prev].dmsf2>dd)/*upper step is also COLD*/
	  return COLD;
      }
      else/*there's not upper potential step*/
	return COLD;
    }
  }
  return lt;
}
/*=========================================================================*/ 
int get_free(void){
return (int)(search.free-search.begin);
}

int get_maxfree(void){
return search.maxfree;
}

void set_maxfree(int a){
search.maxfree=a;
}
/*=========================================================================*/
/*alist **neighbors stores pointers to all HOT,WARM, and COLD
  collisions in which p participates. ap stores the atom number of the
  collision partner*/
int init_neighbor_list (alist **neighbors,enum list_type *neig_lt,
			atom_number p, int *atomp){
  static int nlta=3;
  enum list_type lta[3]={HOT,WARM,COLD};
  enum list_type lt;
  atom_number p2,q;
  alist *pt;
  int i,np=0;
  /*create complete list of neighbors of p that may be updated*/

  for(i=0;i<nlta;i++){
    lt=lta[i];

    pt=plists[lt][p];
    while(pt->pnext){
      q=pt->q;
      neig_lt[q]=lt;  neighbors[q]=pt; atomp[np++]=q;
      pt=pt->pnext;
    }

    pt=qlists[lt][p];
    while(pt->qnext){
      p2=pt->p;
      neig_lt[p2]=lt;  neighbors[p2]=pt; atomp[np++]=p2;
      pt=pt->qnext;
    }
  }

  return np;
}
/*========================================================================*/
/*remove p from atom list pertaining to old cell, and add to list pertaining
  to new cell*/
void update_cell_atoms (atom_number p, cell_index old, cell_index new){
  clist *pt,*prev;

  prev=search.cell_atoms[old];
  if(prev->n==p){/*p turned out to be first atom listed*/
    pt=prev;
    search.cell_atoms[old]=pt->next;
  }
  else{
    pt=prev->next;
    while(pt->n != p){
      prev=pt;
      pt=pt->next;
    }
    prev->next=pt->next;
  }
  pt->next=search.cell_atoms[new];
  search.cell_atoms[new]=pt;
}
/*=========================================================================*/
/*search for new neighbors by new surfaces probed, and remove old neighbors
  by old surfaces not probed anymore. Beware of bonded atoms in the new and
  old surfaces*/
void update_neighbor_lists (atom_number p,atom_number w){
  enum list_type lt,*neig_lt;
  atom_type p_type,q_type;
  shell_order *ab;/*all boundaries*/
  shell_order mso;/*particular boundary probed*/
  shell_order *mso_tptq;
  int nb,bi,iq,nq,A,B,C,A0,B0,address;
  int f=1,a,b,c,x,y,a00,b00,c00,a0,b0,c0,a1,b1,c1,*i,*j,*k;
  int *atomp=search.atomp, *atomq=search.atomq;
  alist **neighbors=search.neighbors, *pt;
  atom_number q,nhalf;
  well_type *ecoll_pq;
  coll_type ct;
  dimensions *bound=get_bounds();
  double t;
  /*get_atom()  instead of just "a" (for global atom *a), because name "a" is
    taken as local variable int a */
  atom *at=get_atom();
  atom *ap=at+p;

  if(w==n1){/*p crosses a boundary perpendicular to the YZ plane*/
    /*i represents the x-axis, j the y-axis, and k the z-axis*/
    i=&a;j=&b;k=&c;
    A=bound[0].period; B=bound[1].period; C=bound[2].period;
    if(ap->v.x<0) f=-1;
    a00=ap->i.x.i;    b00=ap->i.y.i;    c00=ap->i.z.i;
  }
  else if(w==n2){ /*printf("w=%d\n",w);*/
    /*i represents the z-axis, j the x-axis, and k the y-axis*/
    j=&a;k=&b;i=&c;
    A=bound[1].period; B=bound[2].period; C=bound[0].period;
    if(ap->v.y<0) f=-1;
    a00=ap->i.y.i;    b00=ap->i.z.i;    c00=ap->i.x.i;
  }
  else{ /*printf("w=%d\n",w);*/
    /*i represents the y-axis, j the z-axis, and k the x-axis*/
    k=&a;i=&b;j=&c;
    A=bound[2].period; B=bound[0].period; C=bound[1].period;
    if(ap->v.z<0) f=-1;
    a00=ap->i.z.i;    b00=ap->i.x.i;    c00=ap->i.y.i;
  }

  A0=bound[0].period; B0=bound[1].period;
  p_type=ap->origc;
  mso_tptq=soc.mso_tptq[p_type]; 
  ecoll_pq=ecoll[p_type];
  nhalf=search.nhalf;
  neig_lt=search.neig_lt;
  ab=soc.boundaries[p_type];/*stores the shell order boundaries*/
  nb=ab[0];/*number of boundaries to probe*/
  search.np=init_neighbor_list(neighbors,neig_lt,p,atomp); 
  
  for(bi=1;bi<=nb;bi++){/*bi shell-boundary index*/

    mso=ab[bi];
    a0=a00-f*(mso+1); /*old cell*/ 
    a1=a00+f*mso; /*new cell*/
    if(a0>=A)a0-=A; else if(a0<0)a0+=A;
    if(a1>=A)a1-=A; else if(a1<0)a1+=A;
    b0=b00-mso; b1=b00+mso;
    c0=c00-mso; c1=c00+mso;

    for(x=b0;x<=b1;x++){/*Can't do (b=b0;b<=b1;b++) because of bound. cond. */
      if(x>=B)b=x-B; else if(x<0)b=x+B; else b=x;
      for(y=c0;y<=c1;y++){
	if(y>=C)c=y-C; else if(y<0)c=y+C; else c=y;

	a=a0; /*old cell, remove atoms*/
	address=(*i)+A0*( (*j) + B0*(*k) );
	if(nq=list_atoms(address,atomq)){/*if there are atoms in the cell*/
	  for(iq=0;iq<nq;iq++){/*go through every atom in the cell*/
	    q=atomq[iq]; /*collision partner*/
	    pt=neighbors[q];
	    if(pt){/*q listed as neighbor*/
	      if(pt->mso==mso){/*not bonded and out of range*/
		if(less_than(p,q)){
		  remove_neighbor_noupdlist(p,neig_lt[q],pt);
		}
		else if(remove_neighbor(q,neig_lt[q],pt)){
		  calendar.update_node_for[calendar.n_updates++]=q;
		}
	      }
	    }
	  }
	}

	a=a1; /*new cell, add atoms*/
	address=(*i)+A*( (*j) + B*(*k) );
	if(nq=list_atoms(address,atomq)){/*if there are atoms in the cell*/
	  for(iq=0;iq<nq;iq++){/*go through every atom in the cell*/
	    q=atomq[iq]; /*collision partner*/
	    q_type=at[q].origc;/*type of q*/
	    if(mso_tptq[q_type]==mso){/*p and q within interaction range*/
	      /*If not bonded, q is new neighbor. Thus, the collision type
	       corresponds to the most external unbond potential step*/
	      if(!isFriend(p,q)){
		ct=ecoll_pq[q_type];
		lt=tball_sf(p,q,ct,&t);/*type of list and coll t*/
		if(less_than(p,q)){
		  insert_collision(lt,t,p,q,ct,mso);/*into lt plist of p*/
		}
		else{/*into lt plist of q (and into lt qlist of p) */
		  pt=insert_collision(lt,t,q,p,ct,mso);
		  if(update_tlists(pt,lt,q)){/*winner coll for q !*/
		    calendar.update_node_for[calendar.n_updates++]=q;
		  }
		}
	      }
	    }
	  }
	}

      }
    }
  }/*Matches for(bi=1;bi<=nb;bi++)*/

  update_tlist(HOT,p);
  update_tlist(WARM,p);
  pt=plists[WALL][p];
  pt->q=twall(p,&(pt->t));
  update_winlist(p);
  calendar.update_node_for[calendar.n_updates++]=p;
  /*reinitialize all neighbors[q] to NULL*/
  for(f=0;f<search.np;f++)
    neighbors[atomp[f]]=NULL;

  /*printf("End of update_neighbor_lists(...)\n");*/
}/*Matches update_neighbor_lists(...)*/
/*========================================================================*/
void mv_coll (alist *pt,enum list_type lt1,enum list_type lt2){
  alist *ppt=pt->pprev, *npt=pt->pnext;
  atom_number p=pt->p,q=pt->q;
  /*remove pt from plists[lt1][p]*/
  if(ppt){
    ppt->pnext=npt;
    npt->pprev=ppt;
  }
  else{/*pt is first collision of plists[lt1][p]*/
    plists[lt1][p]=npt;
    npt->pprev=0;
  }
  /*remove pt from qlists[lt1][q]*/
  ppt=pt->qprev; npt=pt->qnext;
  if(ppt){
    ppt->qnext=npt;
    npt->qprev=ppt;
  }
  else{
    qlists[lt1][q]=npt;
    npt->qprev=0;
  }
  /*insert-at-front pt in plists[lt2][p] and in qlists[lt2][q]*/
  npt=plists[lt2][p];
  pt->pprev=NULL;
  pt->pnext=npt;
  npt->pprev=pt;
  plists[lt2][p]=pt;
  
  npt=qlists[lt2][q];
  pt->qprev=NULL;
  pt->qnext=npt;
  npt->qprev=pt;
  qlists[lt2][q]=pt;
}/*Matches mv_coll(..)*/
/*========================================================================*/
/*DUPLICATE FOR DEBUGGING PURPOSES*/
void mv_coll2 (alist *pt,enum list_type lt1,enum list_type lt2){
  alist *ppt=pt->pprev, *npt=pt->pnext;
  atom_number p=pt->p,q=pt->q;
  printf("Beginning of mv_coll2(lt1,lt2)\n");
  printf("lt1 ");print_list_type(lt1); printf("lt2 "); print_list_type(lt2);
  printf("remove pt from plists[lt1=%d][p=%d]\n",lt1,p);
  if(ppt){
    ppt->pnext=npt;
    npt->pprev=ppt;
  }
  else{
    printf("pt is first collision of plists[lt1=%d][p=%d]\n",lt1,p);
    plists[lt1][p]=npt;
    npt->pprev=0;
  }
  printf("remove pt from qlists[lt1=%d][q=%d]\n",lt1,q);
  ppt=pt->qprev; npt=pt->qnext;
  if(ppt){
    ppt->qnext=npt;
    npt->qprev=ppt;
  }
  else{
    qlists[lt1][q]=npt;
    npt->qprev=0;
  }
  printf("insert-at-front in plists[lt2=%d][p=%d] and qlists[lt2=%d][q=%d]\n",
	 lt2,p,lt2,q);
  npt=plists[lt2][p];
  pt->pprev=NULL;
  pt->pnext=npt;
  npt->pprev=pt;
  plists[lt2][p]=pt;
  
  npt=qlists[lt2][q];
  pt->qprev=NULL;
  pt->qnext=npt;
  npt->qprev=pt;
  qlists[lt2][q]=pt;
}/*Matches mv_coll2(..)*/
/*========================================================================*/
/*requires previous update of time and position of atom p and
  safe-limit origin. This is actually done in collision_sf(..)*/
void update_sf (atom_number p){
  alist *pt,*opt;
  int nupd=0,i;
  atom_number r;

  /*printf("update_sf(..)\n",);*/
  /*update safe-limit collision*/
  pt=plists[SAFE][p];  pt->q=tsf(p,&(pt->t));

  pt=plists[WARM][p];
  while(pt->pnext){
    if(find_lt(p,pt->q,pt->c)==COLD){
      /*printf("moving WARM coll %d %d to COLD\n",p,pt->q);*/
      opt=pt;
      pt=pt->pnext;
      nupd++;
      mv_coll(opt,WARM,COLD); /*insert-at-front in COLD list*/
      /*coll enters COLD list of p with pt->t < INF, but we don't care*/
    }
    else pt=pt->pnext;
  }

  pt=plists[COLD][p];
  for(i=0;i<nupd;i++) pt=pt->pnext; /*don't check recently inserted*/
  while(pt->pnext){
    if(find_lt(p,pt->q,pt->c)==WARM){
      /*printf("moving COLD coll %d %d to WARM\n",p,pt->q);*/
      opt=pt;
      pt=pt->pnext;
      mv_coll(opt,COLD,WARM);
      tball(p,opt->q,opt->c,&(opt->t));/*find collision time*/
    }
    else pt=pt->pnext;
  }

  nupd=0;

  pt=qlists[WARM][p];
  while(pt->qnext){
    if(find_lt(pt->p,p,pt->c)==COLD){
      /*printf("moving COLD coll %d %d to WARM\n",pt->p,p);*/
      nupd++;
      opt=pt;
      pt=pt->qnext;
      mv_coll(opt,WARM,COLD); /*insert-at-front in COLD list*/
      /*coll enters COLD list of pt->p with pt->t < INF, but we don't care*/
    }
    else pt=pt->qnext;
  }

  pt=qlists[COLD][p];
  for(i=0;i<nupd;i++) pt=pt->qnext; /*don't check recently inserted*/
  while(pt->qnext){
    r=pt->p;
    if(find_lt(r,p,pt->c)==WARM){
      /*printf("moving COLD coll %d %d to WARM\n",pt->p,p);*/
      opt=pt;
      pt=pt->qnext;
      mv_coll(opt,COLD,WARM); 
      tball(r,p,opt->c,&(opt->t));
      if(update_tlists(opt,WARM,r))
	calendar.update_node_for[calendar.n_updates++]=r;
    }
    else pt=pt->qnext;
  }

  update_tlist(WARM,p);
  /*even though nupd=0, we have updated plists[SAFE][p], thus update_winlist*/
  update_winlist(p); /*find the plists of p with the smallest time*/
  calendar.update_node_for[calendar.n_updates++]=p;
  /*printf("End of update_sf(..)\n");*/
}
/*========================================================================*/
void update_collision_times (atom_number p, atom_number q,coll_type ct){
  enum list_type wlt,lt,olt;
  alist *pt,*opt;
  atom_number s,sa[2],r;
  int i,j;
  /*printf("update_collision_times(%d,%d)\n",p,q);*/
  wlt=winlist[p]; 
  pt=tlists[wlt][p];
  pt->c=ct;/*update collision type*/
  if(wlt!=coll[ct].lt)/*HOT <==> WARM transition after change in coll type*/
    mv_coll(pt,wlt,coll[ct].lt);

  /*update collision times of p*/
  pt=plists[WALL][p];     
  pt->q=twall(p,&(pt->t));
  pt=plists[SAFE][p];     
  pt->q=tsf(p,&(pt->t));  

  pt=plists[HOT][p]; /*update plists[HOT][p]*/ 
  while(pt->pnext){
    tball(p,pt->q,pt->c,&(pt->t));
    pt=pt->pnext;
  }
  update_tlist(HOT,p);/*search for first coll in HOT list*/ 

  pt=plists[WARM][p];/*update plists[WARM][p]*/
  while(pt->pnext){
    tball(p,pt->q,pt->c,&(pt->t));
    pt=pt->pnext;
  }
  update_tlist(WARM,p);

  pt=qlists[HOT][p];/*update qlists[HOT][p]*/
  while(pt->qnext){
    r=pt->p;
    tball(r,p,pt->c,&(pt->t));
    /*We have to check if r smaller than p and in that case, if
      collision (r,p) is a winner collision for plists[HOT] of r. If
      in addition, collision (r,p) is the absolute winner for all
      plists of r, then we have to update the node in the binary-tree
      for r*/
    if(update_tlists(pt,HOT,r))
      calendar.update_node_for[calendar.n_updates++]=r;
    pt=pt->qnext;
  }

  pt=qlists[WARM][p];/*update qlists[WARM][p]*/
  while(pt->qnext){
    r=pt->p;
    tball(r,p,pt->c,&(pt->t));
    /*We have to check if r smaller than p and in that case, if
      collision (r,p) is a winner collision for plists[WARM] of r. If
      in addition, collision (r,p) is the absolute winner for all
      plists of r, then we have to update the node in the binary-tree
      for r*/
    if(update_tlists(pt,WARM,r))
      calendar.update_node_for[calendar.n_updates++]=r;
    pt=pt->qnext;
  }

  update_winlist(p); /*find the plists of p with the smallest time*/
  calendar.update_node_for[calendar.n_updates++]=p;

  /*We don't update COLD list, since particles are too far away to collide*/
  /*Update collision times of q*/
  pt=plists[WALL][q];
  pt->q=twall(q,&(pt->t));
  pt=plists[SAFE][q];
  pt->q=tsf(q,&(pt->t));

  pt=plists[HOT][q];/*update plists[HOT][q]*/
  while(pt->pnext){
    tball(q,pt->q,pt->c,&(pt->t));
    pt=pt->pnext;
  }
  update_tlist(HOT,q);

  pt=plists[WARM][q];/*update plists[WARM][q]*/
  while(pt->pnext){
    tball(q,pt->q,pt->c,&(pt->t));
    pt=pt->pnext;
  }
  update_tlist(WARM,q);

  pt=qlists[HOT][q];/*update qlists[HOT][q]*/
  while(pt->qnext){
    r=pt->p;
    if(r!=p){/*r might be p, which we already updated*/
      tball(r,q,pt->c,&(pt->t));
      if(update_tlists(pt,HOT,r))
	calendar.update_node_for[calendar.n_updates++]=r;
    }
    pt=pt->qnext;
  }

  pt=qlists[WARM][q];/*update qlists[lt][s]*/
  while(pt->qnext){
    r=pt->p;
    if(r!=p){
      tball(r,q,pt->c,&(pt->t));
      if(update_tlists(pt,WARM,r))
	calendar.update_node_for[calendar.n_updates++]=r;
    }
    pt=pt->qnext;
  }
  
  update_winlist(q); /*find the plists[lt][p] with the smallest time*/
  calendar.update_node_for[calendar.n_updates++]=q;

  /*printf("End of update_collision_times(...)\n");*/
}/*Matches void update_collision_times(...)*/
/*========================================================================
DUPLICATE FOR DEBUGGING PURPOSES
*/
void update_collision_times2 (atom_number p, atom_number q,coll_type ct){
  enum list_type wlt,lt,olt;
  alist *pt,*opt;
  atom_number s,sa[2],r;
  int i,j;
  printf("update_collision_times2(p=%d,q=%d,ct=%d)\n",p,q,ct);
  wlt=winlist[p];
  pt=tlists[wlt][p];
  pt->c=ct;/*update collision type*/
  if(wlt!=coll[ct].lt){/*HOT <==> WARM transition after change in coll type*/
    /*print_list_type(wlt); print_list_type(coll[ct].lt);*/
    mv_coll(pt,wlt,coll[ct].lt);
  }
  printf("update collision times of p=%d\n",p);
  pt=plists[WALL][p];     
  pt->q=twall(p,&(pt->t));
  pt=plists[SAFE][p];     
  pt->q=tsf(p,&(pt->t));  

  pt=plists[HOT][p]; /*update plists[HOT][p]*/ 
  while(pt->pnext){
    tball(p,pt->q,pt->c,&(pt->t));
    pt=pt->pnext;
  }
  update_tlist(HOT,p);/*search for first coll in HOT list*/ 

  pt=plists[WARM][p];/*update plists[WARM][p]*/
  while(pt->pnext){
    tball(p,pt->q,pt->c,&(pt->t));
    pt=pt->pnext;
  }
  update_tlist(WARM,p);

  pt=qlists[HOT][p];/*update qlists[HOT][p]*/
  while(pt->qnext){
    r=pt->p;
    tball(r,p,pt->c,&(pt->t));
    /*We have to check if r smaller than p and in that case, if
      collision (r,p) is a winner collision for plists[HOT] of r. If
      in addition, collision (r,p) is the absolute winner for all
      plists of r, then we have to update the node in the binary-tree
      for r*/
    if(update_tlists(pt,HOT,r))
      calendar.update_node_for[calendar.n_updates++]=r;
    pt=pt->qnext;
  }

  pt=qlists[WARM][p];/*update qlists[WARM][p]*/
  while(pt->qnext){
    r=pt->p;
    tball(r,p,pt->c,&(pt->t));
    /*We have to check if r smaller than p and in that case, if
      collision (r,p) is a winner collision for plists[WARM] of r. If
      in addition, collision (r,p) is the absolute winner for all
      plists of r, then we have to update the node in the binary-tree
      for r*/
    if(update_tlists(pt,WARM,r))
      calendar.update_node_for[calendar.n_updates++]=r;
    pt=pt->qnext;
  }

  update_winlist(p); /*find the plists of p with the smallest time*/
  calendar.update_node_for[calendar.n_updates++]=p;

  /*We don't update COLD list, since particles are too far away to collide*/
  printf("Update collision times of q=%d\n",q);
  pt=plists[WALL][q];
  pt->q=twall(q,&(pt->t));
  pt=plists[SAFE][q];
  pt->q=tsf(q,&(pt->t));

  pt=plists[HOT][q];/*update plists[HOT][q]*/
  while(pt->pnext){
    tball(q,pt->q,pt->c,&(pt->t));
    pt=pt->pnext;
  }
  update_tlist(HOT,q);

  pt=plists[WARM][q];/*update plists[WARM][q]*/
  while(pt->pnext){
    tball(q,pt->q,pt->c,&(pt->t));
    pt=pt->pnext;
  }
  update_tlist(WARM,q);

  pt=qlists[HOT][q];/*update qlists[HOT][q]*/
  while(pt->qnext){
    r=pt->p;
    if(r!=p){/*r might be p, which we already updated*/
      tball(r,q,pt->c,&(pt->t));
      if(update_tlists(pt,HOT,r))
	calendar.update_node_for[calendar.n_updates++]=r;
    }
    pt=pt->qnext;
  }

  pt=qlists[WARM][q];/*update qlists[lt][s]*/
  while(pt->qnext){
    r=pt->p;
    if(r!=p){
      tball(r,q,pt->c,&(pt->t));
      if(update_tlists(pt,WARM,r))
	calendar.update_node_for[calendar.n_updates++]=r;
    }
    pt=pt->qnext;
  }
  
  update_winlist(q); /*find the plists[lt][p] with the smallest time*/
  calendar.update_node_for[calendar.n_updates++]=q;

  printf("End of update_collision_times2(...)\n");
}/*Matches void update_collision_times2(...)*/
/*========================================================================*/
/*update the first collision for list lt. If list empty, then "first" and "pt"
  both point to search.storage, whose first element has t=DBL2*/
void update_tlist (enum list_type lt, atom_number p){
  alist *pt,*first;
  double t;
  pt=first=plists[lt][p];
  t=first->t;
  while(pt->pnext){
    if(pt->t < t){ first=pt; t=first->t;}
    pt=pt->pnext;
  }
  tlists[lt][p]=first;
}
/*=================================================================*/
/*All possible scenarios after update of a collision. Returns 1 if collision
  becomes the first collision of all p-collisions of p*/
int update_tlists (alist *pt,enum list_type lt,atom_number p){
  enum list_type wlt;
  alist *ptlt,*ptwlt;

  if(lt!=COLD){
    wlt=winlist[p];
    ptwlt=tlists[wlt][p]; /*absolute winner of all plists of p*/
    ptlt=tlists[lt][p];

    if(pt==ptwlt){/*pt was the winning collision of all plists of p*/
      update_tlist(lt,p);/*update the particular list*/
      update_winlist(p); /*find the winning list*/
      return 1;
    }
    
    else if(pt==ptlt){/*pt was the winning collision of lt plists of p*/
      update_tlist(lt,p); /*update the particular list*/
      if(tlists[lt][p]->t < ptwlt->t){/*compare to winning list*/
	winlist[p]=lt;
	return 1;
      }
    }
    /*pt is not ptlt*/
    else if(pt->t < ptlt->t){/*pt becomes winning collision of lt list*/
      tlists[lt][p]=pt;
      if(pt->t < ptwlt->t){/*pt becomes winning coll of all plists of p*/
	winlist[p]=lt;
	return 1;
      }
    }
  }
  return 0;
}
/*=======================================================================*/
/*Find the list with the smallest collision time*/
void update_winlist (atom_number p){
  enum list_type lt,wlt;
  double t0,t;
  int i;
  /*printf("update_winlist(..)\n");*/
  wlt=winlist[p]; /*winlist has been initialized in allocsearch(..)*/
  t0=tlists[wlt][p]->t;
  for(lt=0;lt<nlt-1;lt++){ /*nlt-1 because last list (COLD) not collisions*/
    t=tlists[lt][p]->t;
    if(t<t0){
      wlt=lt;
      t0=tlists[wlt][p]->t;
    }
  }
  winlist[p]=wlt;
  /*printf("End of update_winlist(..)\n");*/
}
/*=======================================================================*/
void update_table (atom_number p, atom_number q, coll_type ct){
  if(q<n1){ /*printf("  interparticle collision...\n");*/
    update_collision_times(p,q,ct);
  }
  else if(q<n4){ /*printf("  wall collision...\n");*/
    newloc(a,p,q);/*update cell address of atom p*/
    update_neighbor_lists(p,q);
  }
  else{ /*printf(" safe-limit collision...\n");*/
    update_sf(p);/*Update COLD, WARM lists. Update safe-limit coll time*/
   }
  renew_calendar();
}/*Matches update_table(..)*/
/*=========================================================================*/
alist *find_coll (atom_number i, atom_number j, enum list_type *lt0){
  enum list_type lt;
  int ilt,nlta=3;
  enum list_type lta[3]={HOT,WARM,COLD};
  alist *pt;
  /*go through plists and qlists of i*/
  for(ilt=0;ilt<nlta;ilt++){
    lt=lta[ilt];/*specify the particular list type*/
    /*go through the plist*/
    pt=plists[lt][i];
    while(pt->pnext){
      if(pt->q==j){ *lt0=lt; return pt; }
      pt=pt->pnext;
    }
    /*go through the qlist*/
    pt=qlists[lt][i];
    while(pt->qnext){
      if(pt->p==j){ *lt0=lt; return pt; }
      pt=pt->qnext;
    }
  }
  return NULL;
}
/*=========================================================================*/
int allocalendar( int n1 ){
  calendar.n1=n1;
  if(!(calendar.storage=(node*)malloc(n1*sizeof(node)))               ||
     !(calendar.update_node_for=(atom_number*)malloc(n1*sizeof(atom_number))))
    return 0;
  return 1;
}
/*=========================================================================*/
void initcalendar (void){
  atom_number p;
  node *pt,*new_node,*current,*prev;
  alist ***tlists=search.tlists;
  enum list_type wlt,*winlist=search.winlist;
  
  double t;
  /*printf("initcalendar(...)\n");*/
  pt=calendar.storage; /*insert first node*/
  for(p=0;p<calendar.n1;p++){/*initialize elements of calendar.storage*/
    wlt=winlist[p];
    pt->parent=NULL;
    pt->t=tlists[wlt][p]->t;
    pt->fast_child=NULL;
    pt->slow_child=NULL;
    pt++;/*move to next element of storage*/
  }

  calendar.root=calendar.storage; /*initialize root to storage[0]*/

  /*we begin with p==1 because we inserted p==0 as root of the tree*/
  for(p=1;p<calendar.n1;p++)/*order events in binary tree scheme*/
    insert_node(p);

  find_new_champion();/*find absolute winner collision*/
}
/*=========================================================================*/
void find_new_champion (void){
  node *current,*prev;

  current=calendar.root;
  while(current){
    prev=current;
    current=current->fast_child;
  }
   /*Pointer arithmetic,returns index of storage such that storage[p]
     is the champion*/
  calendar.p=prev-calendar.storage;
}                            
/*=========================================================================*/
/*obtain winner pair p1,q1 and winning time timea*/
coll_type squeeze_table (atom_number *p1,atom_number *q1,double *timea){ 
  alist *inpt;
  enum list_type wlt;
  atom_number p=calendar.p;
  /*printf("squeeze_table(...)\n");*/
  wlt=winlist[p];
  inpt=tlists[wlt][p];
  *timea=inpt->t;
  *p1=inpt->p;
  *q1=inpt->q;
  return (int)(inpt->c);
  /*printf("End of squeeze_table(...)\n");*/
}/*Matches squeeze_table(..)*/
/*=========================================================================*/  
int * get_collp(void){return search.collp;}
int * get_collq(void){return search.collq;}
int * get_atomp(void){return search.atomp;}
int * get_atomq(void){return search.atomq;}
int get_np(void){return search.np;}
int get_nq(void){return search.nq;}

int pairs(int (*do_something)(int,int,int))
{ 
  int n_bonds=0;
  int i;
  alist *inpt1;
  for (inpt1=search.begin[1];inpt1<*search.free;inpt1++)
    if(inpt1->q<search.n)
      n_bonds+=do_something(inpt1->p,inpt1->q,inpt1->c);
  return n_bonds;
}  
/*========================================================================*/
char *atom_name (atom_type atc){
  static char name[8];
  if(atc==-1) strcpy(name,"SF");
  else if(atc==0) strcpy(name,"WALL");
  else if(atc==1||atc==3) strcpy(name,"N");
  else if(atc==2||atc==4) strcpy(name,"C");
  else if(atc==5) strcpy(name,"CA");
  else strcpy(name,"CB");
  return name;
}
/*=========================================================================*/
void print_CollisionData(coll_type ct){
  printf("prev=%3d next=%3d\n",coll[ct].prev,coll[ct].next);
}
/*=========================================================================*/
void print_coll (alist *pt){
  atom_type pc,qc;
  char namep[8],nameq[8];
  pc=a[pt->p].origc;
  if(pt->q >=n1 && pt->q<=n3) qc=0;
  else if(pt->q==n4) qc=-1;
  else qc=a[pt->q].origc;
  strcpy(namep,atom_name(pc));
  strcpy(nameq,atom_name(qc));

  if(pt==NULL){ printf("NULL collision!\n"); return ;}
  printf("pt->p=%d(%s) pt->q=%d(%s) pt->c=%d pt->mso=",
	 pt->p,namep,pt->q,nameq,pt->c);
  if(pt->mso==INF_SO) printf("INF_SO ");  else printf("%d ",pt->mso);
  printf("pt->t=");
  if(pt->t<DBL1) printf("%lf\n",pt->t);  else printf("INF\n");
  printf("  pt=%d pt->pprev=%d pt->pnext=",pt,pt->pprev);
  if(pt->pnext==search.storage) printf("search.storage ");
  else printf("%d ",pt->pnext);
  printf("pt->qprev=%d pt->qnext=",pt->qprev);
  if(pt->qnext==search.storage) printf("search.storage ");
  else printf("%d ",pt->qnext);
  printf("\n");
}
/*========================================================================*/
void print_coll_info_plist(enum list_type lt,atom_number p,atom_number q){
  alist *pt=plists[lt][p];
  while(pt->pnext){
    if(pt->q==q){print_list_type; printf(" ");  print_coll(pt); }
    pt=pt->pnext;
  }
}
/*========================================================================*/
void print_coll_info_qlist(enum list_type lt,atom_number p,atom_number q){
  alist *pt=qlists[lt][p];
  while(pt->qnext){
    if(pt->p==q){print_list_type; printf(" ");  print_coll(pt); }
    pt=pt->qnext;
  }
}
/*========================================================================*/
void print_coll_info_list(enum list_type lt, atom_number p, atom_number q){
  print_coll_info_plist(lt,p,q);
  print_coll_info_qlist(lt,p,q);
}
/*========================================================================*/
void print_coll_info(atom_number p, atom_number q){
  enum list_type lt;
  for(lt=0;lt<nlt;lt++){
    print_coll_info_list(lt,p,q);
  }
}
/*========================================================================*/
void print_plist_coll(enum list_type lt,atom_number p){
  alist *pt;
  if(p>=n1){ printf("ERROR: atom number too big!\n");exit(0);}
  print_list_type(lt);
  pt=plists[lt][p]; printf("plist:\n");
  while(pt->pnext){ print_coll(pt); pt=pt->pnext; }
}
/*========================================================================*/
void print_qlist_coll(enum list_type lt,atom_number q){
  alist *pt;
  if(q>=n1){ printf("ERROR: atom number too big!\n");exit(0);}
  print_list_type(lt);
  pt=qlists[lt][q]; printf("qlist:\n");
  while(pt->qnext){ print_coll(pt); pt=pt->qnext; }
}
/*========================================================================*/
void print_colls (enum list_type lt, atom_number p){
  alist *pt;
  if(p>=n1){ printf("ERROR: atom number too big!\n");exit(0);}
  print_list_type(lt);
  pt=plists[lt][p]; printf("plist:\n");
  while(pt->pnext){ print_coll(pt); pt=pt->pnext; }
  pt=qlists[lt][p]; printf("qlist:\n");
  while(pt->qnext){ print_coll(pt); pt=pt->qnext; }
}
/*========================================================================*/
void print_all_colls (atom_number p){
  enum list_type lt;
  if(p>=n1){ printf("ERROR: atom number too big!\n");exit(0);}
  for(lt=0;lt<nlt;lt++)
    print_colls(lt,p);
}
/*========================================================================*/
void check_plist_pointers (enum list_type lt, atom_number p){
  alist *pt=plists[lt][p],*pt2;
  char l[10];
  if(p>=n1){ printf("ERROR: atom number too big!\n");exit(0);}
  strcpy(l,print_char_list_type(lt));
  if(pt->pnext){
    if(pt->pprev!=NULL){
      printf("ERROR:first pt->pprev!=NULL for plists[%s][%d]\n",l,p);
      printf("pointed collision pt->pprev is:\n");
      print_coll(pt->pprev);
      printf("List of collisions plists[%s][%d]:\n",l,p);
      print_plist_coll(lt,p);
      exit(0);
    }
    while(pt->pnext){
      pt2=pt;
      pt=pt->pnext;
      if(pt2==pt){
	printf("ERROR:circular list for plists[%s][%d]\n",l,p);
	exit(0);
      }
    }
    if(pt2->pnext!=search.storage){
      printf("ERROR:last pt->pnext!=search.storage for plists[%s][%d]\n",l,p);
      printf("pointed collision pt->pnext is:\n");
      print_coll(pt2->pnext);
      printf("List of collisions plists[%s][%d]:\n",l,p);
      print_plist_coll(lt,p);
      exit(0);
    }
  }
}
/*========================================================================*/
void check_qlist_pointers (enum list_type lt, atom_number q){
  alist *pt=qlists[lt][q],*pt2;
  char l[10];
  if(pt->qnext){
    strcpy(l,print_char_list_type(lt));
    if(pt->qprev!=NULL){
      printf("ERROR:first pt->qprev!=NULL for qlists[%s][%d]\n",l,q);
      printf("pointed collision pt->qprev is:\n");
      print_coll(pt->qprev);
      printf("List of collisions qlists[%s][%d]:\n",l,q);
      print_qlist_coll(lt,q);
      exit(0);
    }
    while(pt->qnext){ 
      pt2=pt; 
      pt=pt->qnext;
      if(pt2==pt){
	printf("ERROR:circular list for qlists[%s][%d]\n",l,q);
	exit(0);
      }
    }
    if(pt2->qnext!=search.storage){
      printf("ERROR:last pt->qnext!=search.storage for qlists[%s][%d]\n",l,q);
      printf("pointed collision pt->qnext is:\n");
      print_coll(pt2->qnext);
      printf("List of collisions qlists[%s][%d]:\n",l,q);
      print_qlist_coll(lt,q);
      exit(0);
    }
  }
}
/*========================================================================*/
void check_list_pointers (enum list_type lt,atom_number p){
  if(p>=n1){ printf("ERROR: atom number too big!\n");exit(0);}
  check_plist_pointers(lt,p);
  check_qlist_pointers(lt,p);
}
/*========================================================================*/
void check_pointers (atom_number p){
  enum list_type ilt,lt,lta[3]={HOT,WARM,COLD};
  if(p>=n1){ printf("ERROR: atom number too big!\n");exit(0);}
  for(ilt=0;ilt<3;ilt++){ 
    lt=lta[ilt];
    check_list_pointers(lt,p); 
  }
}
/*========================================================================*/
void check_all_pointers (void){
  atom_number p;
  /*printf("check_all_pointers\n");*/
  for(p=0;p<n1;p++)
    check_pointers(p); 
  /*printf("End of check_all_pointers\n");*/
}
/*========================================================================*/
void renew_calendar (void){
  int i;
  /*printf("renew_calendar(...)\n");*/

  for(i=0;i<calendar.n_updates;i++){
    delete_node(calendar.update_node_for[i]);
    insert_node(calendar.update_node_for[i]);
  }

  find_new_champion();
  /*printf("End of renew_calendar(...)\n");*/
}
/*=========================================================================*/
/* Recipe to remove node "deleted" when has both left and right children:

      predecessor                                  predecessor
          |         /                    \              |
       deleted     | Recipe:             |           successor
        /  \       | (1)replace "deleted"|             / \
       /    \      |    by "successor"   |            /   \
      L      \     | (2)link "r" to "z"  |           L     \
     /\       R     \                    /          / \     R
    ?  ?     / \                                   ?   ?   / \   
            .   ?                                         .   ?
           .                                             .
          /                                             /
         r                                             r
        / \                                           / \
       /   \                                         z   ?
 successor  ?                                       / \    
    / \                                            ?   ? 
 NULL  z               
      / \
     ?   ?
*/
void delete_node (atom_number p){
  node *predecessor,*deleted,*successor;
  node *L,*R,*r,*z; /*see drawing above for meaning*/

  deleted=calendar.storage+p;
  predecessor=deleted->parent;

  if(deleted->slow_child==NULL){/*case I*/
    /*printf("CASE I...");*/
    successor=deleted->fast_child;/*"successor" could also be NULL */
  }
  else{
    if(deleted->fast_child==NULL){/*case II*/
      successor=deleted->slow_child;
    }
    else{
      L=deleted->fast_child;
      R=deleted->slow_child;
      successor=R->fast_child;
      if(successor==NULL){/*case III*/
	successor=R;
      }
      else{/*case IV, complicated situation depicted above*/
	while(successor->fast_child)
	  successor=successor->fast_child;
	r=successor->parent;/*init r*/
	z=successor->slow_child;/*init z, can be NULL*/
	r->fast_child=z;
	if(z) z->parent=r;
	successor->slow_child=R;
	R->parent=successor;
      }
      successor->fast_child=L;
      L->parent=successor;
    }
  }
  /*link now the successor to the predecessor*/
  if(predecessor){
    if(predecessor->fast_child==deleted)
      predecessor->fast_child=successor;
    else
      predecessor->slow_child=successor;
    if(successor)/*"successor" could also be NULL*/
      successor->parent=predecessor;
  }
  else{/*predecessor is NULL, thus "root" was pointing to "deleted"*/
    calendar.root=successor;
    successor->parent=NULL;
  }
  /*printf("End of delete_node(%d)\n",p);*/
}
/*========================================================================*/
void insert_node(atom_number p){
  node *new_node,*current,*prev;
  enum list_type wlt;
  double t;
  /*printf("insert_node(%d)\n",p);*/
  wlt=winlist[p];
  t=tlists[wlt][p]->t;
  new_node=calendar.storage+p;
  new_node->t=t;
  new_node->fast_child=NULL;
  new_node->slow_child=NULL;

  current=calendar.root;
  while(current){/*we assume there's allways at least one node in the tree*/
    prev=current;    
    if(t < current->t)
      current=current->fast_child;
    else
      current=current->slow_child;
  }
  new_node->parent=prev;
  if(t < prev->t)
    prev->fast_child=new_node;
  else
    prev->slow_child=new_node;
  /*printf("End of insert_node(..)\n");*/
}/*Matches void insert_node(atom_number p)*/
/*========================================================================*/
double d2 (atom_number p,atom_number q){
  double dr[3],half,dr2=0.0;
  int i;
  dimensions *bound=get_bounds();

  /*printf("\nfunction d2(...) ...\n");*/
  moveatom(a+p);
  moveatom(a+q);
  dr[0]=a[p].r.x-a[q].r.x;
  dr[1]=a[p].r.y-a[q].r.y;
  dr[2]=a[p].r.z-a[q].r.z; 
  for(i=0;i<3;i++){
    half=bound[i].length/2.0;
    if( dr[i]>half ) dr[i]-=bound[i].length;
    else if( dr[i]<-half ) dr[i]+=bound[i].length;
    dr2+=dr[i]*dr[i];
  }

  return dr2;  
}
/*========================================================================*/
void dump_coordinates(void){
  atom_number p;
  for(p=0;p<search.n;p++)
    printf("%d %5f %5f %5f\n",p+1,a[p].r.x,a[p].r.y,a[p].r.z);
}
/*======================================================================*/
void dump_calendar (void){
  atom_number p;
  node *pt;

  pt=calendar.storage;
  for(p=0;p<calendar.n1;p++){
    printf("p=%d t=%7.4f parent=%d fast_child=%d slow_child=%d\n",
	   p,pt->t,pt->parent,pt->fast_child,pt->slow_child);
    pt++;
  }
}
/*=====================================================================*/
void dump_calendarII(void){
  atom_number p;
  node *pt;

  pt=calendar.storage;
  for(p=0;p<calendar.n1;p++){
    if(pt->parent || pt==calendar.root){
      printf("p=%d t=%7.4f parent=%d fast_child=%d slow_child=%d\n",
	     p,pt->t,pt->parent,pt->fast_child,pt->slow_child);
    }
    pt++;
  }
}
/*=========================================================================*/
int check_cell_addresses(atom_number p){
  int i,j,k,address;
  clist *pt;
  dimensions * bound=get_bounds();
  
  i=a[p].i.x.i;   j=a[p].i.y.i;   k=a[p].i.z.i; 
  address=(k*bound[1].period+j)*bound[0].period+i;
  if(address != a[p].add){
    printf("address=%d, a[%d].add=%d\n",address,a[p].add);
    exit(0);
  }
  for(pt=search.cell_atoms[address];pt;pt=pt->next)
    if(pt->n==p) return 1;
  printf("%d is not in cell with address %d\n",p,address);
  printf("Details: a[%d].add=%d  a[%d].i.x.i=%d  a[%d].i.y.i=%d  a[%d].i.z.i=%d\n",p,a[p].add,p,i,p,j,p,k);
  printf("                              A=%d            B=%d            C=%d          \n",
	 bound[0].period,bound[1].period,bound[2].period);
  printf("Atoms in cell %d %d %d:",i,j,k); 
  for(pt=search.cell_atoms[address];pt;pt=pt->next)
    printf(" %d",pt->n);
  printf("\n");
  exit(0);
}

void chell_all_atom_cell_addresses(void){
  atom_number p;
  for(p=0;p<search.n;p++)
    check_cell_addresses(p);
}
/*=========================================================================*/
void dump_info_of (atom_number p){
  atom *ap=a+p;
  char l[10];
  strcpy(l,atom_name(ap->origc));
  printf("info for p=%d %s\n",p,l);
  printf("t=%18.15lf coord: %18.13lf %18.13lf %18.13lf\n",
	 ap->t,ap->r.x,ap->r.y,ap->r.z);
  printf("veloc: %lf %lf %lf\n",ap->v.x,ap->v.y,ap->v.z);
  printf("cell : %d %d %d\n",ap->i.x.i,ap->i.y.i,ap->i.z.i);
  printf("add: %d c: %d origc: %d\n",ap->add,ap->c,ap->origc);
  printf("sfr: %18.13lf %18.13lf %18.13lf\n",ap->sfr.x,ap->sfr.y,ap->sfr.z);
}
/*=========================================================================*/
void print_list_type(enum list_type lt){
  printf("list_type=");
  switch(lt){
  case HOT:  printf("HOT");  break;
  case WARM: printf("WARM"); break;
  case WALL: printf("WALL"); break;
  case SAFE: printf("SAFE"); break;
  case COLD: printf("COLD"); break;
  default: printf("ERROR: NO list_type found"); exit(0);
  }
  printf("\n");
}
/*=========================================================================*/ 
char *print_char_list_type(enum list_type lt){
  static char list[10];
  switch(lt){
  case HOT:  strcpy(list,"HOT");  break;
  case WARM: strcpy(list,"WARM"); break;
  case WALL: strcpy(list,"WALL"); break;
  case SAFE: strcpy(list,"SAFE"); break;
  case COLD: strcpy(list,"COLD"); break;
  default: printf("ERROR: NO list_type found"); exit(0);
  }
  return list;
}
/*=========================================================================*/
void check_inf_times(){
  atom_number p;
  alist *pt;
  enum list_type lt,lta[2]={HOT,WARM};/*COLD friends have t=INF*/

  for(p=0;p<n1;p++){
    for(lt=0;lt<2;lt++){

      pt=plists[lt][p];
      while(pt->pnext){
	if(pt->mso==INF_SO && pt->t>=DBL1){
	  printf("ERROR: infinite time for bounded neighbors!\n");
	  print_coll(pt); exit(0);
	}
	pt=pt->pnext;
      }
      /*Don't check qlists because plists sweeps all collisions*/
    }
    pt=plists[WALL][p];
    if(pt->t>=DBL1){
      printf("ERROR: infinite time for wall collision!\n");
      print_coll(pt); exit(0);
    }
    
    pt=plists[SAFE][p];
    if(pt->t>=DBL1){
      printf("ERROR: infinite time for safe crossing!\n");
      print_coll(pt); exit(0);
    }
  }
}/*Matches check_inf_times*/
/*=======================================================================*/
void check_zero_times(void){
  atom_number p;
  alist *pt;
  enum list_type lt,lta[2]={HOT,WARM};/*COLD friends have t=INF*/

  for(p=0;p<n1;p++){
    for(lt=0;lt<2;lt++){

      pt=plists[lt][p];
      while(pt->pnext){
	if(pt->t<ZERO){
	  printf("ERROR: zero time for coll!\n");
	  print_coll(pt); exit(0);
	}
	pt=pt->pnext;
      }
      /*Don't check qlists because plists sweeps all collisions*/
    }
    pt=plists[WALL][p];
    if(pt->t<ZERO){
      printf("ERROR: zero time for wall collision!\n");
      print_coll(pt); exit(0);
    }
    
    pt=plists[SAFE][p];
    if(pt->t<ZERO){
      printf("ERROR: zero time for safe crossing!\n");
      print_coll(pt); exit(0);
    }
  }
}
/*=======================================================================*/
int bond(atom_number i, atom_number j){   
  alist *pt;
  atom_number p,q;
  int k;
  enum list_type lt,lta[3]={HOT,WARM,COLD};

  if(less_than(i,j)){ p=i; q=j; }
  else{ p=j; q=i; }
  /*By finding the smaller of i and j, we only need to go through the plists*/
  for(k=0;k<3;k++){
    lt=lta[k];
    pt=plists[lt][p];
    while(pt->pnext){
      if(pt->q==q) return pt->c;
      pt=pt->pnext;
    }
  }
  /*i and j do not have any scheduled collision, thus return the most
    exterior potential well*/
  return -ecoll[a[i].origc][a[j].origc];
}
/*=======================================================================*/
