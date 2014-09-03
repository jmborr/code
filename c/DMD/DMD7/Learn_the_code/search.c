#include <math.h>
#include <stdio.h>
#include <float.h>
#include "bcp.h"
#include "search.h"
#include "controls.h"


static tsearch search; /*huge structure, defined in search.h*/
static atom * a;
static well_type ** ecoll;
static FILE * fp;

double maxrprobed;

void initsearch(void) /*47 lines of code*/
{ 
  alist * inpt;
  int i;
  int n=search.n;/*number of atoms*/
  alist * tc=search.storage;
  alist ** ptc=search.begin;

  for (i=0;i<n;i++)
    {
      search.collp[i]=-1;/*no atoms scheduled for collision calculation with atom "p"*/
      search.collq[i]=-1;
      search.atom_change[i]=0;
      search.atom_storage[i].n=i;
      search.atom_storage[i].next=NULL;
    }
  search.np=0;
  search.nq=0;
  search.ncha=0;

  /*begin members point to each element of storage*/
  for(i=0;i<search.maxfree;i++) ptc[i]=&tc[i];
  ptc[search.maxfree]=NULL;
  /*free points to this particular member of begin*/
  search.free=(search.begin)+2*search.n;
  for(i=0;i<2*search.n;i++)
    {/*initialize array search.storage (up to index 2*search.n-1), to which search.begin points to*/
      inpt=search.begin[i];
      inpt->t=DBL2; /*bcp.h:10:#define DBL2 (DBL_MAX*10e-10)*/
      inpt->p=-1;
      inpt->q=-1; 
      inpt->c=-1;
      inpt->pprev=NULL;
      inpt->qprev=NULL;
      inpt->pnext=NULL;
      inpt->qnext=NULL;
    }

  /*atom_collisionp points to first n1 members of storage (because  begin is pointing to storage. Analogously, atom_collisionq  points to second n1 members of storage*/
  for(i=0;i<n;i++)
    {
      search.atom_collisionp[i]=search.begin[i];
      search.atom_collisionq[i]=search.begin[i+search.n];
    }

  for(i=0;i<=search.z;i++)/*search.z=maxadd*/
    {
      search.cell_atoms[i]=NULL;
      search.cell_champion[i]=NULL;
    } 
}/*Matches void initsearch(void)*/


int realloc_search (void)
{
  int i,j,k,level;
  int maxadd;
  dimensions * bound=get_bounds();
  if ((search.x==bound[0].period)&&(search.y==search.x*bound[1].period)
      &&(search.z==search.y*bound[2].period))return 1;
  search.x=bound[0].period;
  search.y=search.x*bound[1].period;
  if(search.z==search.y*bound[2].period)return 1;
  search.z=search.y*bound[2].period;
  maxadd=search.z;
   k=0;
  j=maxadd;
  level=0;
  do
    { j=(j>>1)+(j&1);
      k+=j;
      level++;
    }while(j>1);
  search.final=level;

if(maxadd>search.maxalloc)
  {
    free(search.cell_atoms);
    free(search.cell_champion);
    if(!(search.cell_atoms=(clist **)malloc((maxadd+1)*sizeof(clist *))))return 0;
    if(!(search.cell_champion=(alist **)malloc((maxadd+1)*sizeof(alist*))))return 0;
    free(search.olymp[0]);
    free(search.olymp);
    search.olymp=(size_al **)malloc((level+1)*sizeof(size_al *));
    if(!(search.olymp))return 0;
    search.olymp[0]=(size_al *)malloc((k+1)*sizeof(size_al));
    if(!(search.olymp[0]))return 0;
    j=maxadd;
    level=0;
    
    while(j>1)
      {
	j=(j>>1)+(j&1);
	search.olymp[level+1]=search.olymp[level]+j;
	level++;
      }
    
    search.change=(size_al *)malloc(46*sizeof(size_al));
    if(!(search.change))return 0;
    search.maxalloc=maxadd; 
  }
  return 1;
}

/*When invoked from make_tables() by the line allocsearch(n1),
we are passing the number of atoms n1 to local variable n. This
is confusing because there is a global variable n used in 
make_system.c which has the value (n1-1) */
int allocsearch(int n)
{
  int i,j,k,level;
  int maxadd; /*biggest possible address, or total number of cells.
The adress specifies the cell index in only one integer*/
  dimensions * bound=get_bounds(); /*returns pointer to global
				     variable named also bound*/
  a=get_atom(); /*returns pointer to global variable named also
"a", which was previously initialized in make_system.c*/
  ecoll=get_ecoll(); /*returns pointer to global variable named
		       also ecoll*/
  search.x=bound[0].period; /*number of cells along x-axis*/
  search.y=search.x*bound[1].period; /*number of cells on a XY plane*/
  search.z=search.y*bound[2].period;  /*total number of cells*/
  maxadd=search.z; /*the maximum address is the total number of cells*/
  search.maxalloc=search.z;
  search.n=n;
  printf("search.n=%d\n",search.n);
  search.nhalf =n>>1;     
  if(maxadd>MAXADD)return 0; /*ERROR: too many cells for computer memory*/
  
  if(!(search.collp=(int *)malloc(n*sizeof(int))))return 0;
  if(!(search.collq=(int *)malloc(n*sizeof(int))))return 0;
  
  if(!(search.atomp=(int *)malloc(n*sizeof(int))))return 0;    
  if(!(search.atomq=(int *)malloc(n*sizeof(int))))return 0;
  
  if(!(search.atom_change=(int *)malloc(n*sizeof(int))))return 0;    
  if(!(search.changed_atoms=(int *)malloc(n*sizeof(int))))return 0;
  
  if(!(search.cell_atoms=(clist **)malloc((maxadd+1)*sizeof(clist *))))return 0;
  if(!(search.cell_champion=(alist **)malloc((maxadd+1)*sizeof(alist*))))return 0;
  if(!(search.atom_storage=(clist *)malloc(n*sizeof(clist))))return 0;
  
  search.maxfree=NFREE; /*bcp.h:1:#define NFREE (1000000)*/
  if(!(search.storage =(alist *)malloc(search.maxfree*sizeof(alist))))return 0;
  if(!(search.begin=(alist **)malloc((search.maxfree+1)*sizeof(alist *))))return 0;
  if(!(search.atom_collisionp=(alist **)malloc((n)*sizeof(alist *))))return 0;
  if(!(search.atom_collisionq=(alist **)malloc((n)*sizeof(alist *))))return 0;
                        /*                      #                        */    
  k=0;                  /*         #------------|-----------#            */
  j=maxadd;             /*         #                 #------|-------#    */ 
  level=0;              /* #-------|--------#  #-----|----#   #-----|---#*/      
                        /*                                               */
  do                    /*The previous diagram shows a binary tree with  */
    { j=(j>>1)+(j&1);   /*four levels. The base is level=0 and the top is*/
    k+=j;               /*level=3, which later is passed to search.final */
    level++;            /*Thus we construct a binary tree where the base */
    }while(j>1);        /*begins with a number of nodes equal to the     */ 
                        /*number of cells                                */
  search.final=level;
  search.olymp=(size_al **)malloc((level+1)*sizeof(size_al *));
  if(!(search.olymp))return 0; 
  search.olymp[0]=(size_al *)malloc((k+1)*sizeof(size_al));/*k is the number of nodes in the tree, excluding the base of the tree. Thus, search.olymp[0][0] to search.olymp[0][maxadd>>1 + maxadd&1] will store all the nodes of the first level, not the ground level*/
  if(!(search.olymp[0]))return 0;
  j=maxadd;
  level=0;
  
  while(j>1)
    {
      j=(j>>1)+(j&1);
      search.olymp[level+1]=search.olymp[level]+j;/*search.olymp[x] will point to the nodes of level x+1*/
      level++;
    }
  /*Note that we have only allocated search.olymp, but not filled it*/
  search.change=(size_al *)malloc(46*sizeof(size_al));/*46 is the maximum number of cells affected by the collision of two atoms. The way in which the cell is affected is that we erase collisions that belong to that particular cell. If "nch" cells are affected during a collision, then search.change[0] to search.change[nch-1] store the adresses of the affected cells.*/
  if(!(search.change))return 0;
  return 1;
}

/* used in find_atoms when atom "p1" moves from one cell to another */
int find_collisions(int p1)/*20 lines of code*/
{
  int np=0;
  alist * pt=search.atom_collisionp[p1];/*point to the first collision in the "p" collision list for atom "p1". This collision is nothing but the crossing of the cell boundary, since it is the champion collision*/
  *(--search.free)=pt;/*claim the space occupied by this collision as free, i.e., we remove this collision from the list*/
  pt=pt->pnext;/*shift the first collision to happen in the "p" collision list for atom "p1" to the next collision*/
  pt->pprev=NULL;/*in a free space all pointers must point to null*/
  search.atom_collisionp[p1]=pt;/*now atom_collisionp[p1] will point to the collision which previously was the second collision in the "p" collision list for atom "p1"*/
  /*Anotate that "p1" has suffered a "collision". This only a effects of calculating the new champion collision*/
  if(!search.atom_change[p1])
    {
      search.atom_change[p1]=1;
      search.changed_atoms[search.ncha++]=p1;
    }
  /*Now we go through the "p" collision list for atom "p1". We do not erase these collisions but we anotate who the atom partners are, so that we do not compute the collision time again*/
  while(pt->p==p1)
    {
      int q=pt->q;
      if(search.collp[q]==-1)search.atomp[np++]=q;
      search.collp[q]=-2;
      pt=pt->pnext;/*go through the next collision in the p-list of p1*/
    }
  pt=search.atom_collisionq[p1];
  /*Do the same now for the "p" collision list for atom "p1"*/
  while(pt->q==p1)
    {
      int q=pt->p;
      if(search.collp[q]==-1)search.atomp[np++]=q;
      search.collp[q]=-2;
      pt=pt->qnext;
    }
  return np;
}/*Matches int find_collisions(int p1)*/

/* kill_collisions called only by squeeze1 */ 
int kill_collisions(int p1, int * collp, int * atomp)/*72 lines of code*/
{
  int np=0;
  alist ** free=search.free;
  alist * pt=search.atom_collisionp[p1];/*points to first scheduled collision between p1 and other atom whose atom index is in [p1+1,p1+n/2] (Note that if p1+n/2>n, then p1+n/2 is p1-n/2, because we circularize the list so that index n+1==0)*/
  int ncha=search.ncha;/*in this point of the code, search.ncha==0*/
  int *changed_atoms=search.changed_atoms;/*array of size n1*/
  int *atom_change=search.atom_change;/*array of size n1. For a particular index "i", atom_change[i]==1 indicates whether this atom has collided with some other atom. Otherwise, atom_change[i]==0*/
  if(!atom_change[p1])
    {
      atom_change[p1]=1;
      changed_atoms[ncha++]=p1;/*list of atoms that have collided*/
    } 
  while(pt->p==p1)/*go through the p-list of p1*/
    {
      alist * qprev=pt->qprev;/*Here we will take care that the q-list of q (not the q-list of p1) is arranged correctly after we kill this collision*/
      int q=pt->q;/*collision partner of p1*/
      if(q<search.n)/*not a cell boundary crossing*/
	{
	  if(collp[q]==-1)/*-1 is initialezed value*/
	    {
	      collp[q]=pt->c;/*store the type of the scheduled collision (that is going to be erased) into collp[q]*/
	      atomp[np++]=q;/*store the atom index of the partner of p1 into atomp, and increase the number of putative partners of p1*/
	    }
	  if(qprev==NULL)/*if this collision wast the last in the q-list of q, then search.atom_collisionq[q] was actually pointing to it. Now we have to change this, so that search.atom_collisionq[q] will point to the next collision in the q-list of q because we are going to erase this collision*/
	    {
	      search.atom_collisionq[q]=pt->qnext;
	      search.atom_collisionq[q]->qprev=NULL;/*Now the next collision is going to be the last collision in the q-list*/
	    }
	  else/*if this collision is not the last i the q-list of q, but stays in the middle, then we have to jump over this collision in the q-list because we are going to erase*/
	    { 
	      qprev->qnext=pt->qnext;
	      pt->qnext->qprev=qprev;
	    }
	}
      *(--free)=pt;/*reclaim the space occupied by this collision as free space*/
      pt=pt->pnext;
      pt->pprev=NULL;
    }
  search.atom_collisionp[p1]=pt;/*the p-list of p1 has shrinked now to the initial fake-collision that signals the end of the p-list, so that now search.atom_collisionp[p1] to this fake-collision*/

  pt=search.atom_collisionq[p1];/*go through the q-list of p1*/
  while(pt->q==p1)
    {
      alist * pprev=pt->pprev;
      int p=pt->p;
      if(collp[p]==-1)
	{
	  collp[p]=pt->c;
          atomp[np++]=p;
	}
	  if(pprev==NULL)
	    {
	      search.atom_collisionp[p]=pt->pnext;
	      search.atom_collisionp[p]->pprev=NULL;
	      if(!atom_change[p])
		{
		  atom_change[p]=1;
		  changed_atoms[ncha++]=p;
		} 
	    }
	  else
	    { 
	      pprev->pnext=pt->pnext;
	      pt->pnext->pprev=pprev;
	    }
      *(--free)=pt;
      pt=pt->qnext;
      pt->qprev=NULL;
    }
  search.atom_collisionq[p1]=pt;
  search.ncha=ncha;
  search.free=free;
  return np;
}/*Matches int kill_collisions(int p1, int * collp, int * atomp)*/

/* used in new_loc in bcp.c */
void find_atoms(int p1, size_al address1)/*77 lines of code*/
{ 
  clist * inpt;
  int np;
  int nch=search.nch;
  int i,j,k,i1,j1,k1,i2,j2,k2,address=a[p1].add;/*address of the cell before  atom "p1" crossed the cell boundary (old cell). In contrast, address1 is the addres of the cell after the atom crossed the cell boundary (new cell)*/
  int addressz,addressy;
  clist * old_rec=search.cell_atoms[address];/*pointer to the list of atoms in the old cell. We have to remove the entry "p1" from this list, as done below*/
  if(old_rec->n==p1)/*if the first atom number entry in the list of atoms for the old cell is actually "p1", then...*/
    search.cell_atoms[address]=search.cell_atoms[address]->next;
  else/*go through the list until p1 is found*/
    {
      clist *  prev_rec=old_rec;
      old_rec=prev_rec->next;
      while(old_rec->n!=p1)
	{
	  prev_rec=old_rec;
	  old_rec=prev_rec->next;
	}
      prev_rec->next=old_rec->next;
    }
  old_rec->next=search.cell_atoms[address1];/*now the first atom entry in the list of atoms for the new cell is going to be "p1"*/
  search.cell_atoms[address1]=old_rec;
  a[p1].add=address1;/*actualize the cell address of atom structure for p1*/
  search.cell_champion[address]=NULL;/*since "p1" colliding with the wall was the cell champion in the old cell, now we need a new champion*/
  search.change[nch]=address;
  nch++;
  search.cell_champion[address1]=NULL;/*it may be now that p1 colliding with some other atom is going to be the cell champion for the new cell*/
  search.change[nch]=address1;
  nch++;
  search.nch=nch;
  np=find_collisions(p1);  

  /*go through all 27 cells looking for partners atoms to collide. From "find_colisions", we already know previous atoms that were to collide. It is strange to go over the 27 cells, since new collisions only can happen in the 9 new cells that are probed when p1 crosses a cell boundary*/
  i1=a[p1].i.x.i-1;
  j1=(a[p1].i.y.i-1)*search.x;
  i2=i1+2;
  j2=j1+(search.x<<1);
  if(search.z==search.y)
  {k1=a[p1].i.z.i*search.y;k2=k1;}
  else
  {k1=(a[p1].i.z.i-1)*search.y;k2=k1+(search.y<<1);}
  
  for(k=k1;k<=k2;k+=search.y)
    { 
      addressz=k;
      if(addressz<0)addressz+=search.z;
      if(addressz==search.z)addressz=0;
      for(j=j1;j<=j2;j+=search.x)
	{ 
	  addressy=j; 
	  if(addressy<0)addressy+=search.y;
	  if(addressy==search.y)addressy=0;
	  addressy+=addressz;
	  for(i=i1;i<=i2;i++)
	    { 
	      address=i; 
	      if(address<0)address+=search.x;
	      if(address==search.x)address=0;
	      address+=addressy;
	      /*go through the list of atoms contained in cell with "address"*/
	      for(inpt=search.cell_atoms[address];inpt;inpt=inpt->next)
	      {
		int p2=inpt->n;
		/*new possible collisions to be scheduled if search.collp[p2]==-1 and if the atom is not p1. This last event will happen when we are looking thorgh the list of atoms contained in the same cell that the cell where p1 is located. Old collisions have search.collp[p2]==-2*/
		if((search.collp[p2]==-1)&&(p1!=p2))
		  { 
		    search.atomp[np++]=p2;
                    search.collp[p2]=ecoll[a[p1].c][a[p2].c];/*the collision type has to be the outmost type, and cannot be with an atom that is a link because the it would have already an scheduled collision*/
		  }

	      }
	    } 

	}
    }  

 search.nq=0;
 search.np=np;
}/*Matches void find_atoms(int p1, size_al address1)*/

/* find atoms in ajacent cells to a given atom; 
   returns number of such atoms */
int find_neighbors(int p, int * neib)
{ 
 clist * inpt;
 int nn=0;
 int i,j,k,i1,j1,k1,i2,j2,k2,address;
 int addressz,addressy;
  
  i1=a[p].i.x.i-1;
  j1=(a[p].i.y.i-1)*search.x;
  i2=i1+2;
  j2=j1+(search.x<<1);
  if(search.z==search.y)
  {k1=a[p].i.z.i*search.y;k2=k1;}
  else
  {k1=(a[p].i.z.i-1)*search.y;k2=k1+(search.y<<1);}
  
  for(k=k1;k<=k2;k+=search.y)
    { 
      addressz=k;
      if(addressz<0)addressz+=search.z;
      if(addressz==search.z)addressz=0;
      for(j=j1;j<=j2;j+=search.x)
	{ 
	  addressy=j; 
	  if(addressy<0)addressy+=search.y;
	  if(addressy==search.y)addressy=0;
	  addressy+=addressz;
	  for(i=i1;i<=i2;i++)
	    { 
	      address=i; 
	      if(address<0)address+=search.x;
	      if(address==search.x)address=0;
	      address+=addressy;
	      for(inpt=search.cell_atoms[address];inpt;inpt=inpt->next)
		    neib[nn++]=inpt->n;
	    } 
	}  
    }
 return nn;
}  


/*return number of atoms in the cell specified by address, and fill
array atomx with the atom numbers of the atoms belonging to the cell*/
int list_atoms(size_al address,int * atomx)
{
	int i=0;
	clist *inpt1;
	for( inpt1=search.cell_atoms[address]; inpt1; inpt1=inpt1->next)
 	{ atomx[i++]=inpt1->n;}
 return i;	
}

int bond(int i, int j)
{   
    alist *inpt1;
	for( inpt1=search.atom_collisionp[i];inpt1->p==i;inpt1=inpt1->pnext)
	if(j==inpt1->q) return inpt1->c;
	for( inpt1=search.atom_collisionq[i];inpt1->q==i;inpt1=inpt1->qnext)
	if(j==inpt1->p) return inpt1->c;
	return -ecoll[a[i].c][a[j].c];
}

/*will link search.atom_storage[i].next to the last search.atom_storage[?] of
the list that belongs to the same cell as search.atom_storage[i]*/
void add_atom2cell(int i)
{
  size_al address=a[i].add;
  clist * pt=search.cell_atoms[address]; 
  search.cell_atoms[address]=&search.atom_storage[i];
  search.cell_atoms[address]->next=pt;
}

alist * local_champion(size_al address)
{
  double t;
  int i;
  clist *  inpt=search.cell_atoms[address];/*remember search.cell_atoms[address] points to the beginning of list of atoms that are contained in cell specified by "address"*/
  if(!inpt)return *search.begin;
  i=inpt->n;/*first atom index of the list*/
  t=search.atom_collisionp[i]->t;/*obtain time of next scheduled collision */
                          /*for atom i in its p-list of scheduled collisions */
  inpt=inpt->next;
  while(inpt)/*go through the list of atoms contained in the cell to find the*/
             /*collision with the minimal time abscribed to this cell*/
    {
      int j=inpt->n;
      double t1=search.atom_collisionp[j]->t;
      if(t1<t)
	{
	  t=t1;
	  i=j;
	}
      inpt=inpt->next;
    }
  return search.atom_collisionp[i];/*we return the address in memory of the */
                                   /*scheduled collision*/
}/*Matches alist * local_champion(size_al address)*/

int init_tables(void)/*204 lines of code*/
{
  int i0,j0,k0,address,ix,iy,iz,level;
  int i,j,k,i1,j1,k1,i2,j2,k2;
  int maxadd=search.z;
  size_al address1; /*typedef unsigned int*/
  size_al address2;
  int n=search.n-1;/*since search.n is number of atoms, this n represents the maximum index, like the n used in bcp.c and make_system.c*/
  int n1=n+1;/*same idea as n1 in bcp.c and make_system.c*/
  double t;
  int p,q,ct;
  
  int addressz,addressy; 
  
  initsearch();

  for(i=0;i<=n;i++)/*this loop schedules collisions of all atoms with the    */
                   /*walls, and insert the collisions in the collisions lists*/
    {
      add_atom2cell(i);/*arrange arrays cell_atoms and atom_storage,so that  */
                       /*cell_atoms[ a[i].add ]  point to atom_storage[i]    */
      q=twall(i,&t);   /*q=n1||n2||n3, calculates when atom will collide with*/
                       /*wall, and with which wall                           */
      bubble(t,i,q,-1);/*-1 indicates a collision type with a wall, and      */
                       /*allocate the collision in array search.storage in   */
                       /*the list of collisions of atom "i"                  */
    }
  
  a[n1].c=0; /*remember atom *a is an array of size number of atoms plus one,*/
             /*thus n1 is precisely the last index                           */
  for (k0=0;k0<maxadd;k0++)/*do for all cells in the system                  */
    if(search.np=list_atoms(k0,search.atomp))/*return number of atoms in cell*/
        /*k0 and fill atomp with the atom numbers of the atoms inside cell k0*/
      {
	i0=search.atomp[0]; /*first atom of the list                         */
	i1 = a[i0].i.x.i -1;/*adress of the cell adjacent to cell k0 along   */
	                    /*the x-axis to the left                         */
	i2=i1+2;            /*address of the cell adjacent to cell k0 along  */
                            /*the x-axis to the right                        */
	j1=(a[i0].i.y.i-1)*search.x;/*same for y-axis                        */
	j2=j1+(search.x<<1);
	if(search.z==search.y)/*if bidimensional system                      */
	  {
	    k1= a[i0].i.z.i * search.y;/*total adress                        */
	    k2=k1;
	  }
	else
	  {
	    k1= (a[i0].i.z.i-1) * search.y;
	    k2= k1 + (search.y<<1);
	  }
	for(k=k1;k<=k2;k+=search.y)/*Find collisions inside the 27 cells     */
	  {                       
	    addressz=k;
	    if(addressz<0)addressz+=search.z;/*Implement the periodic        */
	    if(addressz==search.z)addressz=0;/*boundary conditions on        */
	    for(j=j1;j<=j2;j+=search.x)      /*the z-axis                    */
	      { 
		addressy=j; 
		if(addressy<0)addressy+=search.y;/*implement the periodic    */
		if(addressy==search.y)addressy=0;/*boundary conditions on    */
		addressy+=addressz;              /*the y-axis                */  
		for(i=i1;i<=i2;i++)
		  { 
		    int iq;
		    int ip;
		    address=i; 
		    if(address<0)address+=search.x;/*implement the periodic  */
		    if(address==search.x)address=0;/*boundary conditions on  */
		    address+=addressy;             /*the x-axis              */
   /*Now obtain list of atoms in cell specified by adress, and store in atomq*/
		    search.nq=list_atoms(address,search.atomq);
		    for (iq=0;iq<search.nq;iq++)  /*Find collisions between  */
		      for (ip=0;ip<search.np;ip++)/*atoms on the two adjacent*/
			{                         /*cells                    */
			  i0=search.atomp[ip];
			  j0=search.atomq[iq];
			  if(i0<j0)/*only store the collision if i0<j0,      */
			    {      /*otherwise, wait till we take the adja-  */
			           /*cent cell as the central one. Then i0   */
			           /*will be j0 and viceversa                */
			      int ct=collision_type(i0,j0);/*remember ct is  */
    			           /*the index of array coll, which list all */
			           /*ellastic/non-ellastic collisions        */
			      if(ct<0) return (int) ct;
			      add_potential(ct);
			      if(tball(i0,j0,ct,&t))/*Calculate collision,   */
				{                 /*returns 0 is no collision*/
				  p=i0;
				  q=j0;
				}
			      else/*fake collision, will not be inserted by  */
				{ /*bubble*/
				  p=-1;
				  q=-1;
				  ct=-1;
				}  
			      bubble(t,p,q,ct);
			    }
			}
		  }
	      }
	  }
      }
  
  for (k0=0;k0<maxadd;k0++)
    search.cell_champion[k0]=local_champion(k0);/*local_champion(k0) returns */
              /*the adress in memory where is located the scheduled collision*/
              /*with the minimal time of all the collisions abscribed to cell*/
              /* specified by adress "k0"                                    */
              /*Now we proceed to fill all search.olymp with the adresses of */
              /*the cells whose cell-times are winning                       */
  k=maxadd;
  k=(k>>1)+(k&1);/*k=(int)(maxadd/2) if k even, k=(int)(maxadd/2)+1 if k odd */
  for(i=0;i<k;i++)
    {
      address1=i<<1; /*add 1 to "i", unless i==0                             */
      address2=address1+1;/*address1 and address2 are two cells with         */
                     /*consecutive cell addresses                            */
      search.olymp[0][i]=address1;/*Remember olymp[0][] is not the ground    */
                     /*level, but the first level, ie. [i] goes from [0] to  */
                     /*[(maxadd>>1)+(maxadd&1)-1]. In this line we assume    */
                     /*first that cell address1 contains smallest time       */
      if((address2<maxadd)&&(search.cell_champion[address1]->t
			     >search.cell_champion[address2]->t))
                     /*maxadd-1 is the address of the last cell.Here we check*/
	             /*if address2 has a smaller time than address1          */
	search.olymp[0][i]=address2;
    }
  
  for(level=1;level<search.final;level++)/*Once the first level is filled, we*/
    {                /*proceed with the binary tree and fill the rest of the */
      j=(k>>1)+(k&1);/*levels                                                */
      for(i=0;i<j;i++)
	{
	  size_al i1=i<<1;
	  size_al i2=i1+1;
	  search.olymp[level][i]=search.olymp[level-1][i1];/*Assume first    */
	             /*search.olymp[level-1][i1] has the smallest time       */
	  if(i2<k)
	    {
	      address1=search.olymp[level-1][i1];
	      address2=search.olymp[level-1][i2];
	      if(search.cell_champion[address1]->t 
		 > search.cell_champion[address2]->t)/*check if [i2] has a smaller time than [i1]*/
		search.olymp[level][i]=address2;
	    }
	}
      k=j;
    }
  return 1;  
}/*matches int init_tables(void)*/

/*bubble assumes that collision to be inserted in the collision list is not */
/*already there                                                             */
int bubble (double t, int p1, int q1, int c)
{
  int p,q;/*next we will require that p<q*/
  alist *npt,*pt,*inpt; /*will point to different cells of search.storage*/
  /*printf("buble called, t=%lf p1=%d q1=%d c=%d\n",t,p1,q1,c) ;*/
  if(q1<p1)
    {
      p=q1;
      q=p1;
    }
  else      /*Invert so that now p>q! The reason is that if we don't do this */
    {       /*inversion, then for instance all scheduled collisions involving*/
      p=p1; /*atom 0 will be stored in the p-list of atom 0,and on the other */
      q=q1; /*hand,no collision involving atom n will be stored in its p-list*/
    }       /*To avoid this imbalance,we "circularize" the atoms so that atom*/
            /*n+1 will be actually atom 0. Now let's put an example. If we   */
            /*have a system with ten atoms, then all scheduled collisions    */
            /*involving atom 0 and atoms with indexes 1,2,3,4 will be stored */
            /*in the p-list of atom 0.All scheduled collisions involving atom*/
            /* 2 and atoms 3,4,5,6 will be stored in the p-list of atom 2.   */
            /*All scheduled collisions involving atom 8 and atoms 9,0,1,2will*/
            /*be stored in the p-list of atom 8, hence the circularization   */
  if((q<search.n)&&(q-p>search.nhalf))/*search.n is number of atoms, i.e q   */
    {       /*is not a wall index. Also, nhalf = n>>1                        */
      int r=q;
      q=p;
      p=r; 
    }	   
  if (p>=0)
    {	   
      npt=*(search.free);/*search.free initially points to                   */
                         /*search.begin[2*search.n], which in turn points to */
                         /*search.storage[2*search.n]. Thus *(search.free) = */
                         /*= &storage[2*search.n];                           */
      if(!npt){writetext(get_text_name());exit(0);}
      search.free++;/*We are goint to fill the element of array              */
      npt->t=t;     /*search.storage where *(search.free) was pointing to.   */
      npt->p=p;     /*Since it is now occupied, we have to shift the pointer */
      npt->q=q;     /*to a non-filled space in memory                        */
      npt->c=c;
      if(q<search.n)/*Walls are not atoms, thus no q-list for them           */
	{
	  inpt=search.atom_collisionq[q];
	  npt->qnext=inpt;
	  inpt->qprev=npt;          
	  npt->qprev=NULL;
	  search.atom_collisionq[q]=npt;
	}
      inpt=search.atom_collisionp[p];/*initially points to search.storage[p]*/
      if(inpt->t >= t)               /*Initially inpt->t=DBL2               */
	{
	  npt->pnext=inpt;
	  inpt->pprev=npt;          
	  npt->pprev=NULL;
	  search.atom_collisionp[p]=npt;
	  return p;
	}
      else 
	{
	  pt=inpt->pnext;
	  while ( pt->t<t)
	    pt=pt->pnext;
	  
         pt->pprev->pnext=npt;
         npt->pprev=pt->pprev;
         npt->pnext=pt;
         pt->pprev=npt;
	}
    }
    return -1;
}


void olymp_sort(void)/*78 lines of code*/
{/*As opposed to init_tables when we filled all the binary tree, here we only modify it with the cells that have been affected by the collision*/
  size_al *change=search.change;/*typedef unsigned int size_al. search.change is an array of 46 elements, which is the maximum number of cells affected by the collision of two atoms. The way in which the cell is affected is that we erase collisions that belong to that particular cell. If "nch" cells are affected during a collision, then search.change[0] to search.change[nch-1] store the adresses of the affected cells.*/
  size_al *curr, *next, *prev;
  size_al i,j,k,address1,address2;
  int ch,nch1,nch=search.nch;
  int level;
  nch1=0;/*number of current processed cells*/
  curr=search.olymp[0];/*pointer to the first level (not the ground) of the tree*/
  next=search.olymp[1];/*pointer to the next level in the tree to which curr is pointing*/
  for(ch=0;ch<nch;ch++)/*for each cell that is affected, process the loop*/
    {
      i=change[ch];/*i contains the address of one of the affected cells*/
      search.cell_champion[i]=local_champion(i);/*local_champion(i) returns the adress in memory where is located the scheduled collision with the minimal time of all the collisions abscribed to cell specified by adress "i"*/
      address1=i>>1; /*divide by 2 and take the integer part.*/    
      if(curr[address1]<MAXADD)
	{
	  curr[address1]=MAXADD;/*label curr[address1]=search.olymp[0][address1] with MAXADD, to denote that */
	  change[nch1]=address1;/*store the indexes of the first level to be tested */
	  nch1++;/*number of indexes of the first level that are to be tested */
	}
    }
  nch=nch1;/*nch is not the number of first level indexes that have been tested*/
  nch1=0;/*number of next (second) level that are processed up to now*/

  for(ch=0;ch<nch;ch++)
    {
      i=change[ch];
      address1=i<<1;/*obtain the index of the ground level, which is the cell address*/
      address2=address1+1;
      curr[i]=address1;
      if((address2<search.z)&&(search.cell_champion[address1]->t>search.cell_champion[address2]->t))/*compare the two cell. At least, one of the two cells correspond to an affected cell because of the collision*/
	curr[i]=address2;
      j=i>>1;/*divide by two. Correspond to the index in the next level*/
      if(next[j]<MAXADD)
	{
	  next[j]=MAXADD;/*MAXADD indicates to be tested*/
	  change[nch1]=j;
	  nch1++;
	} 
    }
  nch=nch1;
  k=search.z; 
  for(level=2;level<search.final;level++) /**/
    { 
      k=(k>>1)+(k&1);
      nch1=0;
      prev=curr;
      curr=next;
      next=search.olymp[level];
      for(ch=0;ch<nch;ch++)
	{ 
	  i=change[ch];
	  address1=i<<1;
	  address2=address1+1;
	  address1=prev[address1];
	  curr[i]=address1;
	  if(address2<k)
	    { 
	      address2=prev[address2];
	      if(search.cell_champion[address1]->t>search.cell_champion[address2]->t)
		curr[i]=address2;
	    }
	  j=i>>1;
	  if(next[j]<MAXADD)
	    {
	      next[j]=MAXADD;
	      change[nch1]=j;
	      nch1++;
	    } 
	  
	}  
      nch=nch1;
    }
  next[0]=(search.cell_champion[curr[0]]->t>search.cell_champion[curr[1]]->t) ? curr[1]:curr[0];/*The champion of the system!*/
  return;  
}/*matches void olymp_sort(void)*/

int get_free(void)
{
return (int)(search.free-search.begin);
}
int get_maxfree(void)
{
return search.maxfree;
}
void set_maxfree(int a)
{
search.maxfree=a;
}


void update_table(int p1, int q1, int ct1)/*101 lines of code*/
{
  int i,p2,ct;
  int add;
  int ncha=search.ncha;/*numbers of atoms that have collided*/
  double t;
  int * changed=search.changed_atoms;/*list of atoms that have collided*/
  int * change=search.atom_change;
  int p=p1;
  int q=twall(p1,&t);/*first schedule collision with wall*/

  if((p=bubble(t,p,q,-1))>=0)
    { 
      if(!change[p])
	{ 
	  change[p]=1;
	  changed[ncha++]=p;
	}
    }
  search.collp[p1]=-1;

  if (q1<search.n)
    {
      search.collp[q1]=ct1; /*to make sure that we compute the collision of atoms p1 and q1, when we will look through atomq which is different from the old type ct1 in the list collq */
      search.collq[p1]=-1;
      search.collq[q1]=-1;
      p=q1;
      ct=-1;
      q=twall(q1,&t);
      if((p=bubble(t,p,q,-1))>=0)
	{ 
	  if(!change[p])
	    { 
	      change[p]=1;
	      changed[ncha++]=p;
	    }
	}   
    }

  for(i=0;i<search.np;i++)
    {
      p2=search.atomp[i];
      ct=search.collp[p2];
      search.collp[p2]=-1;
      if(ct>=0) /* the collision is not recalculated if ct is -2 which means that the atoms will collide as before */
	if(tball(p2,p1,ct,&t))
	  {  
	    if((p=bubble(t,p1,p2,ct))>=0)
              { 
		if(!change[p])
		  { 
		    change[p]=1;
		    changed[ncha++]=p;
		  }   
	      }
	  }
    }
  if (q1<search.n)
    for(i=0;i<search.nq;i++)
      {
	p2=search.atomq[i];
	ct=search.collq[p2];
	search.collq[p2]=-1; 
	if(ct>=0)
	  if(tball(p2,q1,ct,&t))
	    {
	      if((p=bubble(t,q1,p2,ct))>=0)
		{ 
		  if(!change[p])
		    { 
		      change[p]=1;
		      changed[ncha++]=p;
		    }  
		}   
	    }
      }
  for(i=0;i<ncha;i++)
    {
      int j=changed[i];
      if(change[j])
	{
	  change[j]=0;
	  if(search.cell_champion[a[j].add])
	    {
	      search.change[search.nch++]=a[j].add;
	      search.cell_champion[a[j].add]=NULL;
	    }
	}
    }
  search.ncha=0;
  olymp_sort();
}/*Matches void update_table(int p1, int q1, int ct1)*/


/* squeeze1 called only by squeeze_table */
int squeeze1(int p1, int * collp, int * atomp)/*59 lines of code*/
{ 
 clist *inpt;
 int i,j,k,i1,j1,k1,i2,j2,k2,address;
 int addressz,addressy;
 int np=kill_collisions(p1,collp,atomp);/*np is number of erased collisions between p1 and partner atoms (not including cell crossing event of p1). atomp[0] to atomp[np-1] stores atom index of partner atoms, and collp[0] to collp[np-1] stores collision types of the erased collisions*/
 i1=a[p1].i.x.i-1; /*"left" index along x-axis of neighbor cell where p1 is located*/
 j1=(a[p1].i.y.i-1)*search.x;/*"left" cell index along y-axis*/
 i2=i1+2;/*"left" cell index along z-axis*/
 j2=j1+(search.x<<1);
 
 if(search.z==search.y){/*if two dimensions*/
   k1=a[p1].i.z.i*search.y;
   k2=k1;
 }
 else{
   k1=(a[p1].i.z.i-1)*search.y;
   k2=k1+(search.y<<1);
 }
 
 for(k=k1;k<=k2;k+=search.y){/*look whithin and neighboring cells*/
   addressz=k;
   if(addressz<0)
     addressz+=search.z;
   if(addressz==search.z)
     addressz=0;
   
   for(j=j1;j<=j2;j+=search.x){ 
     addressy=j; 
     if(addressy<0)
       addressy+=search.y;
     if(addressy==search.y)
       addressy=0;
     addressy+=addressz;
     
     for(i=i1;i<=i2;i++){ 
       address=i; 
       if(address<0)
	 address+=search.x;
       if(address==search.x)
	 address=0;
       address+=addressy;
       inpt=search.cell_atoms[address];/*Pointer to the begining of the list */
       while(inpt) /*of atom indexes contained in cell specified by "address"*/
	 { 
	   int p=inpt->n;
	   if(collp[p]==-1)/*if there was no previously scheduled collision*/
	     {
	       atomp[np++]=p;/*we will schedule a collision*/
	       collp[p]=ecoll[a[p].c][a[p1].c];/*store collision type*/
	     }  
	   inpt=inpt->next;
	 }	   	
     }
   }
 }  
 return np;/*return number of atoms for which a collision with p1 will be    */
           /*scheduled                                                       */
}/*Matches int squeeze1(int p1, int * collp, int * atomp)*/

int squeeze_table(int * p1, int * q1, double * timea)/*25 lines of code*/
{ 
  alist *inpt=search.cell_champion[search.olymp[search.final-1][0]];
  search.nch=0;
  search.ncha=0;
  *timea=inpt->t;
  *p1=inpt->p;
  *q1=inpt->q;
  if((*q1)>=search.n)/*if *q1 is a wall collision                            */
  {
   return (int)(inpt->c);
  }
 search.np=squeeze1(*p1,search.collp, search.atomp);
 search.nq=squeeze1(*q1,search.collq, search.atomq);
  return (int)(inpt->c);
}/*Matches int squeeze_table(int * p1, int * q1, double * timea)*/

  
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
  /*  for (i=0;i<search.n;i++)
    for( inpt1=search.atom_collisionp[i];inpt1->p==i;inpt1=inpt1->pnext)
      if(inpt1->q<search.n)
      n_bonds+=do_something(inpt1->p,inpt1->q,inpt1->c);*/
  for (inpt1=search.begin[1];inpt1<*search.free;inpt1++)
    if(inpt1->q<search.n)
      n_bonds+=do_something(inpt1->p,inpt1->q,inpt1->c);
  return n_bonds;
}  

void print_collisions(int p1,int q1)
{
  int np=0;
  alist * pt=search.atom_collisionp[p1];
      fprintf(fp,"collision %lf %d %d %d\n",get_time(),get_ll(),p1,q1);
  while(pt->p==p1)
    {
      fprintf(fp,"%10.6lf %3d %3d %2d\n",pt->t,pt->p,pt->q,pt->c);
      pt=pt->pnext;
    }
  pt=search.atom_collisionq[p1];
  fprintf(fp,"--------\n");
  while(pt->q==p1)
    {
      fprintf(fp,"%10.6lf %3d %3d %2d\n",pt->t,pt->p,pt->q,pt->c);
      pt=pt->qnext;
    }
fflush(fp);
}

void set_search_volume(double maxrb){
  maxrprobed = maxrb * 1.20 ;
}
