#ifndef _SEARCH_
#define _SEARCH_
#include "bcp.h"
#include "cell_size.h"

/*each atom has 5 types of p-collisions (listed below) and q-collisions.

WALL collision schedules the time when atom will cross one of the cell
boundaries

SAFE collision schedules the time when atom will diffuse a distance
higher than half the safe-limit, with respect to the position when we
last updated its safe-limit origin-position (SFOP, stored in "crd sfr"
of structure "atom"). At the beginning of the simulation, the SFOP
coincides with the position of the atom.

HOT and WARM collisions schedule a time when the relative distance
between the interacting pair or atoms, (p,q), will coincide with a
distance showing a discontinuity in the potential energy of (p,q).

A collision is HOT if the collision type is of HOT type. A collision
may be WARM or COLD if the collision type is of WARM type. A collision
type is of HOT type if it corresponds to the inner distance of a
potential well, and the width of the potential well is smaller than
twice the safe-limit distance. A collision type is of WARM type if it
corresponds to the inner distance of a potential well, and the width
of the potential well is bigger than twice the safe-limit distance. A
collision type is also of WARM type if it corresponds to the outermost
potential discontinuity.

A collision between atoms (p,q) is WARM if the collision type is of
WARM type, and the distance between the SFOP of p and the SFOP of q is
smaller than the safe-limit distance plus the distance associated to
the collision type (the distance at which the discontinuity of the
potential happens). If the collision type corresponds to the inner
barrier of a potential well, and the distance between the SFOP of p
and the SFOP of q is higher than the safe-limit distance plus the
distance associated to the collision type, then the collision can
still be WARM if the distance between the SFOP of p and the SFOP of q
is higher than the distance of the outer wall minus the safe-limit
distance.

A collision between atoms (p,q) is COLD if the collision type is of
WARM type and the distance between the SFOP of p and the SFOP of q
does not satisfy the previous requirements to be defined as WARM.  We
are guaranteed that p and q will not undergo a collision unless either
p or q diffuses a distance bigger than half the safe-limit from its
SFOP. Thus, we never have to calculate collision times for COLD
collisions, which saves computing time

Once atom p diffuses over half the safe-limit, we have to check all
the COLD collisions to see if they have become WARM, that is, if p has
become "dangerously" close to other atoms. Similarly, we have to check
the WARM collisions to see if any of these collisions have become
COLD. If a COLD collision is assigned as WARM, then we switch it to
the WARM list of collisions and we calculate its collision time. If a
WARM collision becomes COLD, then we swithch it to the COLD list (and 
there's no need to calculate any collision time).

Every time atom p undergoes a collision with other atoms, we have to update
the times for the HOT, WARM, SAFE, and WALL lists.

Every time atom p crosses a cell wall, we have to look for new
neighbor atoms, and discard neighbors that are now far away. Thus, we
may end up adding and/or removing atoms from the HOT, WARM and COLD
lists of p. For example, say atom p crosses to a cell which is to the
right of the old cell, and along the x-axis. Say also that atom p has
only one type of interaction potential with any other type of atom,
with interaction range D. Say the cell_size is A. Then atoms that are
in the cell 1+(int)(D/A) to the right of the new cell will be new
neighbors of atom p because p may collide with them. Similarly, atoms
that are in the cell 1+(int)(D/A) to the left of the old cell will
not interact with p unless either p crosses back to the old cell or
any of these atoms cross one cell to the right. These atoms are now
old neighbors of p and must be discarded. Thus, we only have to sample
the "surface" of the interaction volume of p every time p crosses a
cell. The situation is more complicated if atom p interacts with
different potentials depending on the type of the interacting
partner. For example, if p is a CB atom, then it has a 7.5Ang
interaction range with other CB's, and it has a 3.4Ang interaction
range with CA's. If this is the case and p crosses a cell like in the
previous example, then we have to check the cell 1+(int)(7.5/A) for
new CB's, and the cell 1+(int)(3.4/A) for new CA's. We don't care if
there are CA atoms in the cell 1+(int)(7.5/A), because they are very
far away.

For example, if atom p (of type CB) is at cell (0,0,0) and atom q (of
type CB) is at cell ( 1+(int)(7.5/A), 1+(int)(7.5/A), 1+(int)(7.5/A)
), then we accept q as neighbor of p. However, note that p and q are
separated by at least a distance of sqrt(3)*7.5, which usually means
that the distance between respective SFOP is also big. If distance
between SFOP's is bigger than 7.5Ang plus the safe-limit distance,
then atom q will be assigned to the COLD list of p and we don't have
to compute collision times for the (p,q) pair. (As a side-note, the
SFOP of atom p does not have to be in the same cell as the cell where
atom p is. The same goes for atom q and its SFOP).

We excluded in this discussion atoms q that are bonded to atom p. The
length of these bonds may well exceed the interaction range of unbound
atoms. For example, the bond between CA_i and CA_i has inner distance
of 3.7Ang, which is higher than the hard-core distance for unbound
CA's. Function isFriend() tells us whether two atoms are bonded. Thus,
if we find atom q in new cell and it turns out to be bonded to p, then
we know it was already accounted for and it is not a new
neighbor. Similarly, we are not to discard bonded atoms q that we find
in the old cells.

At the very beginnig of the simulation, when we're looking for the
neighbors of atom p, then we look for all CA's that are close, ie.,
CA's in cell whose cell indexes are not bigger or smaller than
1+(int)(3.4/A). Thus we search ( 1+(int)(3.4/A) )^3 cells around p. We
will search ( 1+(int)(7.5/A) )^3 cells around atom p when looking for
neighbor CB's.

Lists HOT, WARM, and COLD are doubled. There are so-called plists and
qlists. plists[HOT][p] is a pointer to the beginning of the HOT list
that lists all collisions in which p is one of the partners, and each
of the other partners has an atom number that is bigger than the atom
number of p. By "bigger", we mean bigger in the sense of a list that
closes onto itself. For example, if there are 10 atoms in the system,
then 10 is bigger than 9,8,7, and 6, but smaller than 1,2,3,4, and
5. This is true because we identify atom 1 with (virtual) atom 11. The
list closes onto itself. To determine who is the smallest of a pair
(p,q) pair, see function less_than().

For every plists of atom p for which we calculate collision times,
plists[HOT][p], plists[WARM][p], plists[SAFE][p], and plists[WALL][p],
we find the collision with the smallest time. Then we point to these
collisions with pointers tlists[HOT][p], tlists[WARM][p],
tlists[SAFE][p], and tlists[WALL][p]. Actually, plists[SAFE][p] and
plists[WALL][p] only contain one collision each, so that
tlists[SAFE][p] and tlists[WALL][p] are kind of redundant. Thus we end
up with four winner collisions, which we scan to find who is the
winner of all. if the winner collision corresponds to winner of list
wlt (wlt is either of HOT,WARM,SAFE, and WALL), then tlists[wlt][p]
points to the winner collision under all plists of p. We store the
value of wlt in winlist[p]. If there are N atoms in the system, then
we have N winner times that we sort into a binary tree, in order to
search for the smallest time in the whole system. Most usually the
winner collision is of HOT type.

By splitting the collisions into HOT and WARM lists, we gain
efficiency when updating collision times. Say atoms (p,q) undergo
collisions. Then we have to recompute all times for collisions in
which either p or q are involved. Say a third atom r had a HOT collision
scheduled with atom p, and that atom number of r is smaller than atom number
of p. It may so happen that collision (r,p) was the winner collision
for all plists of atom r. What we have to do is first find the
collision in the list plists[HOT][r] with the smallest time, point
tlists[HOT][r] to this collision, and then compare to tlist[WARM][r],
tlists[WALL][r] and tlists[SAFE][r] to find the winner collision for
all plists of r. Note that we did not have to look in plists[WARM][r]
because collision (r,p) was a HOT collision. In a tipical (p,q)
collision, we find about three other atoms for which we have to update the
winner time because their winner collision involved either p or q.

After we know all the atoms for which the winner time has changed, we
go to the binary tree that sorts all winner times, erase the nodes
related to these atoms, and insert new nodes at the proper places. For
a description of how to build and update the binary tree, read 
"Smith et al. J. Comp. Phys. 134, 16-30 (1997)"


Table of collision events and types of updates for lists
                    |   WALL  |   HOT  |  WARM  |  COLD  |  SAFE |
------------------------------------------------------------------
wall collision      |    T    |    N   |    N   |    N   |       |
                    |         |        |        |        |       |
particle collision  |    T    |    T <===>  T   |        |   T   |
                    |         |        |        |        |       |
hydrogen bond react |    T    |    T <===>  T <===>  T   |       |
                    |         |        |        |        |       |
safe limit crossing |         |        |    N <===>  N   |   T   |
------------------------------------------------------------------
  Legend:   T    update list due to collision time changes
            N    update list to add and/or remove neighbors
          <===>  possibility of switching collisions between lists

NOTE: T's and N's for the hydrogen bond reaction correspond to more than one
      pair of interacting atoms (in total, there are 5 interacting pairs)
*/

#define nlt 5 /*number of list types*/
#define nft 4 /*number of list types with first times*/
enum list_type {HOT=0,WARM,WALL,SAFE,COLD};

typedef struct cell_list
{ 
  int n;
  struct cell_list * next;
} clist;
/* atom list of collisions */

typedef struct atom_list{
  double t;
  int p;
  int q;
  int c;
  shell_order mso; /*maximum shell order between the type of atom p and the
		    the type of atom q*/
  struct atom_list * pnext;
  struct atom_list * qnext;
  struct atom_list * pprev;
  struct atom_list * qprev;
} alist;  

typedef struct{
  clist *atom_storage;
  clist **cell_atoms;
  alist **begin;
  alist **free;
  alist **neighbors;
  enum list_type *neig_lt;
  alist ***plists; /*pointers to all p-collisions lists*/
  alist ***qlists;
  alist ***tlists; /*pointers to first collision of all p-collisions lists, */ 
  /*winlist[p] is the list type of the collision with the winner time
    for all plists*/
  enum list_type *winlist;
  alist * storage; 
  int maxfree;
  int x;
  int y;
  int z;
  int * atomp;
  int * atomq;
  int * collp;
  int * collq;
  int np;
  int nq;
  int n;
  int nhalf;
} tsearch;

/*fast_child(slow_child) points to node with a smaller(bigger)
  collision time than its parent*/
typedef struct node_tree
{ 
  struct node_tree *parent;
  double t;
  struct node_tree *fast_child, *slow_child;
} node;

typedef struct 
{
  int n1;
  int n_updates;
  atom_number p; /*p-partner of the champion collision*/
  atom_number *update_node_for; /*says which nodes have changed their "t"*/
  struct node_tree *storage;
  struct node_tree *root;
}tree_scheme;

extern int allocsearch(int n);
extern void initsearch(void);
extern int init_tables(void);
extern void set_maxfree(int a);
extern int get_maxfree(void);
extern int get_free(void);
extern void find_atoms(int p1, size_al address1);
extern coll_type squeeze_table(atom_number *p1,atom_number *q1,double *timea);
extern int check_atoms(int add);

extern int * get_collp(void);
extern int * get_collq(void);
extern int * get_atomp(void);
extern int * get_atomq(void);
extern int get_np(void);
extern int get_nq(void);
extern int pairs(int (*)(int,int,int));

extern void add_atom2cell(int i);
int less_than(atom_number i,atom_number j);
alist* insert_collision(enum list_type lt, double t, atom_number p,
			atom_number q, coll_type c, shell_order mso);
alist* insert_collision2(enum list_type lt, double t, atom_number i,
			atom_number j, coll_type c, shell_order mso);
int remove_neighbor(atom_number p, enum list_type lt, alist *pt);
void remove_neighbor_noupdlist(atom_number p, enum list_type lt, alist *pt);
enum list_type find_lt(atom_number p, atom_number q, coll_type ct);
enum list_type find_lt2(atom_number p, atom_number q, coll_type ct);
enum list_type return_lt(double dd,coll_type ct);
enum list_type return_lt2(double dd,coll_type ct);
void update_tlist(enum list_type lt, atom_number p);
int update_tlists(alist *pt,enum list_type lt, atom_number p);
void update_winlist(atom_number p);
int init_neighbor_list(alist **neighbors,enum list_type *neig_lt,
		       atom_number p, int *atomp);
void update_cell_atoms(atom_number p, cell_index old, cell_index new);
void update_neighbor_lists(atom_number p, atom_number w);
void mv_coll(alist *pt,enum list_type lt1,enum list_type lt2);
void update_sf(atom_number p);
void update_table(atom_number p, atom_number q, coll_type ct);
alist *find_coll(atom_number i, atom_number j, enum list_type *lt0);
int allocalendar( int n );
void initcalendar(void);
void find_new_champion(void);
void renew_calendar(void);
void delete_node(atom_number p);
void insert_node(atom_number p);
double d2(atom_number p,atom_number q);
void dump_coordinates(void);
void dump_calendar(void);
void dump_calendarII(void);
int check_cell_addresses(atom_number p);
void chell_all_atom_cell_addresses(void);
void dump_info_of(atom_number p);
void print_list_type(enum list_type lt);
char *print_char_list_type(enum list_type lt);
void check_inf_times();
void check_zero_times(void);
void dump_cell_occupants(cell_index address);
char *atom_name (atom_type atc);
void print_CollisionData(coll_type ct);
void print_coll(alist *pt);
void print_coll_info_plist(enum list_type lt,atom_number p,atom_number q);
void print_coll_info_qlist(enum list_type lt,atom_number p,atom_number q);
void print_coll_info_list(enum list_type lt,atom_number p,atom_number q);
void print_coll_info(atom_number p, atom_number q);
void print_plist_coll(enum list_type lt,atom_number p);
void print_qlist_coll(enum list_type lt,atom_number p);
void print_colls(enum list_type lt, atom_number p);
void print_all_colls(atom_number p);
void check_plist_pointers(enum list_type lt, atom_number p);
void check_qlist_pointers(enum list_type lt, atom_number p);
void check_list_pointers(enum list_type lt,atom_number p);
void check_pointers(atom_number p);
void check_all_pointers(void);

int bond(atom_number i, atom_number j);
#endif
