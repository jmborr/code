#ifndef _CS_
#define _CS_
#include "bcp.h"

#define NDIM (3)
#define INF_SO (100000) /*by definition, shell order between two atoms permanently linked*/

void set_cell_size(double);

typedef struct
{
  double input_cell_size;
  double sf; /*safe limit, we update sf of each atom after it diffuses sf/2*/
  double dl;
  int nat;
  shell_order amso; /*Absolute Maximum Shell Order*/
  /*mso_tptq is the maximum shell order between type of atom p and
  type of atom q, for p and q unbound, that is, for cells beyond
  mso_tptq[cp][cq], atom p of type cp has no chance to interact with q
  of type cq*/
  shell_order **mso_tptq; 
  /*mso_tp is the maximum shell order for type of atom p, that is, for
    cells beyond mso_tp[cp], atom p of type cp has no chance of
    interacting with any atom to which is not linked*/
  shell_order *mso_tp;
  /*boundary[p][0] indicates the number of surface-shell-orders to
    check when p crosses a cell wall. boundary[p][i]
    (i<boundary[p][0]) is the particular surface shell-order to
    check. For example, with cell-size=1.5A, CA atoms have only one
    surface shell order (so=3) and CB atoms have two (so=5 for CB-CB
    interactions and so=3 for the rest of interactions. Thus, if our
    CB atom was in cell (0,0,0) and then undergoes a wall collision
    whereby it is now in cell (1,0,0), we have to check cells
    (1+3,-3..3,-3..3) for new CA neighbors, and check cells
    (1+5,-5..5,-5..5) for new CB neighbors. Similarly, we check cells
    (1-3,-3..3,-3..3) for old CA neighbors and (1-5,-5..5,-5..5) for
    old CB neighbors. Thus we look to 7*7+11*11=170 cells for new
    neighbors and the same number for old neighbors (note that because
    of periodic boundary conditions, cell (-2,-3,-3) is actually cell
    (N-2,N-3,N-3), and so on).*/
  shell_order **boundaries;
} shell_order_scheme;

void dump_input_cell_size();
void set_input_cell_size(double);
double return_input_cell();
void set_input_cell_size(double cs);
void set_input_safe_limit(double sf);
void set_dl(double dl);
void set_amso(double interaction_range);
int allocate_and_init_shell_order_scheme(double ***coldata,atom *sam,int nat);
shell_order return_mso(atom_number p, atom_number q);
shell_order return_so(atom_number p, atom_number q);
void dump_boundaries(void);
#endif
