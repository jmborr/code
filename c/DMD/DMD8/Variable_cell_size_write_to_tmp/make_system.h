#ifndef _MS_
#define _MS_

#include "cell_size.h"

enum keys {SYS_SIZE=1,DELTA_RES,CELL_SIZE,SAFE_LIMIT,NUM_ATOMS,TYPE_ATOMS,NONEL_COL,EL_COL,LINK_PAIRS,REACT,LIST_ATOMS,LIST_BONDS,LIST_PERM_BONDS,NUM_BONDS,LIST_PARAM,COL_TABLE,HB_LIST};


extern int write_key_coord(FILE *path);
extern int writebonds(FILE * path);
extern int startup(void);
extern atom * get_sample_atoms();
extern int get_atom_types();
extern void set_write_param(int n0,int n1,int yes);
extern void dump_coll(CollisionData c);
#endif
