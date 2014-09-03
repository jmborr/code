enum keys {SYS_SIZE=1,NUM_ATOMS,TYPE_ATOMS,NONEL_COL,EL_COL,LINK_PAIRS,REACT,LIST_ATOMS,LIST_BONDS,LIST_PERM_BONDS,NUM_BONDS,LIST_PARAM,COL_TABLE,HB_LIST,GID_LIST};


extern int write_key_coord(FILE *path);
extern int writebonds(FILE * path);
extern int startup(void);
extern atom * get_sample_atoms();
extern int get_atom_types();
extern void set_write_param(int n0,int n1,int yes);
extern void set_hb(int p, int q);
extern void break_hb(int p, int q);
extern void close_hb(atom* a);
