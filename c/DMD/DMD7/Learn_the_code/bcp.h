#define NFREE (1000000)
#define DELTALL (100)
#define noErr (0)

#define MOVIE (3)
#define TEXT (2)
#define FILE_IS_OPEN (1)
#define FILE_ERROR (-1)
#define MAXADD (4194304)
#define DBL2 (DBL_MAX*10e-10)
#define DBL1 (DBL_MAX*0.5e-10)
typedef unsigned int size_al;
typedef int well_type;
typedef union {
  int i;
  double d;
} idouble;
typedef struct {
  idouble x;
  idouble y;
  idouble z;
} icrd;
typedef struct {
  double x;
  double y;
  double z;
} crd;
typedef struct {
  double e;      /*Negative of the energy increase listed in the initial     */
	         /*configuration file                                        */
  double eo;     /*Same as "e"                                               */
  double etot;   /*Difence in potential energy from an infinite distance     */
		 /*(zero potential) to the present distance                  */
  double dd;     /*Square of the interparticle distance                      */
  double dm;     /*1/((mi+mj)*coll[j].dd)  or 0 if mi==DBL1 || mj==DBL1      */
  double edm;    /*2*coll[j].e/(mi*mj*coll[j].dm) if non-ellastic or 0 if    */
                 /*ellastic or other expression if static mass               */
  double edmo;   /*Same as "emd"                                             */
  double mij;    /*mi+mj or DBL1 if mi==DBL1 || mj==DBL1                     */
  well_type next;/*                                                          */
  well_type prev;/*For the outer wells of the bonds -1                       */
  well_type react;/*Shows record of reaction and bond status; negative -- if */
                  /*bond less than -1 -- backward reaction, 0 if no reaction */
                  /*and positive if forward reaction                         */
} CollisionData;
/* hard core collisions are sampled at the end of collision list             */
/* soft core colisions are sampled at the end of collision list              */


typedef struct {
  double eo; /*corrected energy surplus, In case of breaking a link*/
/*  double e1; energy surplus: the change in energy after reaction is
computed as surplus + old potential energy - new potential energy, after
energies of all pairs are taken into account the dd and dm are taken from the
collision data */
  double dd;/*square of the maximum distance that the bond can be stretched*/
  int old1;/*type of atom 1 before reaction happens*/
  int old2;
  int new1;/*type of atom 1 after reaction happens*/
  int new2;
  int bond; /* for sake of simplicity of filling the reaction data. Set to 1 when initialized*/
/* since the reaction record is old1 old2 new1 new2 bond [d, eo] 
if [d and eo are missing, they are taken from non-ellastic collision list
for old1 and old2 */
  well_type in; /* collision type after reaction. Set to -1 when initialized*/
  well_type out;/* collision type after reverse reaction. Set to -1 when initialized*/ 
} ReactionData;


typedef struct {
  crd r;
  crd v;
  icrd i;/*cell index*/
  crd q;
  crd u;
  double t; /*last time when atom position was updated*/
  double w; /*time to hit the closest wall of the enclosing cell*/
  double m; /*mass*/
  double b; /*interaction radius*/
  double s; /*ellastic radius*/
  int add;
  int c; /*atom type*/
} atom;
typedef struct {
  double r[3];
  double v[3];
  idouble i[3];
  double q[3];
  double u[3];
  double t;
  double w;
  double m;
  double b;
  double s;
  int add;
  int c;
} iatom;
typedef struct {
  crd q;
  crd v;
  icrd i;
  crd r;
  crd u;
  double t;
  double w;
  double m;
  double b;
  double s;
  int add;
  int c;
} moved_atom;
typedef struct {
  double q[3];
  double v[3];
  idouble i[3];
  double r[3];
  double u[3];
  double t;
  double w;
  double m;
  double b;
  double s;
  int add;
  int c;
} moved_iatom;

typedef struct {
  double r;
  double dummy1;
  double dummy2;
  double v;
  double dummy3;
  double dummy4;
  idouble i;
  idouble dummy5;
  idouble dummy6;
} triad;


typedef struct {
  double length;/*size of the system along this dimension*/
  double dl; /*lentgh of each cell along this dimension*/
  int period;/*number of cells along this dimension*/
} dimensions;

extern int get_delta_ll(void);
extern void set_delta_ll(int new_deltall);
extern double get_rate(void);
extern void set_rate(double rate);
extern int open_movie_file(int is_open,char * fname);
void set_frate(double frate);
double get_frate(void);
void set_mfrate(double frate);
double get_mfrate(void);

extern void moveatoms(void);
extern double countenergy(void);
extern int open_echo_file(int is_open,char * fname);
extern int writetext(char * fname);

extern int twall(int i,double * t);
extern int tball(int i1,int j1,int ct,double * t);
extern int collision_type(int i, int k);

extern double dist(crd r, crd s);


extern void set_maxfree(int a);
int set_text_name(int is_open,  char * fname);
extern int get_maxfree(void);
extern int get_free(void);
extern void set_new_bounds(double * L, double maxrb,int ndim);
extern double get_time(void);
extern void set_temp(double temp);
extern void set_time(double time);
extern double get_mes_time(void);
extern double get_temp_limit(void);
extern double get_temperature(void);
extern void set_temp_limit(double t);
extern double get_coeff(void);
extern void set_coeff(double c);
extern void init_parameters(void);
extern void init_update_param(int * is_x);
extern double get_pressure(void);
extern double get_temp(void);
extern void rescale(void);
extern void add_potential(int ct);
extern int readfile (void);
extern atom * get_atom(void);
extern int get_atom_number(void);
extern int get_dimension(void);
extern dimensions * get_bounds(void);
extern double get_movie_dt(void);
extern int is_reaction(well_type k);
extern int is_bond(well_type k);
extern double get_avePot(void);
extern int write_echo(void);
extern int is_internal(well_type k);
extern double etot(well_type k);
extern double get_corr(void);
extern int get_ll(void);
extern void stop_atoms(moved_iatom * a, int n);
extern void corr_vel(void);
extern void update_atoms(void);
extern void reset_colldata(void);
extern int cleanup(void);
extern void set_timeb(double t);
extern void advance_timeb();
extern double get_timeb(void);
extern double get_timec(void);
extern char * get_text_name(void);
extern well_type ** get_ecoll(void);
extern double atom_dist(int i, int j);


extern void set_press_limit(double P0);
extern void set_press_coeff(double dv_dp);
extern void set_L_limit(double * L);
extern void set_max_gap_mes(int n_coll_V);
extern void set_max_time(double newmaxtime);
extern double get_max_time(void);

void event_loop(void);
