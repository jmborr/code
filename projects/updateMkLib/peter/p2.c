
/*
 * PULCHRA version 2
 *
 * PowerfUL CHain Restoration Algorithm
 * by Piotr Rotkiewicz, 1999-2005
 *
 * piotr@pirx.com
 *
 * cc -O3 p2.c p2_data.c -lm -o p2
 */

#define COMPILE_BB
#define COMPILE_ROT

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>

#include "p2_common.h"

#define PULCHRA_VERSION 2.16

#define uchar char
#define uint unsigned int
#define real double

#define MAX_BUF_SIZE 1000


#define FILE_SUCCESS     0
#define FILE_NOT_FOUND  -1
#define FILE_WARNING    -2

#define FATAL_MAE -1

#define FLAG_BACKBONE  1
#define FLAG_CALPHA    2
#define FLAG_SIDECHAIN 4
#define FLAG_SCM       8

#define FLAG_PROTEIN  1
#define FLAG_DNA      2
#define FLAG_RNA      4
#define FLAG_CHYDRO   8

#define RADDEG 180./M_PI
#define DEGRAD M_PI/180.

int _VERBOSE = 0;
int _CA_OPTIMIZE = 1;
int _CA_RANDOM = 0;
int _CA_ITER = 100;
int _CA_TRAJECTORY = 0;
int _CISPRO = 0;
int _CENTER_CHAIN = 0;
int _REBUILD_BB = 1;
int _REBUILD_SC = 1;
int _REBUILD_H = 1;
int _PDB_SG = 0;
int _TIME_SEED = 0;
int _XVOLUME = 1;
int _XVOL_ITER = 3;
real _CA_START_DIST = 3.0;
real _CA_XVOL_DIST = 4.0;

#define CALC_C_ALPHA
#define CALC_C_ALPHA_ANGLES
#define CALC_C_ALPHA_START
#define CALC_C_ALPHA_XVOL

real CA_K=1.0;
real CA_ANGLE_K=2.0;
real CA_START_K=0.01;
real CA_XVOL_K=0.5;

#define CA_DIST 3.8
#define CA_DIST_TOL 0.1
#define CA_DIST_CISPRO 2.9
#define CA_DIST_CISPRO_TOL 0.1
#define E_EPS 1e-10

#define GRID_RES 6.0

int chain_length = 0;

char AA_NAMES[21][4] =  
  { "GLY", "ALA", "SER", "CYS", "VAL",
    "THR", "ILE", "PRO", "MET", "ASP",
    "ASN", "LEU", "LYS", "GLU", "GLN",
    "ARG", "HIS", "PHE", "TYR", "TRP",
    "UNK" };

char SHORT_AA_NAMES[22] = { "GASCVTIPMDNLKEQRHFYWX" };

int AA_NUMS[256];

int nheavy[20] = { 0, 1, 2, 2, 3, 3, 4, 3, 4, 4, 4, 4, 5, 5, 5, 7, 6, 7, 8, 10};

char *backbone_atoms[4] = { "N  ", "CA ", "C  ", "O  " };

char *heavy_atoms[200]= {
/* GLY */  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ALA */ "CB ", NULL,   NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* SER */ "CB ", "OG ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL, 
/* CYS */ "CB ", "SG ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* VAL */ "CB ", "CG1", "CG2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* THR */ "CB ", "OG1", "CG2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ILE */ "CB ", "CG1", "CG2", "CD1",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* PRO */ "CB ", "CG ", "CD ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* MET */ "CB ", "CG ", "SD ", "CE ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ASP */ "CB ", "CG ", "OD1", "OD2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ASN */ "CB ", "CG ", "OD1", "ND2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* LEU */ "CB ", "CG ", "CD1", "CD2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* LYS */ "CB ", "CG ", "CD ", "CE ", "NZ ",  NULL,  NULL,  NULL,  NULL,  NULL, 
/* GLU */ "CB ", "CG ", "CD ", "OE1", "OE2",  NULL,  NULL,  NULL,  NULL,  NULL,
/* GLN */ "CB ", "CG ", "CD ", "OE1", "NE2",  NULL,  NULL,  NULL,  NULL,  NULL,
/* ARG */ "CB ", "CG ", "CD ", "NE ", "CZ ", "NH1", "NH2",  NULL,  NULL,  NULL,
/* HIS */ "CB ", "CG ", "ND1", "CD2", "CE1", "NE2",  NULL,  NULL,  NULL,  NULL,
/* PHE */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ",  NULL,  NULL,  NULL,
/* TYR */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ", "OH ",  NULL,  NULL,
/* TRP */ "CB ", "CG ", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"};

/* reads full-atom pdb file */

struct _res_type;

typedef struct _atom_type {
  struct _atom_type *next;
  real x, y, z;
  char *name;
  int num, locnum;
  int flag;
  char cispro;
  int gx, gy, gz;
  struct _res_type *res;
  struct _atom_type *prev;
} atom_type;

typedef struct _res_type {
  struct _res_type *next;
  atom_type *atoms;
  int num, locnum, natoms;
  int type;
  char pdbsg;
  char *name;
  char chain;
  real sgx, sgy, sgz;
  real cmx, cmy, cmz;
  struct _res_type *prev;
} res_type;

typedef struct _mol_type {
  struct _mol_type *next;
  res_type *residua;
  int nres;
  unsigned char *r14;
  char *name;
  uchar *seq;
  char **contacts;
  real **cutoffs;
  struct _mol_type *prev;
} mol_type;

mol_type *chain = NULL;

atom_type *new_atom(void)
{
  atom_type *tmpatom;

    tmpatom = (atom_type*) calloc(sizeof(atom_type),1);
    if (tmpatom) {
      tmpatom->x=tmpatom->y=tmpatom->z=0.;
      tmpatom->name=NULL;
      tmpatom->num=tmpatom->locnum=tmpatom->flag=0;
      tmpatom->next=tmpatom->prev=NULL;
    }

  return tmpatom;
}

res_type* new_res(void)
{
  res_type *tmpres;

    tmpres = (res_type*) calloc(sizeof(res_type),1);
    if (tmpres) {
      tmpres->num=0;
      tmpres->name=NULL;
      tmpres->atoms=NULL;
      tmpres->chain=' ';
      tmpres->next=tmpres->prev=NULL;
    }

  return tmpres;
}

mol_type *new_mol(void)
{
  mol_type *tmpmol;
  
    tmpmol = (mol_type*) calloc(sizeof(mol_type),1);
    if (tmpmol) {
      tmpmol->name=NULL;
      tmpmol->residua=NULL;
      tmpmol->next=tmpmol->prev=NULL;
    }  
    
  return tmpmol;
}

void add_atom(atom_type* atomlist, atom_type* newatom)
{
  atom_type *tmpatom;

    if (!atomlist)
      atomlist=newatom;
    else {
      tmpatom=atomlist->next;
      atomlist->next=newatom;
      newatom->prev=atomlist;
      newatom->next=tmpatom;
      if (tmpatom) tmpatom->prev=newatom;
    }
}

void add_res(res_type* reslist, res_type* newres)
{
  res_type *tmpres;

    if (!reslist) 
      reslist=newres;
    else {
      tmpres=reslist->next;
      reslist->next=newres;
      newres->prev=reslist;
      newres->next=tmpres;
      if (tmpres) tmpres->prev=newres;
    }  
}

void add_mol(mol_type* mollist, mol_type* newmol)
{
  mol_type *tmpmol;

    if (!mollist)
      mollist=newmol;
    else {
      tmpmol=mollist->next;
      mollist->next=newmol;
      newmol->prev=mollist;
      newmol->next=tmpmol;
      if (tmpmol) tmpmol->prev=newmol;
    }
}

void delete_atom(atom_type* atom)
{
  atom_type *tmpatom;

    if (atom->prev) atom->prev->next=atom->next;
    if (atom->next) atom->next->prev=atom->prev;
    if (atom->name) free(atom->name);
    free(atom);
    atom=NULL;
}

void delete_res(res_type* res)
{
  res_type *tmpres;
  atom_type *tmpatom;

    if (res->prev) res->prev->next=res->next;
    if (res->next) res->next->prev=res->prev;
    if (res->name) free(res->name);
    if (res->atoms) {
      while (res->atoms) {
        tmpatom = res->atoms->next;
        delete_atom(res->atoms);
        res->atoms=tmpatom;
      }
    }
    free(res);
    res=NULL;
}

void delete_mol(mol_type* mol)
{
  mol_type *tmpmol;
  res_type *tmpres;
  int i;

    if (mol->prev) mol->prev->next=mol->next;
    if (mol->next) mol->next->prev=mol->prev;
    if (mol->name) free(mol->name);
    if (mol->residua) {
      while (mol->residua) {
        tmpres = mol->residua->next;
        delete_res(mol->residua);
        mol->residua=tmpres;
      }
    }
    if (mol->contacts) {
      for (i=0; i<mol->nres; i++) free(mol->contacts[i]);
      free(mol->contacts);
    }
    if (mol->cutoffs) {
      for (i=0; i<mol->nres; i++) free(mol->cutoffs[i]);
      free(mol->cutoffs);
    }
    free(mol);
    mol=NULL;
}


atom_type* get_last_atom(atom_type* atom)
{
    while (atom->next) atom=atom->next;

  return atom;
}

res_type* get_last_res(res_type* res)
{
    while (res->next) res=res->next;

  return res;
}

mol_type *get_last_mol(mol_type* mol)
{
    while (mol->next) mol=mol->next;

  return mol;
}

char setseq(char* aaname)
{
  int i;

    for (i=0; i<21; i++)
      if ((aaname[0]==AA_NAMES[i][0]) &&
          (aaname[1]==AA_NAMES[i][1]) &&
          (aaname[2]==AA_NAMES[i][2]))
         break;
      if (i==21) {
        if (!strcmp(aaname, "GLX"))
          return 'E';
        if (!strcmp(aaname, "ASX"))
          return 'D';
        if (!strcmp(aaname, "HID"))
          return 'H';
        i--;
    }

  return SHORT_AA_NAMES[i];
}

/***
int orient(res_type *res1, res_type *res2)
{
  real x1, y1, z1;
  real x2, y2, z2;
  real len, vect;

    if (!res1 || !res2) return 0;

    if (!res1->prev || !res1->next) {
      x1=0.; y1=0.; z1=0.;
    } else {
      x1=res1->next->sgx-2*res1->sgx+res1->prev->sgx;
      y1=res1->next->sgy-2*res1->sgy+res1->prev->sgy;
      z1=res1->next->sgz-2*res1->sgx+res1->prev->sgz;
    }
    if (!res2->prev || !res2->next) {
      x2=0.; y2=0.; z2=0.;
    } else {
      x2=res2->next->sgx-2*res2->sgx+res2->prev->sgx;
      y2=res2->next->sgy-2*res2->sgy+res2->prev->sgy;
      z2=res2->next->sgz-2*res2->sgz+res2->prev->sgz;
    }

    vect = x1*x2+y1*y2+z1*z2;
    len = sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2+z2*z2);
    if (len) vect /= len;

    if (vect>0.5 && vect<=1.0) return 1; 
    if (vect<=0.5 && vect>-0.5) return 2;

  return 3; 
}
***/

int orient(res_type *res1, res_type *res2)
{
  real x1, y1, z1;
  real x2, y2, z2;
  real cax, cay, caz;
  real len, vect, angle;
  atom_type *atom;

    if (!res1 || !res2) return 0;

    atom=res1->atoms;
    cax=cay=caz=0.;
    while (atom) {
      if (!strncmp(atom->name,"CA",2)) {
        cax=atom->x; cay=atom->y; caz=atom->z;
      }
      atom=atom->next;
    }
    x1=res1->sgx-cax; y1=res1->sgy-cay; z1=res1->sgz-caz;
    if (x1==0. && y1==0. && z1==0.) x1+=1.0;

    atom=res2->atoms;
    cax=cay=caz=0.;
    while (atom) {
      if (!strncmp(atom->name,"CA",2)) {
        cax=atom->x; cay=atom->y; caz=atom->z;
      }
      atom=atom->next;
    }
    x2=res2->sgx-cax; y2=res2->sgy-cay; z2=res2->sgz-caz;
    if (x2==0. && y2==0. && z2==0.) x2+=1.0;

    vect = x1*x2+y1*y2+z1*z2;
    len = sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2+z2*z2);
    if (len) vect /= len;

    angle=RADDEG*acos(vect);

    if (angle>120.) return 1; /*anti*/
    if (angle>60.) return 2;  /*mid*/

//    if (vect>0.5 && vect<=1.0) return 1; /* anti */
//    if (vect<=0.5 && vect>-0.5) return 2; /* mid */

  return 3; /*par*/
}



int res_contact(res_type *res1, res_type *res2) {
  atom_type *atoms1, *atoms2;
  real dx, dy, dz;

    atoms1 = res1->atoms;
    while (atoms1) {
      atoms2 = res2->atoms;
      while (atoms2) {
        dx=atoms1->x-atoms2->x;
        dy=atoms1->y-atoms2->y;
        dz=atoms1->z-atoms2->z;
        if ((atoms1->flag & FLAG_SIDECHAIN) && (atoms2->flag & FLAG_SIDECHAIN) && (dx*dx+dy*dy+dz*dz<20.25)) {
          return 1;
        }
        atoms2=atoms2->next;
      }
      atoms1=atoms1->next;
    }

  return 0;
}



int read_pdb_file(char* filename, mol_type* molecules, char *realname)
{
  FILE *inp;
  char buffer[1000];
  char atmname[10];
  char resname[10];
  char version;
  int prevresnum, resnum, atmnum, locatmnum, num, locnum=0, i, j;
  atom_type *prevatom1, *prevatom2, *prevatom3, *prevatom4;
  int sgnum, cc, nres, ok, natom;
  real sgx, sgy, sgz;
  res_type *res, *test1, *test2;
  atom_type *atom;
  real x, y, z;
  real dist;
  unsigned char bin;
  int warn=0;
  real cutoff;
    
    if (_VERBOSE) printf("Reading input file %s...\n", filename);

    inp = fopen(filename, "r");
    if (!inp) {
      if (_VERBOSE) printf("ERROR: can't open %s !!!\n", filename);
      return FILE_NOT_FOUND;
    }

    molecules->nres=0;
    molecules->name=(char*)malloc(strlen(realname)+1);
    strcpy(molecules->name, realname);

    atmname[3]=0;
    resname[3]=0;
    prevresnum=-666;
    locatmnum=0;
    sgnum=0;
    sgx=sgy=sgz=0.;
    res=NULL;
    while (!feof(inp)) {
      if (fgets(buffer, 1000, inp)!=buffer) break;
      if (!strncmp(buffer, "END", 3) || !strncmp(buffer, "TER", 3)) break; // end of file; only singel molecule read
      if (!strncmp(buffer, "ATOM", 4)) { 
        if (buffer[16]!=' ' && buffer[16]!='A') continue;
        sscanf(&buffer[22], "%d", &resnum);
        strncpy(resname, &buffer[17], 3);
        strncpy(atmname, &buffer[13], 3);
        if (resnum==prevresnum && !strncmp(atmname, "N ", 2)) {        	
          if (_VERBOSE) printf("WARNING: fault in numeration at residuum %s[%d]\n", resname, resnum);
          warn=1;
        }  
        if (resnum!=prevresnum || !strncmp(atmname, "N ", 2)) {
          prevresnum=resnum;
          if (res)
            if (sgnum) {
              res->sgx=sgx/(real)sgnum;
              res->sgy=sgy/(real)sgnum;
              res->sgz=sgz/(real)sgnum;
            } else {
              res->sgx=res->sgy=res->sgz=0.;
            }
          locatmnum=0;
          version=' ';
          res = new_res();
          sgnum=0;
          sgx=sgy=sgz=0.;
          molecules->nres++;
          res->name = malloc(strlen(resname)+1);
          res->type = AA_NUMS[setseq(resname)];
          
          res->locnum=locnum++;
          res->num = resnum;
          res->natoms=0;
          res->chain = buffer[21];
          strcpy(res->name, resname);
          if (molecules->residua) {
            add_res(get_last_res(molecules->residua), res);
          } else {
            molecules->residua = res;
          }
        }
        atom = new_atom();
        atom->res = res;
        res->natoms++;
        locatmnum++;
        sscanf(&buffer[7], "%d", &atmnum);
        sscanf(&buffer[31], "%lf%lf%lf", &x, &y, &z);
        version = buffer[16];
        atom->name = malloc(strlen(atmname)+1);
        strcpy(atom->name, atmname);
        atom->x=x; atom->y=y; atom->z=z;
        atom->num = atmnum;     
        atom->locnum = locatmnum;
        if ((atmname[0]=='S' && atmname[1]=='C')||(atmname[0]=='C' && atmname[1]=='M')) {
          res->cmx = x;
          res->cmy = y;
          res->cmz = z;
          res->pdbsg=1;
        } else
        if (!( ((atmname[0]=='C' || atmname[0]=='N' || atmname[0]=='O') && atmname[1]==' ') ||
               (atmname[0]=='H') ||
               (atmname[0]=='C' && atmname[1]=='A') ||
               (atmname[0]=='O' && atmname[1]=='X' && atmname[2]=='T') ) ) {
          sgx+=x;
          sgy+=y;
          sgz+=z;
          sgnum++;
          atom->flag |= FLAG_SIDECHAIN;
        } else
          atom->flag |= FLAG_BACKBONE;
        if (atmname[0]=='C' && atmname[1]=='A') {
          atom->flag |= FLAG_BACKBONE;
          if (!res->pdbsg) {
            res->cmx = x;
            res->cmy = y;
            res->cmz = z;
          }  
        }           
        if (atmname[0]=='C' && atmname[1]=='M') {
          atom->flag |= FLAG_SCM;
        }           
        if (atmname[0]=='S' && atmname[1]=='C') {
          atom->flag |= FLAG_SCM;
        }           
        if (res->atoms) {
          add_atom(get_last_atom(res->atoms), atom);
        } else {
          res->atoms = atom;
        }
      }
    }

    if (res)
      if (sgnum) {
        res->sgx=sgx/(real)sgnum;
        res->sgy=sgy/(real)sgnum;
        res->sgz=sgz/(real)sgnum;
      } else {
        res->sgx=res->sgy=res->sgz=0.;
      }

    fclose(inp);

    molecules->seq = (uchar*) malloc(sizeof(uchar)*molecules->nres);
    res=molecules->residua; i=0;
    while (res) {    
      molecules->seq[i++]=(uchar)AA_NUMS[(int)setseq(res->name)];
      res=res->next;
    }
   
  if (!warn) return FILE_SUCCESS; else return FILE_WARNING;
}

#define bool int
#define true 1
#define false 0

real calc_ca_energy(atom_type **c_alpha, real **new_c_alpha, real **init_c_alpha, real **gradient, real alpha, real *ene, bool calc_gradient)
{
  int i, j;
  real dx, dy, dz;
  real dist, ddist, ddist2;
  real new_e_pot;
  real theta0, tdif, th, aa, bb, ab;
  real ff0, ff2, dth, m0, m2, grad, f0[3], f2[3];
  real adiff[3], bdiff[3];
  real deriv, theta, dtheta, len1, len2, cos_theta, sin_theta;
  real dx1, dy1, dz1;
  real dx2, dy2, dz2;
  real dx3, dy3, dz3;
  real vx1, vy1, vz1;
  real vx2, vy2, vz2;
  real vx3, vy3, vz3;

  real r12x, r12y, r12z;
  real r32x, r32y, r32z;
  real d12, d32, d12inv, d32inv, c1, c2, diff;
  real f1x, f1y, f1z;
  real f2x, f2y, f2z;
  real f3x, f3y, f3z;
  
        for (i=0; i<chain_length; i++) {
          new_c_alpha[i][0]=c_alpha[i]->x+alpha*gradient[i][0];
          new_c_alpha[i][1]=c_alpha[i]->y+alpha*gradient[i][1];
          new_c_alpha[i][2]=c_alpha[i]->z+alpha*gradient[i][2];
        }   

        new_e_pot = 0.0;
        
        ene[0]=ene[1]=ene[2]=ene[3]=0.0;        
        
        for (i=0; i<chain_length; i++) {

#ifdef CALC_C_ALPHA
          if (i>0) {
            dx=new_c_alpha[i][0]-new_c_alpha[i-1][0];
            dy=new_c_alpha[i][1]-new_c_alpha[i-1][1];
            dz=new_c_alpha[i][2]-new_c_alpha[i-1][2];
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (c_alpha[i]->cispro) {
              ddist=CA_DIST_CISPRO-dist;            
//              if (fabs(ddist)<CA_DIST_CISPRO_TOL) ddist=0.0;
            } else {
              ddist=CA_DIST-dist;            
//              if (fabs(ddist)<CA_DIST_TOL) ddist=0.0;
            }  
            ddist2=ddist*ddist;
            new_e_pot+=CA_K*ddist2;
            ene[0] += CA_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_K)/dist;
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
              gradient[i-1][0]+=grad*dx;
              gradient[i-1][1]+=grad*dy;
              gradient[i-1][2]+=grad*dz;
            }  
          }
#endif

#ifdef CALC_C_ALPHA_START

          dx=new_c_alpha[i][0]-init_c_alpha[i][0];
          dy=new_c_alpha[i][1]-init_c_alpha[i][1];
          dz=new_c_alpha[i][2]-init_c_alpha[i][2];
          dist=sqrt(dx*dx+dy*dy+dz*dz);
          ddist = -dist;
          if (dist>_CA_START_DIST) {
            ddist2=dist*dist;
            new_e_pot+=CA_START_K*ddist2;
            ene[1] += CA_START_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_START_K)/dist;
//printf("re grad: %f\n", grad);
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
            }  
          }
#endif


#ifdef CALC_C_ALPHA_XVOL

          for (j=0;j<chain_length;j++) {
            if (abs(i-j)>3) {
              dx=new_c_alpha[i][0]-new_c_alpha[j][0];
              dy=new_c_alpha[i][1]-new_c_alpha[j][1];
              dz=new_c_alpha[i][2]-new_c_alpha[j][2];
              dist=sqrt(dx*dx+dy*dy+dz*dz);
              ddist = dist;
              if (dist<_CA_XVOL_DIST) {
                ddist2=dist*dist;
                new_e_pot+=CA_XVOL_K*ddist2;
                ene[3] += CA_XVOL_K*ddist2;
                if (calc_gradient) {
                  grad = ddist * (-2.0*CA_XVOL_K)/dist;
//printf("re grad: %f\n", grad);
                  gradient[i][0]-=grad*dx;
                  gradient[i][1]-=grad*dy;
                  gradient[i][2]-=grad*dz;
                  gradient[j][0]+=grad*dx;
                  gradient[j][1]+=grad*dy;
                  gradient[j][2]+=grad*dz;
                }  
              }
            }  
          }
#endif

#ifdef CALC_C_ALPHA_ANGLES

        if (i>0 && i<chain_length-1) {  
          r12x=new_c_alpha[i-1][0]-new_c_alpha[i][0];
          r12y=new_c_alpha[i-1][1]-new_c_alpha[i][1];
          r12z=new_c_alpha[i-1][2]-new_c_alpha[i][2];
          r32x=new_c_alpha[i+1][0]-new_c_alpha[i][0];
          r32y=new_c_alpha[i+1][1]-new_c_alpha[i][1];
          r32z=new_c_alpha[i+1][2]-new_c_alpha[i][2];
          d12 = sqrt(r12x*r12x+r12y*r12y+r12z*r12z);
          d32 = sqrt(r32x*r32x+r32y*r32y+r32z*r32z);
          cos_theta = (r12x*r32x+r12y*r32y+r12z*r32z)/(d12*d32);
          if (cos_theta>1.0) 
            cos_theta = 1.0;
          else
          if (cos_theta<-1.0) 
            cos_theta = -1.0;
          sin_theta = sqrt(1.0-cos_theta*cos_theta);
          theta = acos(cos_theta);

          if (RADDEG*theta<80.) 
            diff = theta-80.*DEGRAD;
          else  
          if (RADDEG*theta>150.) 
            diff = theta-150.*DEGRAD;
          else
            diff = 0.0;
                      
          new_e_pot += CA_ANGLE_K*diff*diff;      
          ene[2] += CA_ANGLE_K*diff*diff;      
          if (calc_gradient) {
            d12inv = 1.0/d12;
            d32inv = 1.0/d32;
            diff *= (-2.0*CA_ANGLE_K)/sin_theta;
            c1 = diff*d12inv;
            c2 = diff*d32inv;
            f1x = c1*(r12x*(d12inv*cos_theta)-r32x*d32inv);
            f1y = c1*(r12y*(d12inv*cos_theta)-r32y*d32inv);
            f1z = c1*(r12z*(d12inv*cos_theta)-r32z*d32inv);
            f3x = c2*(r32x*(d32inv*cos_theta)-r12x*d12inv);
            f3y = c2*(r32y*(d32inv*cos_theta)-r12y*d12inv);
            f3z = c2*(r32z*(d32inv*cos_theta)-r12z*d12inv);
            f2x = -f1x-f3x;
            f2y = -f1y-f3y;
            f2z = -f1z-f3z;
            gradient[i-1][0]+=f1x;
            gradient[i-1][1]+=f1y;
            gradient[i-1][2]+=f1z;
            gradient[i][0]+=f2x;
            gradient[i][1]+=f2y;
            gradient[i][2]+=f2z;
            gradient[i+1][0]+=f3x;
            gradient[i+1][1]+=f3y;
            gradient[i+1][2]+=f3z;
          }
        }
        
#endif
          
        }
        
  return new_e_pot;    
}

#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

real rnd(void)
{
  return 0.001*(real)(rand()%1000);
}


/*
 *  Steepest gradient optimization using v=k*(r0-r)^2
 *  k = CA_K, r0 = CA_DIST
 */
void ca_optimize(char *tname, char *iname)
{
  char buf[1000];
  int i, j, hx, my_iter;
  real dx, dy, dz, dd, dist, dist2, dist3, ddist, ddist2;
  real e_pot, new_e_pot, grad, alpha, e_pot1, e_pot2, e_pot3;
  real adiff[3], bdiff[3];
  real ff0, ff2, aa, ab, bb, th, tdif, dth, m0, m2;
  real theta0, deg_th, maxgrad, sum;
  real f0[3], f2[3];
  real x, y, z;
  int numsteps, numsteps2, msteps;
  int *sec;
  real **new_c_alpha, **gradient, **init_c_alpha, last_alpha, tmp, last_good_alpha, d_alpha, last_e_pot;
  atom_type *atom, **c_alpha;
  res_type *res;
  FILE *inp, *out;
  int mnum, init, ok;
  real alpha1, alpha2, alpha3, a0;
  real ene1, ene2, ene3, e0;  
  real energies[4];      
  real w1, w2, w3, eps;
  real gnorm, last_gnorm;
  int mode, fcnt;
    

    if (_CA_TRAJECTORY) { 
      out = fopen(tname,"w");
      if (out) fclose(out);
    }  
 
    new_c_alpha = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    init_c_alpha = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    for (i=0;i<=chain_length;i++) {
      new_c_alpha[i] = (real*)calloc(sizeof(real)*3,1);
      init_c_alpha[i] = (real*)calloc(sizeof(real)*3,1);
    }  
    gradient = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    for (i=0;i<=chain_length;i++) {
      gradient[i] = (real*)calloc(sizeof(real)*3,1);
    }
              
    c_alpha = (atom_type**)calloc(sizeof(atom_type*)*(chain_length+1),1);
    i = 0;
    res = chain->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A') {
          if (i<chain_length) {
            c_alpha[i] = atom;
            i++;
            break;
          } else {
            if (_VERBOSE) printf("WARNING: number of C-alpha atoms exceeds the chain length!\n");  
            break;
          }
        }
        atom = atom->next;
      }
      res = res->next;
    }

    for (i=0; i<chain_length; i++) {
      init_c_alpha[i][0] = c_alpha[i]->x;
      init_c_alpha[i][1] = c_alpha[i]->y;
      init_c_alpha[i][2] = c_alpha[i]->z;
    }

    if (_CISPRO) {
      for (i=1; i<chain_length; i++) {
        dx = c_alpha[i]->x-c_alpha[i-1]->x;
        dy = c_alpha[i]->y-c_alpha[i-1]->y;
        dz = c_alpha[i]->z-c_alpha[i-1]->z;
        dd = sqrt(dx*dx+dy*dy+dz*dz);
        if ((setseq(c_alpha[i]->res->name)=='P') && (dd>CA_DIST_CISPRO-5*CA_DIST_CISPRO_TOL) && (dd<CA_DIST_CISPRO+5*CA_DIST_CISPRO_TOL)) {
          if (_VERBOSE) printf("Probable cis-proline found at postion %d\n", c_alpha[i]->res->num);
          c_alpha[i]->cispro = 1;
        }  
      }
    }  

    if (_CA_RANDOM) {
      if (_VERBOSE) printf("Generating random C-alpha coordinates...\n");  
      c_alpha[0]->x = 0.0;
      c_alpha[0]->y = 0.0;
      c_alpha[0]->z = 0.0;
      for (i=1;i<chain_length;i++) {
        dx = 0.01*(100-rand()%200);
        dy = 0.01*(100-rand()%200);
        dz = 0.01*(100-rand()%200);
        dd = 3.8/sqrt(dx*dx+dy*dy+dz*dz);
        dx *= dd;
        dy *= dd;
        dz *= dd;        
        c_alpha[i]->x = c_alpha[i-1]->x+dx;
        c_alpha[i]->y = c_alpha[i-1]->y+dy;
        c_alpha[i]->z = c_alpha[i-1]->z+dz;
      }      
    }

    if (iname) {
      inp = fopen(iname,"r");
      if (inp) {
        if (_VERBOSE) printf("Reading initial structure %s...\n", iname);  
        i = 0;
        while (!feof(inp)) {
          if (fgets(buf,1000,inp)==buf && buf[13]=='C' && buf[14]=='A') {
            if (i<chain_length) {
              if (sscanf(&buf[31],"%lf%lf%lf",&x,&y,&z)==3) {
                c_alpha[i]->x = x;
                c_alpha[i]->y = y;
                c_alpha[i]->z = z;
                i++;              
              }  
            } else {
              if (_VERBOSE) printf("WARNING: number of ini-file C-alpha atoms exceeds the chain length!\n");  
              break;
            }
          }
        }
        fclose(inp);
      } else
        if (_VERBOSE) printf("WARNING: can't read initial corrdinates %s\n", iname);
    }
    
    mnum = 1;
    mode = 0;
    init = 0;
    numsteps=numsteps2=0;
    last_alpha = 0.0;
     






//printf("BEFORE:\n");

/*********      

     if (_VERBOSE) {
        for (i=0; i<chain_length; i++) {

#ifdef CALC_C_ALPHA
          if (i>0) {
            dx=c_alpha[i]->x-c_alpha[i-1]->x;
            dy=c_alpha[i]->y-c_alpha[i-1]->y;
            dz=c_alpha[i]->z-c_alpha[i-1]->z;
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (c_alpha[i]->cispro) {
              ddist=CA_DIST_CISPRO-dist;            
              if (fabs(ddist)<CA_DIST_CISPRO_TOL) ddist=0.0;
            } else {
              ddist=CA_DIST-dist;            
              if (fabs(ddist)<CA_DIST_TOL) ddist=0.0;
            }  
            ddist2=ddist*ddist;
       	    if (fabs(ddist)>=CA_DIST_TOL) printf("WARNING: distance %d = %.3lf A\n", i, dist);
          }
#endif
        }
        
        for (i=0; i<chain_length; i++) {
#ifdef CALC_C_ALPHA_ANGLES
          if (i>0 && i<chain_length-1) {
            aa=ab=bb=0.0;
            adiff[0]=c_alpha[i-1]->x-c_alpha[i]->x;
            bdiff[0]=c_alpha[i+1]->x-c_alpha[i]->x;
            aa+=adiff[0]*adiff[0];
            ab+=adiff[0]*bdiff[0];
            bb+=bdiff[0]*bdiff[0];	
            adiff[1]=c_alpha[i-1]->y-c_alpha[i]->y;
            bdiff[1]=c_alpha[i+1]->y-c_alpha[i]->y;
            aa+=adiff[1]*adiff[1];
            ab+=adiff[1]*bdiff[1];
            bb+=bdiff[1]*bdiff[1];	
            adiff[2]=c_alpha[i-1]->z-c_alpha[i]->z;
            bdiff[2]=c_alpha[i+1]->z-c_alpha[i]->z;
            aa+=adiff[2]*adiff[2];
            ab+=adiff[2]*bdiff[2];
            bb+=bdiff[2]*bdiff[2];	

            th=ab/sqrt(aa*bb);
            if (th<-1.0) th=-1.0;
            if (th>1.0) th=1.0;
            th=acos(th);
            deg_th=RADDEG*th;
            if (deg_th>150.) theta0=DEGRAD*150.; else             
            if (deg_th<75.) theta0=DEGRAD*75.; else             
            theta0=th;
       	    if (fabs(deg_th-RADDEG*theta0)>1.0) printf("WARNING: angle %d = %.3lf degrees\n", i, deg_th);
          }
#endif
        }
      }



*********/










    if (_VERBOSE) printf("Optimizing alpha carbons...\n");
 
    eps = 0.5;
    












    fcnt=0;
    
    last_gnorm = 1000.;
    
    do {

      last_e_pot = e_pot;

//printf("my_iter: %d, ca_start_k: %f\n", my_iter, CA_START_K);

      if (_CA_TRAJECTORY) {
        out = fopen(tname,"a");
        if (out) {
          fprintf(out,"MODEL  %d\n",mnum++);
          for (i=0; i<chain_length; i++) {
            fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                    i+1, "CA ", c_alpha[i]->res->name, ' ', c_alpha[i]->res->num,                      
  				          c_alpha[i]->x, c_alpha[i]->y, c_alpha[i]->z);       
          
          }        
          fprintf(out,"ENDMDL\n");
          fclose(out);
        }  
      }

// calculate gradients 

      e_pot=e_pot1=e_pot2=e_pot3=0.;

      for (i=0; i<chain_length; i++)
        gradient[i][0]=gradient[i][1]=gradient[i][2]=0.;

      e_pot = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, 0.0, energies, true);
      
      if (_VERBOSE && !init) {
        printf("Initial energy: bond=%.5lf angle=%.5f restraints=%.5f xvol=%.5f total=%.5f\n", energies[0], energies[2], energies[1], energies[3], e_pot);
      }  
      
      if (!init) init=1;
      
// LINE SEARCH 

      alpha1 = -1.0;
      alpha2 = 0.0;
      alpha3 = 1.0;
        
      ene1 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
      ene2 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha2, energies, false);
      ene3 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);

      while (ene2>MIN(ene1,ene3)) {
        alpha1 *= 2.0;
        alpha3 *= 2.0;
        ene1 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
        ene3 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);
      }
                      
      msteps = 0;      
      do {
        if (alpha3-alpha2>alpha2-alpha1) {
          a0 = 0.5*(alpha2+alpha3);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0-1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0+1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          if (e0<ene2) {
            alpha1 = alpha2;
            alpha2 = a0;
            ene1 = ene2;
            ene2 = e0;
          } else {
            alpha3 = a0;
            ene3 = e0;
          }
        } else {
          a0 = 0.5*(alpha1+alpha2);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0-1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0+1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          if (e0<ene2) {
            alpha3 = alpha2;
            alpha2 = a0;
            ene3 = ene2;
            ene2 = e0;
          } else {
            alpha1 = a0;
            ene1 = e0;
          }          
        }                
        msteps++;
      } while (alpha3-alpha1>1e-6 && msteps<20);

      last_alpha = alpha2;
      e_pot = ene2;
            
//printf("OK! (%d) %.5f : %.50lf\n", msteps, e_pot, last_alpha);
//printf("LAST   : %.5f\n", last_e_pot);

      for (i=0; i<chain_length; i++) {
        c_alpha[i]->x=c_alpha[i]->x+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][0];
        c_alpha[i]->y=c_alpha[i]->y+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][1];
        c_alpha[i]->z=c_alpha[i]->z+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][2];
      }


      e_pot = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, 0.0, energies, false);

      eps *= 0.75;
      if (eps<1e-3) eps=0.0;
      
      numsteps++;
        
        
      gnorm = 0.0;
      for (i=0; i<chain_length; i++) {
        gnorm += gradient[i][0]*gradient[i][0] + gradient[i][1]*gradient[i][1] + gradient[i][2]*gradient[i][2];
      }
      gnorm = sqrt(gnorm/(double)chain_length);
        
//printf("g: %f\n", gnorm);
    
      if (last_gnorm-gnorm<1e-3) fcnt++;
      
      last_gnorm = gnorm;
        
    } while ( (fcnt<3) &&  (gnorm>0.01) && (numsteps<_CA_ITER));


     if (_VERBOSE) {
        for (i=0; i<chain_length; i++) {

#ifdef CALC_C_ALPHA
          if (i>0) {
            dx=c_alpha[i]->x-c_alpha[i-1]->x;
            dy=c_alpha[i]->y-c_alpha[i-1]->y;
            dz=c_alpha[i]->z-c_alpha[i-1]->z;
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (c_alpha[i]->cispro) {
              ddist=CA_DIST_CISPRO-dist;            
              if (fabs(ddist)<CA_DIST_CISPRO_TOL) ddist=0.0;
            } else {
              ddist=CA_DIST-dist;            
              if (fabs(ddist)<CA_DIST_TOL) ddist=0.0;
            }  
            ddist2=ddist*ddist;
       	    if (fabs(ddist)>=CA_DIST_TOL) printf("WARNING: distance %d = %.3lf A\n", i, dist);
          }
#endif
        }
        
        for (i=0; i<chain_length; i++) {
#ifdef CALC_C_ALPHA_ANGLES
          if (i>0 && i<chain_length-1) {
            aa=ab=bb=0.0;
            adiff[0]=c_alpha[i-1]->x-c_alpha[i]->x;
            bdiff[0]=c_alpha[i+1]->x-c_alpha[i]->x;
            aa+=adiff[0]*adiff[0];
            ab+=adiff[0]*bdiff[0];
            bb+=bdiff[0]*bdiff[0];	
            adiff[1]=c_alpha[i-1]->y-c_alpha[i]->y;
            bdiff[1]=c_alpha[i+1]->y-c_alpha[i]->y;
            aa+=adiff[1]*adiff[1];
            ab+=adiff[1]*bdiff[1];
            bb+=bdiff[1]*bdiff[1];	
            adiff[2]=c_alpha[i-1]->z-c_alpha[i]->z;
            bdiff[2]=c_alpha[i+1]->z-c_alpha[i]->z;
            aa+=adiff[2]*adiff[2];
            ab+=adiff[2]*bdiff[2];
            bb+=bdiff[2]*bdiff[2];	

            th=ab/sqrt(aa*bb);
            if (th<-1.0) th=-1.0;
            if (th>1.0) th=1.0;
            th=acos(th);
            deg_th=RADDEG*th;
            if (deg_th>150.) theta0=DEGRAD*150.; else             
            if (deg_th<75.) theta0=DEGRAD*75.; else             
            theta0=th;
       	    if (fabs(deg_th-RADDEG*theta0)>1.0) printf("WARNING: angle %d = %.3lf degrees\n", i, deg_th);
          }
#endif
        }
      }
    
    if (_VERBOSE) printf("Optimization done after %d step(s).\nFinal energy: bond=%.5lf angle=%.5f restraints=%.5f xvol=%.5f total=%.5f\n", numsteps, energies[0], energies[2], energies[1], energies[3], e_pot);

    if (_CA_TRAJECTORY) {
      out = fopen(tname,"a");
      if (out) {
        fprintf(out,"END\n");
      }
    }  

    for (i=0;i<chain_length+1;i++) {
      free(init_c_alpha[i]);
      free(new_c_alpha[i]);
      free(gradient[i]);
    }
    free(new_c_alpha);
    free(gradient);
    free(c_alpha);
    free(init_c_alpha);
}

void center_chain(mol_type *mol)
{
  real cx, cy, cz;
  int natom;
  res_type *res;
  atom_type *atom;
  
    cx = cy = cz = 0.0;
    natom = 0;
    
    res = mol->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {     
        cx += atom->x;
        cy += atom->y;
        cz += atom->z;
        natom++;
  			atom=atom->next;
  		}		          
      res = res->next;
    } 
    
    cx /= (real)natom;
    cy /= (real)natom;
    cz /= (real)natom;
    
    if (_VERBOSE) printf("Molecule center: %8.3f %8.3f %8.3f -> 0.000 0.000 0.000\n", cx, cy, cz);

    res = mol->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {     
        atom->x -= cx;
        atom->y -= cy;
        atom->z -= cz;
        natom++;
  			atom=atom->next;
  		}		          
      res = res->next;
    } 
    
}

void write_pdb(char *name, mol_type *mol)
{
  FILE *out;
  res_type *res;
  atom_type *atom;
  int anum;
    
    out = fopen(name,"w");
    if (!out) {
      if (_VERBOSE) printf("Can't open output file!\n");
      return;
    }
    fprintf(out,"REMARK  REBUILT BY PULCHRA V.%.2f\n", PULCHRA_VERSION);
    anum=1;
    res = mol->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {     
        if (!(atom->name[0]=='D' && atom->name[1]=='U') &&
            !(atom->name[0]=='S' && atom->name[1]=='C') &&
            !(atom->name[0]=='C' && atom->name[1]=='M'))
          fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                        anum++, atom->name, res->name, ' ', res->num,                      
  	  				          atom->x, atom->y, atom->z);       
  			atom=atom->next;
  		}		          
      res = res->next;
    } 
    fprintf(out,"TER\nEND\n");   
    fclose(out);
}


void write_pdb_sg(char *name, mol_type *mol)
{
  FILE *out;
  res_type *res;
  atom_type *atom;
  int anum;
    
    out = fopen(name,"w");
    if (!out) {
      if (_VERBOSE) printf("Can't open output file!\n");
      return;
    }
    fprintf(out,"REMARK  REBUILT BY PULCHRA V.%.2f\n", PULCHRA_VERSION);
    anum=1;
    res = mol->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {     
        if ((atom->name[0]=='C' && atom->name[1]=='A'))
          fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                        anum++, atom->name, res->name, ' ', res->num,                      
  	  				          atom->x, atom->y, atom->z);       
  			atom=atom->next;
  		}		          
      fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                    anum++, "CM ", res->name, ' ', res->num,                      
  				          res->cmx, res->cmy, res->cmz);       
      res = res->next;
    } 
    fprintf(out,"TER\nEND\n");   
    fclose(out);
}



real calc_distance(real x1, real y1, real z1, 
							 		  real x2, real y2, real z2)
{							 		  
  real dx,dy,dz;
  real dist2;

    dx = (x1) - (x2);
    dy = (y1) - (y2);
    dz = (z1) - (z2);
    if (dx || dy || dz ) {
      dist2 = dx*dx+dy*dy+dz*dz;
      return (sqrt(dist2));
    } else
      return 0.0;
}

real calc_r14(real x1, real y1, real z1, 
							 real x2, real y2, real z2,
							 real x3, real y3, real z3,
							 real x4, real y4, real z4)
{
  real r, dx, dy, dz;
  real vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3;
  real hand;
  
    dx = x4-x1;
    dy = y4-y1;
    dz = z4-z1;

    r = sqrt(dx*dx+dy*dy+dz*dz);

    vx1=x2-x1;
    vy1=y2-y1;
    vz1=z2-z1;
    vx2=x3-x2;
    vy2=y3-y2;
    vz2=z3-z2;
    vx3=x4-x3;
    vy3=y4-y3;
    vz3=z4-z3;
	 
    hand = (vy1*vz2-vy2*vz1)*vx3+
           (vz1*vx2-vz2*vx1)*vy3+
           (vx1*vy2-vx2*vy1)*vz3;
  
    if (hand<0) r=-r;

  return r;
}

real superimpose2(real **coords1, real **coords2, int npoints, real **tpoints, int ntpoints)
{
  real mat_s[3][3], mat_a[3][3], mat_b[3][3], mat_g[3][3];
  real mat_u[3][3], tmp_mat[3][3];
  real val, d, alpha, beta, gamma, x, y, z;
  real cx1, cy1, cz1, cx2, cy2, cz2, tmpx, tmpy, tmpz;
  int i, j, k, n;

    cx1=cy1=cz1=cx2=cy2=cz2=0.;

    for (i=0; i<npoints; i++) {
      cx1+=coords1[i][0];
      cy1+=coords1[i][1];
      cz1+=coords1[i][2];
      cx2+=coords2[i][0];
      cy2+=coords2[i][1];
      cz2+=coords2[i][2];
    }

    cx1/=(real)npoints;
    cy1/=(real)npoints;
    cz1/=(real)npoints;

    cx2/=(real)npoints;
    cy2/=(real)npoints;
    cz2/=(real)npoints;

    for (i=0; i<npoints; i++) {
      coords1[i][0]-=cx1;
      coords1[i][1]-=cy1;
      coords1[i][2]-=cz1;
      coords2[i][0]-=cx2;
      coords2[i][1]-=cy2;
      coords2[i][2]-=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]-=cx2;
      tpoints[i][1]-=cy2;
      tpoints[i][2]-=cz2;
    }

    for (i=0; i<3; i++)
      for (j=0; j<3; j++) {
        if (i==j)
          mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=1.0;
        else
          mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=0.0;
        mat_u[i][j]=0.;
      }

    for (n=0; n<npoints; n++) {
      mat_u[0][0]+=coords1[n][0]*coords2[n][0];
      mat_u[0][1]+=coords1[n][0]*coords2[n][1];
      mat_u[0][2]+=coords1[n][0]*coords2[n][2];
      mat_u[1][0]+=coords1[n][1]*coords2[n][0];
      mat_u[1][1]+=coords1[n][1]*coords2[n][1];
      mat_u[1][2]+=coords1[n][1]*coords2[n][2];
      mat_u[2][0]+=coords1[n][2]*coords2[n][0];
      mat_u[2][1]+=coords1[n][2]*coords2[n][1];
      mat_u[2][2]+=coords1[n][2]*coords2[n][2];
    }

    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        tmp_mat[i][j]=0.;

    do {
      d=mat_u[2][1]-mat_u[1][2];
      if (d==0) alpha=0; else alpha=atan(d/(mat_u[1][1]+mat_u[2][2]));
      if (cos(alpha)*(mat_u[1][1]+mat_u[2][2])+sin(alpha)*(mat_u[2][1]-mat_u[1][2])<0.0)       alpha+=M_PI;	
      mat_a[1][1]=mat_a[2][2]=cos(alpha);
      mat_a[2][1]=sin(alpha);
      mat_a[1][2]=-mat_a[2][1];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_a[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_a[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[0][2]-mat_u[2][0];
      if (d==0) beta=0; else beta=atan(d/(mat_u[0][0]+mat_u[2][2]));
      if (cos(beta)*(mat_u[0][0]+mat_u[2][2])+sin(beta)*(mat_u[0][2]-mat_u[2][0])<0.0) beta+=M_PI;
      mat_b[0][0]=mat_b[2][2]=cos(beta);
      mat_b[0][2]=sin(beta);
      mat_b[2][0]=-mat_b[0][2];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_b[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_b[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[1][0]-mat_u[0][1];
      if (d==0) gamma=0; else gamma=atan(d/(mat_u[0][0]+mat_u[1][1]));
      if (cos(gamma)*(mat_u[0][0]+mat_u[1][1])+sin(gamma)*(mat_u[1][0]-mat_u[0][1])<0.0)             
        gamma+=M_PI;	
      mat_g[0][0]=mat_g[1][1]=cos(gamma);
      mat_g[1][0]=sin(gamma);
      mat_g[0][1]=-mat_g[1][0];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_g[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_g[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      val=fabs(alpha)+fabs(beta)+fabs(gamma);
    } while (val>0.001);

    val=0.;
    for (i=0; i<npoints; i++) {
      x=coords2[i][0];
      y=coords2[i][1];
      z=coords2[i][2];
      tmpx=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tmpy=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tmpz=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
      x=coords1[i][0]-tmpx;
      y=coords1[i][1]-tmpy;
      z=coords1[i][2]-tmpz;
      val+=x*x+y*y+z*z;
    }

    for (i=0; i<ntpoints; i++) {
      x=tpoints[i][0];
      y=tpoints[i][1];
      z=tpoints[i][2];
      tpoints[i][0]=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tpoints[i][1]=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tpoints[i][2]=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
    }
      
    for (i=0; i<npoints; i++) {
      coords1[i][0]+=cx1;
      coords1[i][1]+=cy1;
      coords1[i][2]+=cz1;
      coords2[i][0]+=cx2;
      coords2[i][1]+=cy2;
      coords2[i][2]+=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]+=cx1;
      tpoints[i][1]+=cy1;
      tpoints[i][2]+=cz1;
    }

  return sqrt(val/(real)npoints);
}

 
void add_replace(res_type *res, char *aname, real x, real y, real z, int flags)
{
  atom_type *atom, *newatom;
  
//    printf("ADD/REPLACE\n");

    atom = res->atoms;
    while (atom) {
      if (atom->name[0]==aname[0] && atom->name[1]==aname[1] && atom->name[2]==aname[2]) {
        atom->x = x; atom->y = y; atom->z = z;
        atom->flag |= flags;
        break;
      }
      atom = atom->next;
    } 
    
    if (!atom) {
//      printf("NOT FOUND\n");
//      printf("ADDING NEW ATOM\n");
      newatom = (atom_type*)calloc(sizeof(atom_type),1);
      newatom->x = x;
      newatom->y = y;
      newatom->z = z;
      newatom->flag |= flags;
      newatom->res = res;
      newatom->name = (char*)calloc(4,1);
      strcpy(newatom->name,aname);

      atom = res->atoms;        
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A')
          break;
        atom = atom->next;
      }  
      if (aname[0]=='N' && aname[1]==' ') {
        newatom->next = res->atoms;
        res->atoms = newatom;
      } else {
        while (atom->next) atom=atom->next;
        atom->next = newatom;
      }  
    }
    
}


int **RBINS;
real **X_COORDS, **C_ALPHA;

#ifdef COMPILE_BB

void rebuild_backbone(void)
{

  res_type *res, *prevres;
  atom_type *atom;
  real **cacoords, **tmpcoords, **tmpstat;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real r13_1, r13_2, r14;
  real besthit, hit;
  int bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real rmsd, total, maxrms;
  FILE *debug, *out;

    if (_VERBOSE) printf("Rebuilding backbone...\n");

    RBINS = (int**)calloc(sizeof(int*)*(chain_length+1),1);
    for (i=0;i<chain_length+1;i++)
      RBINS[i] = (int*)calloc(sizeof(int)*3,1);

    X_COORDS = (real**)calloc(sizeof(real*)*(chain_length+10),1);
    for (i=0;i<chain_length+10;i++)
      X_COORDS[i] = (real*)calloc(sizeof(real)*3,1);
      
    i = 5;
    res = chain->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A') {
          X_COORDS[i][0] = atom->x;
          X_COORDS[i][1] = atom->y;
          X_COORDS[i][2] = atom->z;
          i++;
        }
        atom = atom->next;
      }
      res = res->next;
    }
    
    cacoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpcoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpstat = (real**)calloc(sizeof(real*)*(8),1);
    for (i=0;i<8;i++) {
      cacoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpcoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpstat[i] = (real*)calloc(sizeof(real)*3,1);;
    }
        
    C_ALPHA = &X_COORDS[5];  

    // rebuild ends...

//printf("adding ends 1...\n");
    
    for (i=0,j=0;i<5;i++,j++) 
      for (k=0;k<3;k++) 
        tmpcoords[j][k] = C_ALPHA[i][k];
    for (i=2,j=0;i<5;i++,j++) 
      for (k=0;k<3;k++) 
        cacoords[j][k] = C_ALPHA[i][k];
    for (i=0,j=0;i<3;i++,j++) 
      for (k=0;k<3;k++) 
        tmpstat[j][k] = C_ALPHA[i][k];

    superimpose2(tmpstat,cacoords,3,tmpcoords,5);

    for (i=-2,j=0;i<0;i++,j++) 
      for (k=0;k<3;k++) 
        C_ALPHA[i][k] = tmpcoords[j][k];

//printf("adding ends 2...\n");

    for (i=chain_length-5,j=0;i<chain_length;i++,j++) 
      for (k=0;k<3;k++) 
        tmpcoords[j][k] = C_ALPHA[i][k];
    for (i=chain_length-5,j=0;i<chain_length-2;i++,j++) 
      for (k=0;k<3;k++) 
        cacoords[j][k] = C_ALPHA[i][k];
    for (i=chain_length-3,j=0;i<chain_length;i++,j++) 
      for (k=0;k<3;k++) 
        tmpstat[j][k] = C_ALPHA[i][k];

    superimpose2(tmpstat,cacoords,3,tmpcoords,5);

    for (i=chain_length-3,j=0;i<chain_length;i++,j++) 
      for (k=0;k<3;k++) 
        C_ALPHA[i+3][k] = tmpcoords[j+3][k];

   
    prevres = NULL;
    res = chain->residua;


    total = maxrms = 0.0;
            
//printf("rebuild...\n");

    for (i=0;i<chain_length+1;i++) {
    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];

    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];

    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];

    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];
      
    	r13_1 = calc_distance(x1, y1, z1, x3, y3, z3);
    	r13_2 = calc_distance(x2, y2, z2, x4, y4, z4);
    	r14 = calc_r14(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

      bin13_1 = (int)((r13_1-4.6)/0.3); 
      bin13_2 = (int)((r13_2-4.6)/0.3); 
      bin14 = (int)((r14+11.)/0.3);
      
      if (bin13_1<0) bin13_1=0;
      if (bin13_2<0) bin13_2=0;
      if (bin14<0) bin14=0;
      if (bin13_1>9) bin13_1=9;
      if (bin13_2>9) bin13_2=9;
      if (bin14>73) bin14=73;

      RBINS[i][0] = bin13_1;
      RBINS[i][1] = bin13_2;
      RBINS[i][2] = bin14;
            
    	cacoords[0][0] = x1;
    	cacoords[0][1] = y1;
    	cacoords[0][2] = z1;

    	cacoords[1][0] = x2;
     	cacoords[1][1] = y2;
     	cacoords[1][2] = z2;

     	cacoords[2][0] = x3;
     	cacoords[2][1] = y3;
     	cacoords[2][2] = z3;

     	cacoords[3][0] = x4;
     	cacoords[3][1] = y4;
     	cacoords[3][2] = z4;

      pro = 0;

      if (prevres && !strncmp(prevres->name,"PRO",3)) {
        j=0;
        besthit=1000.;      
        bestpos=0;
        do {
          hit = fabs(nco_stat_pro[j].bins[0]-bin13_1)+fabs(nco_stat_pro[j].bins[1]-bin13_2)+0.2*fabs(nco_stat_pro[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat_pro[j].bins[0]>=0 && hit>1e-3);
        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  }
     	  }
        for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  }
     	  }
      } else {
        j=0;
        besthit=1000.;      
        bestpos=0;
        do {
          hit = fabs(nco_stat[j].bins[0]-bin13_1)+fabs(nco_stat[j].bins[1]-bin13_2)+0.2*fabs(nco_stat[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat[j].bins[0]>=0 && hit>1e-3);
        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat[bestpos].data[j][k];
       	  }
     	  }
        for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat[bestpos].data[j][k];
       	  }
     	  }
      }

//printf("besthit for %d: %f (%d)\n", i, besthit, bestpos);      

     	rmsd=superimpose2(cacoords, tmpstat, 4, tmpcoords, 8);

     	total += rmsd;
     	if (rmsd>maxrms) maxrms=rmsd;
     	
// add-or-replace
       
      if (prevres) {
        add_replace(prevres, "C  ", tmpcoords[4][0], tmpcoords[4][1], tmpcoords[4][2], FLAG_BACKBONE);
        add_replace(prevres, "O  ", tmpcoords[5][0], tmpcoords[5][1], tmpcoords[5][2], FLAG_BACKBONE);
      }
        
      if (res) {
        add_replace(res, "N  ", tmpcoords[6][0], tmpcoords[6][1], tmpcoords[6][2], FLAG_BACKBONE);
      }

      prevres = res;
      if (res) 
        res = res->next;
    }
    
    if (_VERBOSE) printf("Backbone rebuilding accuracy: average = %.3f, max = %.3f\n", total/(real)chain_length, maxrms);
}

#endif


#ifdef COMPILE_ROT

typedef struct _rot_struct {
  int r13_1, r13_2, r14;
  int nc;
  real ***coords;
  struct _rot_struct *next;
} rot_struct;

rot_struct *rotamers[20];


int read_rotamers(void)
{
  FILE *inp;
  char buf[1000];
  char dum[100];
  int aa, i, j, k, l, n;
  rot_struct *new_rot, *last_rot;
  real x, y, z;
  
    if (_VERBOSE) printf("Reading rotamer library...\n");
    
    inp = fopen("NEWROT","r");
    last_rot=NULL;    
    while (!feof(inp)) {
      if (fgets(buf,1000,inp)==buf) {
        if (buf[0]=='A') {
          sscanf(buf,"%s %d", dum, &aa);
          if (last_rot) last_rot->next = NULL;
          last_rot = NULL;
          if (fgets(buf,1000,inp)!=buf) break;
        }        
//        printf("aa: %d\n", aa);
        if (aa==20) break;
        sscanf(buf,"%d %d %d %s %d", &i, &j, &k, dum, &l);
        new_rot = (rot_struct*)calloc(sizeof(rot_struct),1);
//        printf("%d %d %d nc: %d\n", i, j, k, l);
        new_rot->r13_1 = i;
        new_rot->r13_2 = j;
        new_rot->r14 = k;
        new_rot->nc = l;
        new_rot->next = NULL;
        new_rot->coords = (real***)calloc(sizeof(real**)*l,1);
        for (i=0;i<l;i++) {
          new_rot->coords[i]=(real**)calloc(sizeof(real*)*(nheavy[aa]+1),1);
          for (j=0;j<(nheavy[aa]+1);j++) {
            new_rot->coords[i][j]=(real*)calloc(sizeof(real)*3,1);
          }          
        }
        for (i=0;i<l;i++) {
          fgets(buf,1000,inp);
          for (j=0;j<(nheavy[aa]+1);j++) {
            fgets(buf,1000,inp);
            sscanf(buf,"%lf%lf%lf",&x, &y, &z);
            new_rot->coords[i][j][0]=x;
            new_rot->coords[i][j][1]=y;
            new_rot->coords[i][j][2]=z;
          }
          if (last_rot) {
            last_rot->next = new_rot;
          } else {
            rotamers[aa] = new_rot;
          }
          last_rot = new_rot;
        }          
      }
    }
    fclose(inp);
}


void cross(real *v1, real *v2, real *v3)
{
  v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
  v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

void norm(real *v)
{
  real d;
  
    d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
}


int check_xvol(res_type *res) 
{
  res_type *res2;
  atom_type *atom1, *atom2;
  real dx, dy, dz, dd;
  
    res2 = chain->residua;
    
    while (res2) {
      atom2 = res2->atoms;
      if (res!=res2) {
        while (atom2) {
          atom1 = res->atoms;
          while (atom1) {
            if (atom1->flag & FLAG_SIDECHAIN) {
              dx = atom1->x-atom2->x;
              dy = atom1->y-atom2->y;
              dz = atom1->z-atom2->z;
              dd = dx*dx+dy*dy+dz*dz;
              if (dd<(1.7*1.7)) {
//printf("xvol violation %s[%d]->%s --- %s[%d]->%s\n", res->name, res->num, atom1->name, res2->name, res2->num, atom2->name);
                return 1;
              }
            }
            atom1=atom1->next;
          }
          atom2=atom2->next;
        }
      }  
      res2=res2->next;
    }  

  return 0;  
}


real ***SORTED_ROTAMERS;

void rebuild_sidechains(void)
{
  FILE *out;
  res_type *res, *prevres, *testres;
  atom_type *atom, *atom1, *atom2;
  real **cacoords, **tmpcoords, **tmpstat;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real x5, y5, z5;
  real r14, r13_1, r13_2;
  real dx, dy, dz, dd;
  real hit, besthit;
  int exvol, bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14;
  real rmsd, total;
  real v1[3], v2a[3], v2b[3], v2[3], v3[3];
  int nsc, nca;
  real cax, cay, caz;
  real **lsys, **vv, **sc;
  char scn[12][4];
  rot_struct *rot;
  int ok, last_a, last_b, last_c, last_d, jpos;
  int jx, jy, jz, jxi, jyi, jzi, b13_1, b13_2, b14, jm;
  int crot, bestrot, minexvol, totexvol, rtried, pos, cpos;    
  real cmx, cmy, cmz, ddx, ddy, ddz, ddd, bestdd;
  real sort_rot[100][2];
  
    if (_VERBOSE) printf("Rebuilding side chains...\n");

    lsys = (real**)calloc(sizeof(real*)*3,1);
    vv = (real**)calloc(sizeof(real*)*3,1);
    sc = (real**)calloc(sizeof(real*)*12,1);
    for (i=0;i<12;i++)
      sc[i] = (real*)calloc(sizeof(real)*3,1);
    for (i=0;i<3;i++) {
      lsys[i] = (real*)calloc(sizeof(real)*3,1);
      vv[i] = (real*)calloc(sizeof(real)*3,1);
    }

    SORTED_ROTAMERS = (real***)calloc(sizeof(real**)*(chain_length+1),1);
    for (i=0;i<chain_length+1;i++) {
      SORTED_ROTAMERS[i] = (real**)calloc(sizeof(real*)*10,1);
      for (j=0;j<10;j++) {
        SORTED_ROTAMERS[i][j] = (real*)calloc(sizeof(real)*2,1);
      }
    }
    
    prevres = NULL;
    res = chain->residua;
    totexvol = 0;
    
    for (i=0;i<chain_length;i++) {
      if (!strncmp(res->name,"GLY",3)) {
        if (res->next) res = res->next;
        continue;
      }

//printf("res %d\n", i);
      
    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];
    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];
    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];
    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];

      bin13_1 = RBINS[i][0];
      bin13_2 = RBINS[i][1];
      bin14 = RBINS[i][2];

      v1[0] = x4-x2;
      v1[1] = y4-y2;
      v1[2] = z4-z2;
      
      v2a[0] = x4-x3;
      v2a[1] = y4-y3;
      v2a[2] = z4-z3;

      v2b[0] = x3-x2;
      v2b[1] = y3-y2;
      v2b[2] = z3-z2;
      
      cross(v2a, v2b, v2);
      cross(v1, v2, v3);
      
      norm(v1);
      norm(v2);
      norm(v3);


// gather 10 closest backbone conformations...

      for (j=0;j<10;j++)
        SORTED_ROTAMERS[i][j][0] = 500.;
        
      j = 0;
      besthit = 1000.;
      bestpos = 0;
      do {
        if (rot_stat_idx[j][0]==res->type) {
          hit = fabs(rot_stat_idx[j][1]-bin13_1)+fabs(rot_stat_idx[j][2]-bin13_2)+0.2*fabs(rot_stat_idx[j][3]-bin14);
          if (hit<SORTED_ROTAMERS[i][9][0]) {
            k = 9;
            while (k>=0 && hit<SORTED_ROTAMERS[i][k][0]) {              
              k--;
            }
            k++;
            // k = hit
            for (l=9;l>k;l--) {
              SORTED_ROTAMERS[i][l][0]=SORTED_ROTAMERS[i][l-1][0];
              SORTED_ROTAMERS[i][l][1]=SORTED_ROTAMERS[i][l-1][1];
            }
            SORTED_ROTAMERS[i][k][0]=hit;
            SORTED_ROTAMERS[i][k][1]=j;
          }
        }
        j++;
      } while (rot_stat_idx[j][0]>=0);


//for (j=0;j<10;j++) {
//  printf("best hits: %d %10d %f\n", j, (int)SORTED_ROTAMERS[i][j][1],SORTED_ROTAMERS[i][j][0]);
//}

      besthit = SORTED_ROTAMERS[i][0][0];
      bestpos = SORTED_ROTAMERS[i][0][1];

// new rebuilding procedure...

      pos = rot_stat_idx[bestpos][5];      
      nsc = nheavy[res->type]+1;
      
      if (_PDB_SG) { // more than one rotamer - check SC
        bestdd = 100.; crot = 0;
        for (l=0;l<2;l++) { // check two closest conformations
          cpos = SORTED_ROTAMERS[i][l][1];
          for (m=0;m<rot_stat_idx[cpos][4];m++) {
            for (j=0;j<3;j++) {
              vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];        
              for (k=0;k<3;k++) {
                if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
              }        
            }
            pos = rot_stat_idx[cpos][5]+nsc*m;
            for (j=0;j<nsc;j++) {
              for (k=0;k<3;k++) {
                sc[j][k] = rot_stat_coords[pos+j][k];
              }
            }         
            superimpose2(vv,lsys,3,sc,nsc);      
            for (j=0;j<nsc;j++) {
              sc[j][0] += x3;
              sc[j][1] += y3;
              sc[j][2] += z3;
            }
            cmx = 0.; cmy = 0.; cmz = 0.;
            for (j=0;j<nsc;j++) {
              cmx += sc[j][0];
              cmy += sc[j][1];
              cmz += sc[j][2];
            }
            cmx /= (real) nsc;
            cmy /= (real) nsc;
            cmz /= (real) nsc;
            ddx = res->cmx-cmx;
            ddy = res->cmy-cmy;
            ddz = res->cmz-cmz;
            ddx *= ddx;
            ddy *= ddy;
            ddz *= ddz;
            ddd = ddx+ddy+ddz;
            if (ddd<bestdd) {
              bestdd = ddd;
              crot = pos; // closest rotamer position
            }                        
          }
        }
        pos = crot;
      } // PDB_SG
      
//printf("coords pos: %d\n", pos);
/******
  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                1, "CA ", "GLY", ' ', res->num-2,                      
  	  				  x1, y1, z1);

  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                1, "CA ", "GLY", ' ', res->num-1,                      
  	  				  x2, y2, z2);

  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                1, "CA ", res->name, ' ', res->num,                      
  	  				  x3, y3, z3);

  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                1, "CA ", "GLY", ' ', res->num+1,                      
  	  				  x4, y4, z4);

cpos = SORTED_ROTAMERS[i][0][1];
for (m=0;m<rot_stat_idx[cpos][4];m++) {
printf("REMARK %s[%d] %d\n", res->name, i+1, m+1);
printf("TER\n");
for (j=0;j<3;j++) {
  vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];        
  for (k=0;k<3;k++) {
    if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
  }        
}
pos = rot_stat_idx[cpos][5]+nsc*m;
for (j=0;j<nsc;j++) {
  for (k=0;k<3;k++) {
    sc[j][k] = rot_stat_coords[pos+j][k];
  }
}         
superimpose2(vv,lsys,3,sc,nsc);
for (j=0;j<nsc;j++) {
  sc[j][0] += x3;
  sc[j][1] += y3;
  sc[j][2] += z3;
}
  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                1, "CA ", res->name, ' ', res->num,                      
  	  				  x3, y3, z3);
for (j=1;j<nsc;j++) {
  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                i+1, heavy_atoms[10*res->type+j-1], res->name, ' ', res->num,                      
  	  				  sc[j][0], sc[j][1], sc[j][2]);

}        
printf("TER\nEND\n");
} 
*********/
      
      for (j=0;j<3;j++) {
        vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];        
        for (k=0;k<3;k++) {
          if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
        }        
      }

      for (j=0;j<nsc;j++) {
        for (k=0;k<3;k++) {
          sc[j][k] = rot_stat_coords[pos+j][k];
        }
      }         

      superimpose2(vv,lsys,3,sc,nsc);

      for (j=0;j<nsc;j++) {
        sc[j][0] += x3;
        sc[j][1] += y3;
        sc[j][2] += z3;
      }

//printf("sc rebuilt, %d atoms\n", nsc);

      for (j=1;j<nsc;j++) {
        add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
      }        

      if (res->next) res = res->next;              

    } // i++, next res

    for (i=0;i<12;i++)
      free(sc[i]);
    for (i=0;i<3;i++) {
      free(lsys[i]);
      free(vv[i]);
    }
    free(sc); free(lsys); free(vv);
        
}


typedef struct _atom_list {
  atom_type *atom;
  struct _atom_list *next;
} atom_list;

int get_conflicts(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid)
{
  atom_list *llist;
  atom_type *atom, *atom2;
  int i, j, k, x, y, z;
  int ii, jj, kk, con, iter, maxcon, merged;
  real dx, dy, dz, dd;

    con = 0;
    atom = res->atoms;
    while (atom) {
      i = atom->gx;
      j = atom->gy;
      k = atom->gz;
      for (ii=i-2;ii<=i+2;ii++)
        for (jj=j-2;jj<=j+2;jj++)
          for (kk=k-2;kk<=k+2;kk++) {
            if (ii>=0 && ii<xgrid && jj>=0 && jj<ygrid && kk>=0 && kk<zgrid) {
              llist = grid[ii][jj][kk];
              while (llist) {
                atom2 = llist->atom;
                if (atom && atom2 && res && atom2->res) {
                  merged=0;
                  if (res==atom2->res) { // self-xvol
                    if (atom->flag & FLAG_SIDECHAIN && atom2->flag & FLAG_SIDECHAIN) merged=1;
                    if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='A' && atom2->name[0]=='C' && atom2->name[1]=='B') merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='B' && atom2->name[0]=='C' && atom2->name[1]=='A') merged=1;
                    if (res->name[0]=='P') {
                      if (atom->name[0]=='C' && atom->name[1]=='D' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                      if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]=='D') merged=1;
                    }
                    
                    if (!merged) {
//                      printf("merged: %s[%d] %s-%s %d %d\n", res->name,res->num,atom->name,atom2->name,atom->flag,atom2->flag);
                    }
                  } else
                  if (res->next==atom2->res || res==atom2->res->next) {
                    if (atom->name[0]=='C' && atom->name[1]==' ' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                    if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]==' ') merged=1;
                  }
                  if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1; // for now
                  if (atom->flag & FLAG_SCM || atom2->flag & FLAG_SCM) merged=1; // for now
                  if (!merged) {
                    dx = atom->x-atom2->x; 
                    dx*=dx;
                    dy = atom->y-atom2->y; 
                    dy*=dy;
                    dz = atom->z-atom2->z; 
                    dz*=dz;
                    dd = dx+dy+dz;
                    if (dd<1.6*1.6) {
                      con++;                        
                    }
                  }  
                }  
                llist = llist->next;
              }
            }
          }            
      atom = atom->next;
    }        
  
  return con;  
}



int display_conflicts(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid)
{
  atom_list *llist;
  atom_type *atom, *atom2;
  int i, j, k, x, y, z;
  int ii, jj, kk, con, iter, maxcon, merged;
  real dx, dy, dz, dd;

    con = 0;
    atom = res->atoms;
    while (atom) {
      i = atom->gx;
      j = atom->gy;
      k = atom->gz;
      for (ii=i-2;ii<=i+2;ii++)
        for (jj=j-2;jj<=j+2;jj++)
          for (kk=k-2;kk<=k+2;kk++) {
            if (ii>=0 && ii<xgrid && jj>=0 && jj<ygrid && kk>=0 && kk<zgrid) {
              llist = grid[ii][jj][kk];
              while (llist) {
                atom2 = llist->atom;
                if (atom && atom2 && res && atom2->res) {
                  merged=0;
                  if (res==atom2->res) { // self-xvol
                    if (atom->flag & FLAG_SIDECHAIN && atom2->flag & FLAG_SIDECHAIN) merged=1;
                    if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='A' && atom2->name[0]=='C' && atom2->name[1]=='B') merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='B' && atom2->name[0]=='C' && atom2->name[1]=='A') merged=1;
                    if (res->name[0]=='P') {
                      if (atom->name[0]=='C' && atom->name[1]=='D' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                      if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]=='D') merged=1;
                    }                    
                    if (!merged) {
//                      printf("merged: %s[%d] %s-%s %d %d\n", res->name,res->num,atom->name,atom2->name,atom->flag,atom2->flag);
                    }
                  } else
                  if (res->next==atom2->res || res==atom2->res->next) {
                    if (atom->name[0]=='C' && atom->name[1]==' ' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                    if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]==' ') merged=1;
                  }
                  if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1; // for now
                  if (atom->flag & FLAG_SCM || atom2->flag & FLAG_SCM) merged=1; // for now
                  if (!merged) {
                    dx = atom->x-atom2->x; 
                    dx*=dx;
                    dy = atom->y-atom2->y; 
                    dy*=dy;
                    dz = atom->z-atom2->z; 
                    dz*=dz;
                    dd = dx+dy+dz;
                    if (dd<1.6*1.6) {
                      printf("STERIC CONFLICT: %s[%d]%s-%s[%d]%s\n", atom->res->name,atom->res->num,atom->name,atom2->res->name,atom2->res->num,atom2->name);
                      con++;                        
                    }
                  }  
                }  
                llist = llist->next;
              }
            }
          }            
      atom = atom->next;
    }        
  
  return con;  
}



void optimize_exvol(void)
{
  real min[3], max[3];
  res_type *res, *worst;
  atom_type *atom, *atom2;
  int xgrid, ygrid, zgrid;
  atom_list ****grid, *llist, *alist;
  int i, j, k, l, m, x, y, z;
  int ii, jj, kk, con, iter, maxcon, totcon;
  int cpos, bestpos, pos, con0;
  real v1[3], v2a[3], v2b[3], v2[3], v3[3];
  int nsc, nca;
  real cax, cay, caz;
  real **lsys, **vv, **sc;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  
    min[0]=1e5;
    min[1]=1e5;
    min[2]=1e5;
    max[0]=-1e5;
    max[1]=-1e5;
    max[2]=-1e5;
    
    lsys = (real**)calloc(sizeof(real*)*3,1);
    vv = (real**)calloc(sizeof(real*)*3,1);
    sc = (real**)calloc(sizeof(real*)*12,1);
    for (i=0;i<12;i++)
      sc[i] = (real*)calloc(sizeof(real)*3,1);
    for (i=0;i<3;i++) {
      lsys[i] = (real*)calloc(sizeof(real)*3,1);
      vv[i] = (real*)calloc(sizeof(real)*3,1);
    }

    res = chain->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (atom->x<min[0]) min[0]=atom->x;
        if (atom->y<min[1]) min[1]=atom->y;
        if (atom->z<min[2]) min[2]=atom->z;
        if (atom->x>max[0]) max[0]=atom->x;
        if (atom->y>max[1]) max[1]=atom->y;
        if (atom->z>max[2]) max[2]=atom->z;
        atom = atom->next;
      }
      res = res->next;
    }
    
    xgrid = (max[0]-min[0])/GRID_RES;
    ygrid = (max[1]-min[1])/GRID_RES;
    zgrid = (max[2]-min[2])/GRID_RES;
    
//    if (_VERBOSE) printf("Allocating grid (%d %d %d)...\n", xgrid, ygrid, zgrid);

    grid = (atom_list****)calloc(sizeof(atom_list***)*(xgrid+1),1);
    for (i=0;i<xgrid+1;i++) {
      grid[i] = (atom_list***)calloc(sizeof(atom_list**)*(ygrid+1),1);
      for (j=0;j<ygrid+1;j++) {
        grid[i][j] = (atom_list**)calloc(sizeof(atom_list*)*(zgrid+1),1);
      }    
    }

    res = chain->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        x = xgrid*(atom->x-min[0])/(max[0]-min[0]);
        y = ygrid*(atom->y-min[1])/(max[1]-min[1]);
        z = zgrid*(atom->z-min[2])/(max[2]-min[2]);
        alist = (atom_list*)calloc(sizeof(atom_list),1);
        alist->atom = atom;
        atom->gx = x;
        atom->gy = y;
        atom->gz = z;
        if (grid[x][y][z]!=NULL) {
          llist = grid[x][y][z];
          while (llist->next) llist=llist->next;
          llist->next = alist;
        } else {
          grid[x][y][z]=alist;
        }  
        atom = atom->next;
      }
      res = res->next;
    }    

    if (_VERBOSE) printf("Finding excluded volume conflicts...\n", xgrid, ygrid, zgrid);
    
  
iter = 0;

do {
//printf("ITER: %d\n", iter);

      maxcon = 0;
      totcon=0;

      res = chain->residua;
      while (res) {
        con = get_conflicts(res, grid, xgrid, ygrid, zgrid);
        if (con>0) {
          totcon+=con;
          if (con>maxcon) {
            maxcon = con;
            worst = res;
          }
        }
        res = res->next;
      }      

      if (_VERBOSE && iter==0) {
        printf("Total number of conflicts: %d\n", totcon);         
      }
      
      if (totcon==0) break;
      
      if (_VERBOSE && iter==0) {
        printf("Maximum number of conflicts: %s[%d] : %d\n", worst->name, worst->num, maxcon);         
      }
      
      totcon=0;

      if (maxcon>0) {
            
// try to fix...


    res = chain->residua;
    for (i=0;i<chain_length;i++) {
      if (!strncmp(res->name,"GLY",3)) {
        if (res->next) res = res->next;
        continue;
      }

      nsc = nheavy[res->type]+1;

    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];
    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];
    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];
    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];

      v1[0] = x4-x2;
      v1[1] = y4-y2;
      v1[2] = z4-z2;
      
      v2a[0] = x4-x3;
      v2a[1] = y4-y3;
      v2a[2] = z4-z3;

      v2b[0] = x3-x2;
      v2b[1] = y3-y2;
      v2b[2] = z3-z2;
      
      cross(v2a, v2b, v2);
      cross(v1, v2, v3);
      
      norm(v1);
      norm(v2);
      norm(v3);

      con = get_conflicts(res, grid, xgrid, ygrid, zgrid);

      if (con>0) {

//printf("fixing %s[%d]: %d\n", res->name, res->num, con);

        bestpos=0;
        con0 = 100;      
        for (l=0;l<10;l++) { // check two closest conformations
          cpos = SORTED_ROTAMERS[i][l][1];
          for (m=0;m<rot_stat_idx[cpos][4];m++) {
            for (j=0;j<3;j++) {
              vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];        
              for (k=0;k<3;k++) {
                if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
              }        
            }
            pos = rot_stat_idx[cpos][5]+nsc*m;
            for (j=0;j<nsc;j++) {
              for (k=0;k<3;k++) {
                sc[j][k] = rot_stat_coords[pos+j][k];
              }
            }         
            superimpose2(vv,lsys,3,sc,nsc);      
            for (j=0;j<nsc;j++) {
              sc[j][0] += x3;
              sc[j][1] += y3;
              sc[j][2] += z3;
            }
            for (j=1;j<nsc;j++) {
              add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
            }        
            con = get_conflicts(res, grid, xgrid, ygrid, zgrid);
//printf("test: %d\n", con);

            if (con<con0) {
              con0 = con;
              bestpos = pos;
            }
            if (con==0) break;
          }
          if (con==0) break;
        }             

//printf("best: %d (%d)\n", con0, bestpos);        
totcon += con0;

        for (j=0;j<3;j++) {
          vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];        
          for (k=0;k<3;k++) {
            if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
          }        
        }
        pos = bestpos;
        for (j=0;j<nsc;j++) {
          for (k=0;k<3;k++) {
            sc[j][k] = rot_stat_coords[pos+j][k];
          }
        }         
        superimpose2(vv,lsys,3,sc,nsc);      
        for (j=0;j<nsc;j++) {
          sc[j][0] += x3;
          sc[j][1] += y3;
          sc[j][2] += z3;
        }
        for (j=1;j<nsc;j++) {
          add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
        }                  
//printf("replaced\n");        
      }

      
      
      res=res->next;

    } // i

//printf("finally: total number of violations: %d\n", totcon);
}

iter++;    

} while (iter<_XVOL_ITER);
    

    if (_VERBOSE) {
      if (totcon>0) 
        printf("WARNING: %d steric conflict(s) are still there.\n", totcon);
      else
        printf("All steric conflicts removed.\n");
    }    
    
    for (i=0;i<12;i++)
      free(sc[i]);
    for (i=0;i<3;i++) {
      free(lsys[i]);
      free(vv[i]);
    }
    free(sc); free(lsys); free(vv);
        

}


#endif

int main(int argc, char **argv)
{
  int i, j, next;
  char buf[100];
  char *name=NULL, *ini_name=NULL;
  char out_name[1000];
  real f;
  mol_type *mol;
  struct timeb time0, time1;
  
    for (i=1; i<argc; i++) {
      if (argv[i][0]=='-') {                    
        next=0;
        for (j=1; j<(int)strlen(argv[i]); j++) {
          switch(argv[i][j]) {
            case 'v': _VERBOSE=1; break;
            case 'c': _CA_OPTIMIZE=0; break;
            case 'r': _CA_RANDOM=1; break;
            case 't': _CA_TRAJECTORY=1; break;
            case 'n': _CENTER_CHAIN=1; break;
            case 'b': _REBUILD_BB=0; break;
            case 's': _REBUILD_SC=0; break;
            case 'i': ini_name = argv[++i]; next=1; break;
            case 'g': _PDB_SG=1; break;
            case 'x': _TIME_SEED=1; break;
            case 'o': _XVOLUME=0; break;
            case 'h': _REBUILD_H=0; break;
            case 'p': _CISPRO=1; break;
            case 'u': 
              if (sscanf(argv[++i],"%lf",&f)==1) {
                _CA_START_DIST = f;
              }
              next=1;
            break;            
            default: {
              printf("Unknown option: %c\n", argv[i][j]);
              return -1;
            }
          }     
          if (next) break;
        }  
      } else {
        if (!name) name=argv[i];
      }
    }

    if (!name) {
      printf("PULCHRA PowerfUL Chain Restoration Algorithm version %4.2f\n", PULCHRA_VERSION);
      printf("Usage: %s [options] <pdb_file>\n", argv[0]);
      printf("The program default input is a PDB file.\n");
      printf("Output file <pdb_file.pdb.rebuilt> will be created as a result.\n");
      printf("Valid options are:\n\n");
      printf("  -v : verbose output (default: off)\n");      
      printf("  -n : center chain (default: off)\n");
      printf("  -x : time-seed random number generator (default: off)\n");
      printf("  -g : use PDBSG as an input format (CA=C-alpha, SC or CM=side chain c.m.)\n\n");

      printf("  -c : skip C-alpha optimizing (default: on)\n");
      printf("  -p : detect cis-prolins (default: off)\n");
      printf("  -r : start from a random chain (default: off)\n");
      printf("  -i pdbfile : read the initial C-alpha coordinates from a PDB file\n");      
      printf("  -t : save chain optimization trajectory to file <pdb_file.pdb.trajectory>\n");
      printf("  -u value : maximum shift from the restraint coordinates (default: 0.5A)\n\n");

#ifdef COMPILE_BB      
      printf("  -b : skip backbone (and side chains) reconstruction (default: on)\n");
#endif      
#ifdef COMPILE_ROT
      printf("  -s : skip side chains reconstruction (default: on)\n");
      printf("  -o : don't attempt to fix excluded volume conflicts (default: on)\n");
#endif      
      return -1;
    }

    for (i=0; i<255; i++) /* prepare hash table*/
      AA_NUMS[i] = 20; /* dummy aa code */
    for (i=0; i<20; i++)
      AA_NUMS[(int)SHORT_AA_NAMES[i]] = i;

    setbuf(stdout,0);

    if (_TIME_SEED) srand(time(NULL)); else srand(1234);

    if (_VERBOSE) printf("PULCHRA Protein Chain Restoration Algorithm version %4.2f\n", PULCHRA_VERSION);

    ftime(&time0);

    chain = new_mol();

    if (read_pdb_file(name,chain,"chain")==FILE_NOT_FOUND) {
      if (_VERBOSE) printf("Can't read the input file!\n");
      return -1;
    }  


    if (_VERBOSE) printf("%d residua read.\n", chain->nres);

    chain_length = chain->nres;
    
    if (_CA_OPTIMIZE) {
      sprintf(out_name,"%s.tra",name);
      ca_optimize(out_name, ini_name);
    }  
    
#ifdef COMPILE_BB    
    if (_REBUILD_BB) {
      rebuild_backbone();
    }  
#endif    

#ifdef COMPILE_ROT
    if (_REBUILD_BB && _REBUILD_SC) {
      rebuild_sidechains();
      if (_XVOLUME)
        optimize_exvol();
    }  
#endif

    if (_CENTER_CHAIN) {
      center_chain(chain);
    }
  
          
    sprintf(out_name,"%s.rebuilt",name);
    if (_VERBOSE) printf("Writing output file %s...\n", out_name);
    write_pdb(out_name, chain);

//    sprintf(out_name,"%s.rebuilt.sg",name);
//    write_pdb_sg(out_name, chain);

    ftime(&time1);

    if (_VERBOSE) printf("Done. Reconstruction finished in %.3f s.\n", (real)0.001*(1000.*(time1.time-time0.time)+(time1.millitm-time0.millitm)));
    
    return 0;
}


