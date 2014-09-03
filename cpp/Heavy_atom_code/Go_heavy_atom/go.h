#ifndef _GO_
#define _GO_

#include "pdbClasses2.h"
#include "atom_param.h"

double d2( double *, double * ) ;

void output_nat_cont( PDBchain &, listPDBatom &, listPDBatom &, const int &min_sep,
		      const double &, ifstream &, int *&, int *&,
		      double *&, int & ) ;

void shift_to_center( double **, const int &, const double & ) ;

double output_fraction( PDBchain &, double **, int *, int *, double *, const int & ) ;


#endif /*_GO_*/
